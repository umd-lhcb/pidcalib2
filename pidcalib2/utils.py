import pickle
import re
from typing import Dict, List

import boost_histogram as bh
import numpy as np
import pandas as pd
import uproot
from logzero import logger as log
from tqdm import tqdm
from XRootD import client as xrdclient

from . import binning

# Dict of mothers for each particle type
mothers = {"pi": ["DSt"], "K": ["DSt"], "mu": ["Jpsi"]}

# Dict of charges for each particle type
charges = {"pi": ["P", "M"], "K": ["P", "M"], "mu": ["P", "M"], "p": ["", "bar"]}


def create_branch_names(prefix: str) -> Dict[str, str]:
    """Return a dict of {var name: branch name in the calib. tuple}.

    Args:
        prefix: A string to be prepended to each particle-specific branch
            name.
    """
    branch_names = {
        "DLLK": f"{prefix}_PIDK",
        "DLLp": f"{prefix}_PIDp",
        "ProbNNpi": f"{prefix}_MC15TuneV1_ProbNNpi",
        "ProbNNk": f"{prefix}_MC15TuneV1_ProbNNk",
        "P": f"{prefix}_P",
        "ETA": f"{prefix}_ETA",
        "nTracks": "nTracks",
        "nTracks_Brunel": "nTracks_Brunel",
        "nSPDhits": "nSPDhits",
        "nSPDhits_Brunel": "nSPDhits_Brunel",
        "sw": f"{prefix}_sWeight",
        "TRCHI2NDOF": f"{prefix}_TRCHI2NDOF",
    }
    return branch_names


def pidcalib_sample_dir(year: int, magnet: str) -> str:
    """Return path to EOS dir with relevant PIDCalib samples.

    Args:
        year: Data-taking year
        magnet: Magnet polarity (up, down)
    """
    # TODO Make this a simple dict instead of function
    dirs = {
        2015: {
            "up": "Collision15/PIDCALIB.ROOT/00064787",
            "down": "Collision15/PIDCALIB.ROOT/00064785",
        },
        2016: {
            "up": "Collision16/PIDCALIB.ROOT/00064793",
            "down": "Collision16/PIDCALIB.ROOT/00064795",
        },
        2017: {
            "up": "Collision17/PIDCALIB.ROOT/00090825",
            "down": "Collision17/PIDCALIB.ROOT/00090823",
        },
        2018: {
            "up": "Collision18/PIDCALIB.ROOT/00082947",
            "down": "Collision18/PIDCALIB.ROOT/00082949",
        },
    }

    assert magnet in ("up", "down")
    assert year in dirs

    return dirs[year][magnet]


def get_relevant_branch_names(
    prefix: str, pid_cuts: List[str], bin_vars: List[str]
) -> Dict[str, str]:
    """Return a list of branch names relevant to the PID cuts and binning vars.

    Args:
        prefix: A prefix for the branch names, i.e., "probe".
        pid_cuts: Simplified user-level cut list, e.g., ["DLLK < 4"].
        bin_vars: Variables used for the binning.
    """
    branch_names = create_branch_names(prefix)

    # Remove sWeight if not a calib sample
    if prefix != "probe":
        del branch_names["sw"]

    # Create a list of vars in the PID cuts
    pid_cuts_vars = []
    whitespace = re.compile(r"\s+")
    for pid_cut in pid_cuts:
        pid_cut = re.sub(whitespace, "", pid_cut)
        pid_cut_var, _, _ = re.split(r"(<|>|==|!=)", pid_cut)
        pid_cuts_vars.append(pid_cut_var)

    # Remove all vars that are not used for binning or PID cuts
    for branch in tuple(branch_names):
        if branch not in [*pid_cuts_vars, *bin_vars, "sw"]:
            del branch_names[branch]

    return branch_names


def get_reference_branch_names(
    ref_pars: Dict[str, List[str]], bin_vars: Dict[str, str]
) -> List[str]:
    """Return a list of relevant branch names in the reference sample.

    Args:
        ref_pars: A dict of {particle branch prefix : [particle type, PID cut]}
    """
    branch_names = []

    for ref_par_name in ref_pars:
        for bin_var, bin_var_branch in bin_vars.items():
            branch_name = get_reference_branch_name(
                ref_par_name, bin_var, bin_var_branch
            )
            # Avoid duplicate entries
            if branch_name not in branch_names:
                branch_names.append(branch_name)
    return branch_names


def get_reference_branch_name(prefix: str, bin_var: str, bin_var_branch: str) -> str:
    """Return a full name of a binning branch in the reference data.

    Args:
        prefix: Branch prefix of the particle in the reference sample.
        bin_var: Variable used for the binning.
        bin_var_branch: Branch name of the variable used for binning.
    """
    # TODO: Review these hardcoded branch names (maybe consolidate them
    # somewhere). Maybe add some checks that the bin_var is known.
    global_branches = ("nTracks", "nTracks_Brunel", "nSPDhits", "nSPDhits_Brunel")
    if bin_var in global_branches:
        return bin_var_branch

    return f"{prefix}_{bin_var_branch}"


def get_eos_paths(year: int, magnet: str, max_files: int = None) -> List[str]:
    """Get EOS paths of calibration files for a given year and magnet."""
    eos_url = "root://eoslhcb.cern.ch/"
    calib_subpath = pidcalib_sample_dir(year, magnet)
    calib_path = f"/eos/lhcb/grid/prod/lhcb/LHCb/{calib_subpath}/0000/"

    eos_fs = xrdclient.FileSystem(eos_url)
    status, listing = eos_fs.dirlist(calib_path)  # type: ignore
    if not status.ok:  # type: ignore
        raise Exception(status)

    paths = [f"{eos_url}{calib_path}{xrdpath.name}" for xrdpath in listing]
    if max_files:
        paths = paths[:max_files]
    return paths


def root_to_dataframe(path: str, tree_name: str, branches: List[str]) -> pd.DataFrame:
    """Return DataFrame with requested branches from tree in ROOT file.

    Args:
        path: Path to the ROOT file; either file system path or URL, e.g.
            root:///eos/lhcb/file.root.
        tree_name: Name of a tree inside the ROOT file.
        branches: Branches to put in the DataFrame.
    """
    tree = uproot.open(path)[tree_name]
    df = tree.arrays(branches, library="pd")  # type: ignore
    return df


def calib_root_to_dataframe(
    paths: List[str], particle: str, branches: Dict[str, str]
) -> pd.DataFrame:
    """Read ROOT files via XRootD, extract branches, and save to a Pandas DF.

    DataFrame columns are not named after the branches in the ROOT trees, but
    rather by their associated 'simple user-level' analogues. E.g., 'DLLK'
    instead of 'probe_PIDK'.

    Args:
        paths: Paths to ROOT files; either file system paths or URLs, e.g.
            ["root:///eos/lhcb/file.root"].
        particle: Particle whose decay tree branches are to be extracted.
        branches: Dict of the branches {simple_name: branch_name} to include
           in the DataFrame.

    Returns:
        Pandas DataFrame with the requested branches from a decay tree of a
        single particle.
    """
    df_tot = pd.DataFrame()
    for path in tqdm(paths, leave=False, desc="Reading files"):
        for mother in mothers[particle]:
            for charge in charges[particle]:
                tree_path = f"{mother}_{particle.capitalize()}{charge}Tuple/DecayTree"
                df = root_to_dataframe(path, tree_path, list(branches.values()))
                df_tot = df_tot.append(df)

    # Rename colums of the dataset from branch names to simple user-level
    # names, e.g., probe_PIDK -> DLLK.
    inverse_branch_dict = {val: key for key, val in branches.items()}
    df_tot = df_tot.rename(columns=inverse_branch_dict)  # type: ignore

    log.info(f"Read {len(paths)} files with a total of {len(df_tot.index)} events")
    return df_tot


def make_hist(df: pd.DataFrame, particle: str, bin_vars: List[str]) -> bh.Histogram:
    """Create a histogram of sWeighted events with appropriate binning

    Args:
        df: DataFrame from which to histogram events
        particle: Particle type (K, pi, etc.)
        bin_vars: Binning variables in the user-convention, e.g., ["P", "ETA"]
    """
    axis_list = []
    vals_list = []

    # Loop over bin dimensions and define the axes
    for bin_var in bin_vars:
        axis_list.append(
            bh.axis.Variable(
                binning.binnings[particle][bin_var], metadata={"name": bin_var}
            )
        )
        vals = df[bin_var].values
        vals_list.append(vals)

    # Create boost-histogram with the desired axes, and fill with sWeight applied
    hist = bh.Histogram(*axis_list)
    hist.fill(*vals_list, weight=df["sw"])

    return hist


def create_eff_histograms(
    df_total: pd.DataFrame, particle: str, pid_cuts: List[str], bin_vars: List[str]
) -> Dict[str, bh.Histogram]:
    """Create efficiency histograms for all supplied PID cuts.

    Args:
        df_total: DataFrame with all events
        particle: Particle type (K, pi, etc.)
        pid_cuts: Simple user-level cut list, e.g., ["DLLK<4"].
        bin_vars: Variables used for binning.

    Returns:
        A dictionary with all the efficiency histograms, with the PID cuts as
        keys.
    """

    hists = {}
    hists["total"] = make_hist(df_total, particle, bin_vars)

    zero_bins = np.count_nonzero(hists["total"].view(flow=False) == 0)
    if zero_bins:
        log.warning(f"There are {zero_bins} empty bins in the total histogram!")
        log.warning(hists["total"].view(flow=False))

    n_total = len(df_total.index)
    for i, pid_cut in enumerate(pid_cuts):
        log.info(f"Processing '{pid_cuts[i]}' cut")
        df_passing = df_total.query(pid_cut)
        n_passing = len(df_passing.index)
        percent_passing = n_passing / n_total * 100
        log.info(
            f"{n_passing}/{n_total} ({percent_passing:.2f}%) events passed the cut"
        )
        hists[f"passing_{pid_cut}"] = make_hist(df_passing, particle, bin_vars)
        log.debug("Created 'passing' histogram")

        hists[f"eff_{pid_cut}"] = hists[f"passing_{pid_cut}"].copy()
        hists[f"eff_{pid_cut}"][...] = hists[f"passing_{pid_cut}"].view(
            flow=False
        ) / hists["total"].view(flow=False)
        log.debug(f"Created 'eff_{pid_cut}' histogram")

    return hists


def dataframe_from_local_file(path: str, branch_names: List[str]) -> pd.DataFrame:
    """Return a dataframe read from a local file (instead of EOS).

    Args:
        path: Path to the local file.
        branch_names: Columns to read from the DataFrame.
    """
    if path.endswith(".pkl"):
        df = pd.read_pickle(path)
    elif path.endswith(".csv"):
        df = pd.read_csv(path, index_col=0)
    else:
        log.error(
            (
                f"Local dataframe file '{path}' "
                f"has an unknown suffix (csv and pkl supported)"
            )
        )
        raise Exception("Only csv and pkl files supported")
    log.info(f"Read {path} with a total of {len(df.index)} events")

    try:
        df = df[branch_names]
    except KeyError:
        log.error("The requested branches are missing from the local file")
        raise

    return df


def log_config(config: dict) -> None:
    """Pretty-print a config/dict."""
    longest_key = len(max(config, key=len))
    log.info("=" * longest_key)
    for entry in config:
        if config[entry] is not None:
            log.info(f"{entry:{longest_key}}: {config[entry]}")
    log.info("=" * longest_key)


# Load calib hists from file
def get_calib_hists(
    hist_dir: str,
    year: int,
    magnet: str,
    ref_pars: Dict[str, List[str]],
    bin_vars: Dict[str, str],
) -> Dict[str, bh.Histogram]:
    """Get calibration efficiency histograms from all necessary files.

    Args:
        hist_dir: Directory where to look for the required files.
        year: Data-taking year.
        magnet: Magnet polarity (up, down).
        ref_pars: Reference particle prefixes with a particle type and PID cut.
        bin_vars: Binning variables ({standard name: reference sample branch name}).

    Returns:
        Dictionary with an efficiency histogram for each reference particle.
        The reference particle prefixes are the dictionary keys.
    """
    hists = {}
    for ref_par in ref_pars:
        particle = ref_pars[ref_par][0]

        pid_cut = ref_pars[ref_par][1]
        whitespace = re.compile(r"\s+")
        pid_cut = re.sub(whitespace, "", pid_cut)

        bin_str = ""
        for bin_var in bin_vars:
            bin_str += f"_{bin_var}"
        calib_name = f"{hist_dir}/effhist_{year}_{magnet}_{particle}_{pid_cut}{bin_str}"

        log.debug(f"Loading efficiency histogram from '{calib_name}'")

        with open(calib_name + ".pkl", "rb") as f:
            hists[ref_par] = pickle.load(f)
    return hists


def add_bin_indices(
    df: pd.DataFrame,
    prefixes: List[str],
    bin_vars: Dict[str, str],
    eff_hists: Dict[str, bh.Histogram],
) -> pd.DataFrame:
    """Return a DataFrame with added indices of bins for each event.

    The binnings of binning variables are taken from efficiency histograms.
    Each event falls into a certain bin in each binning variable. This bin's
    index is added to the DataFrame. The same procedure is repeated for each
    variable. Finally, a global index of the N-dimensional bin of the
    efficiency histogram where the event belongs is added. Multiple
    efficiency histograms can be specified since the total PID efficiency for
    the event can depend on multiple particles.

    Args:
        df: Input dataframe.
        prefixes: Branch prefixes of the particles in the reference sample.
        bin_vars: Variables used for binning.
        eff_hists: Efficiency histograms for each prefix/particle.
    """
    df_new = df.copy()
    for prefix in prefixes:
        axes = [
            get_reference_branch_name(
                prefix, axis.metadata["name"], bin_vars[axis.metadata["name"]]
            )
            for axis in eff_hists[prefix].axes
        ]
        for bin_var, branch_name in bin_vars.items():
            ref_branch_name = get_reference_branch_name(prefix, bin_var, branch_name)
            bins = []
            for axis in eff_hists[prefix].axes:
                if axis.metadata["name"] == bin_var:
                    bins = axis.edges
            df_new[f"{ref_branch_name}_index"] = pd.cut(
                df_new[ref_branch_name],
                bins,
                labels=False,
                include_lowest=True,
                right=False,
            )

        df_nan = df_new[df_new.isna().any(axis=1)]  # type: ignore
        df_new.dropna(inplace=True)
        index_names = [f"{axis}_index" for axis in axes]
        indices = np.ravel_multi_index(
            df_new[index_names].transpose().to_numpy().astype(int),  # type: ignore
            eff_hists[prefix].axes.size,
        )
        df_new[f"{prefix}_index"] = indices
        df_new = pd.concat([df_new, df_nan]).sort_index()
    log.debug("Bin indices assigned")
    return df_new  # type: ignore


def add_efficiencies(
    df: pd.DataFrame, prefixes: List[str], eff_hists: Dict[str, bh.Histogram]
) -> pd.DataFrame:
    """Return a DataFrame with added efficiencies for each event.

    Each particle correspondig to a prefix is assigned an efficiency. The
    total event efficiency is the product of individual efficiencies.

    Args:
        df: Input dataframe.
        prefixes: Branch prefixes of the particles in the reference sample.
        eff_hists: Efficiency histograms for each prefix/particle.
    """
    df_new = df.copy()
    df_nan = df_new[df_new.isna().any(axis=1)]
    df_new.dropna(inplace=True)
    df_new["eff"] = 1
    for prefix in prefixes:
        efficiency_table = eff_hists[prefix].view().flatten()
        df_new[f"{prefix}_eff"] = np.take(efficiency_table, df_new[f"{prefix}_index"])
        df_new["eff"] = df_new["eff"] * df_new[f"{prefix}_eff"]

    df_new = pd.concat([df_new, df_nan]).sort_index()
    log.debug("Particle efficiencies assigned")

    num_outside_range = len(df_nan.index)
    num_outside_range_frac = len(df_nan.index) / len(df_new.index)
    log.debug(
        f"Events out of range: {num_outside_range} ({num_outside_range_frac:.2%})"
    )
    return df_new
