import pickle
import re
from typing import Dict, List

import boost_histogram as bh
import numpy as np
import pandas as pd
import uproot4
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
        prefix: A string to be prepended to each branch name, except nTracks.
    """
    branch_names = {
        "DLLK": f"{prefix}_PIDK",
        "DLLp": f"{prefix}_PIDp",
        "ProbNNpi": f"{prefix}_MC15TuneV1_ProbNNpi",
        "ProbNNk": f"{prefix}_MC15TuneV1_ProbNNk",
        "P": f"{prefix}_P",
        "ETA": f"{prefix}_ETA",
        "nTracks": "nTracks_Brunel",
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
    if "nTracks" == bin_var:
        return bin_var_branch
    if "nSPDHits" == bin_var:
        return bin_var_branch
    else:
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
    tree = uproot4.open(path)[tree_name]
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

    hist_total = make_hist(df_total, particle, bin_vars)
    zero_bins = np.count_nonzero(hist_total.view(flow=False) == 0)
    if zero_bins:
        log.warning(f"There are {zero_bins} empty bins in the total histogram!")
        log.warning(hist_total.view(flow=False))

    eff_hists = {}
    n_total = len(df_total.index)
    for i, pid_cut in enumerate(pid_cuts):
        log.info(f"Processing '{pid_cuts[i]}' cut")
        df_passing = df_total.query(pid_cut)
        n_passing = len(df_passing.index)
        percent_passing = n_passing / n_total * 100
        log.info(
            f"{n_passing}/{n_total} ({percent_passing:.2f}%) events passed the cut"
        )
        hist_passing = make_hist(df_passing, particle, bin_vars)
        log.debug("Created 'passing' histogram")

        eff_hists[f"eff_{pid_cut}"] = hist_passing.copy()
        eff_hists[f"eff_{pid_cut}"][...] = hist_passing.view(
            flow=False
        ) / hist_total.view(flow=False)
        log.debug(f"Created 'eff_{pid_cut}' histogram")

    return eff_hists


def dataframe_from_local_file(path: str, branch_names: List[str]) -> pd.DataFrame:
    """Return a dataframe read from a local file (instead of EOS).

    Args:
        path: Path to the local file.
        branch_names: Columns to read from the DataFrame.
    """
    if path.endswith("pkl"):
        df = pd.read_pickle(path)
    elif path.endswith("csv"):
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
        year: Data-taking year
        magnet: Magnet polarity (up, down)
        ref_pars (Dict[str, List[str]]): [description]
        bin_vars (Dict[str, str]): [description]
        TODO: Finish the docstring

    Returns:
        Dict[str, bh.Histogram]: [description]
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


# Calculate per-event effs for a given sample
def get_per_event_effs(df_ref, ref_pars, bin_vars, hists):
    # Create new column to hold the eff values
    df_ref["eff"] = -1.0
    # Per-particle effs
    for ref_par in ref_pars:
        df_ref[f"{ref_par}_eff"] = -1.0

    log.info("Calculating per event efficiencies...")
    # Loop over events and calculate per-event efficiency as the product of
    # individual track efficiencies
    for index, row in tqdm(
        df_ref.iterrows(), total=len(df_ref.index), leave=False, desc="Events"
    ):
        # Loop over tracks
        isAcc = True
        vals = {}
        for ref_par in ref_pars:
            # print(f"Track : {p}")
            # Get branch names of binning variables
            event_vals = {}
            axis_range = {}
            # Loop over hist axes (works for any number of bin dims)
            for i, bin_var in enumerate(bin_vars):
                branch_name = get_reference_branch_name(
                    ref_par, bin_var, bin_vars[bin_var]
                )
                event_vals[i] = row[branch_name]
                # print(f"{event_vals[i]}")
                for axe in hists[ref_par].axes:
                    if axe.metadata is None:
                        log.error("No axe metadata found!")
                        raise Exception("No axe metadata found!")
                    if axe.metadata["name"] == bin_var:
                        axis_range[i] = axe.edges
                # Check if track falls within hist range for this axis
                if (
                    event_vals[i] >= axis_range[i][0]
                    and event_vals[i] < axis_range[i][-1]
                ):
                    pass
                else:
                    isAcc = False
            event_vals_list = []
            if isAcc:
                for i, bin_var in enumerate(bin_vars):
                    event_vals_list.append(event_vals[i])
                # Get global index of the bin the track falls in
                index_num = hists[ref_par].axes.index(*event_vals_list)
                # Get bin content (eff value)
                vals[ref_par] = hists[ref_par].view()[index_num]
                # log.debug(f"Eff = {vals[ref_par]}")

        # Determine efficiency for the event
        if isAcc:
            eff = 1.0
            for ref_par in vals:
                # Fill track efficiency branch
                df_ref.loc[index, f"{ref_par}_eff"] = vals[ref_par]
                eff *= vals[ref_par]
            # log.debug(f"Event eff: {eff}")
        else:
            eff = -1.0
        # Fill total event efficiency branch
        df_ref.loc[index, "eff"] = eff

    # Return df with additional efficiency column
    return df_ref


def get_per_event_effs2(
    df_ref: pd.DataFrame,
    prefixes: List[str],
    bin_vars: Dict[str, str],
    eff_hists: Dict[str, bh.Histogram],
) -> pd.DataFrame:
    """Return input DataFrame with added 'eff' column with PID efficiency.

    Args:
        df_ref: DataFrame for which to calculate per-event efficiencies.
        prefixes: Prefixes for binning vars for each particle involved in PID cuts.
        bin_vars: Variables used for the binning.
        eff_hists: Efficiency histograms for each particle/prefix.
    """
    log.info("Calculating per-event efficiencies...")
    tqdm.pandas(desc="Events", leave=False)
    df_ref["eff"] = df_ref.progress_apply(
        calc_event_efficiency, axis=1, args=(prefixes, bin_vars, eff_hists)
    )

    num_in_range = len(df_ref[df_ref["eff"] == -1].index)
    num_in_range_frac = len(df_ref[df_ref["eff"] == -1].index) / len(df_ref.index)
    log.debug(f"Events outside range: {num_in_range} ({num_in_range_frac:.2%})")
    return df_ref


def calc_event_efficiency(
    row: pd.Series,
    prefixes: List[str],
    bin_vars: Dict[str, str],
    eff_hists: Dict[str, bh.Histogram],
) -> float:
    """Return the PID efficiency of the event.

    This function is intended to be used by the Pandas' apply() function.
    The event efficiency is a product of individual particle efficiencies.
    The involved particles are defined by the 'prefixes' list. The particle
    efficiencies are looked up in the relevant histogram in 'eff_hists'.

    Args:
        row: A row with a single event from a DataFrame.
        prefixes: Prefixes for binning vars for each particle involved in PID cuts.
        bin_vars: Variables used for the binning.
        eff_hists: Efficiency histograms for each particle/prefix.
    """
    efficiency = 1
    for prefix in prefixes:
        # Get an *ordered* list of branch names. For boost_histogram's
        # index() to work, the row[axes] values must come in the same order
        # as the axes in the histogram.
        axes = [
            get_reference_branch_name(
                prefix, axis.metadata["name"], bin_vars[axis.metadata["name"]]
            )
            for axis in eff_hists[prefix].axes
        ]
        # Get global index of the bin the track falls in
        index_num = eff_hists[prefix].axes.index(*row[axes])

        # Event efficiency is calculated as a product of every involved
        # particle efficiency. IndexError indicates the event is outside of the
        # efficiency histogram's range.
        try:
            efficiency *= eff_hists[prefix].view()[index_num]
        except IndexError:
            return -1

    return efficiency
