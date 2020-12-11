import re
from typing import Dict, List, Tuple

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

    assert magnet == "up" or magnet == "down"
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
        pid_cut_var, _, _ = re.split(r"(<|>)", pid_cut)
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

    Returns:
        List[str]: [description]
    """
    branch_names = []

    # TODO: Review these hardcoded branch names (maybe consolidate them
    # somewhere)
    if "nTracks" in bin_vars:
        branch_names.append(bin_vars["nTracks"])
    if "nSPDHits" in bin_vars:
        branch_names.append(bin_vars["nTracks"])

    for ref_par_name in ref_pars:
        for bin_var, bin_var_branch in bin_vars.items():
            if bin_var != "nTracks" and bin_var != "nSPDHits":
                branch_names.append(f"{ref_par_name}_{bin_var_branch}")
    return branch_names


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


def root_to_dataframe(
    path: str, tree_name: str, branches: Dict[str, str]
) -> pd.DataFrame:
    """Return DataFrame with requested branches from tree in ROOT file.

    Args:
        path: Path to the ROOT file; either file system path or URL, e.g.
            root:///eos/lhcb/file.root.
        tree_name: Name of a tree inside the ROOT file.
        branches: Branches to put in the DataFrame.
    """
    tree = uproot4.open(path)[tree_name]
    df = tree.arrays(branches.values(), library="pd")  # type: ignore
    # Rename colums of the dataset from branch names to simple user-level
    # names, e.g., probe_PIDK -> DLLK.
    inverse_branch_dict = {val: key for key, val in branches.items()}
    df = df.rename(columns=inverse_branch_dict)
    return df


def calib_root_to_dataframe(
    paths: List[str], particle: str, branches: Dict[str, str]
) -> pd.DataFrame:
    """Read ROOT files via XRootD, extract branches, and save to a Pandas DF.

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
                df = root_to_dataframe(path, tree_path, branches)
                df_tot = df_tot.append(df)

    log.info(f"Read {len(paths)} files with a total of {len(df_tot.index)} events")
    return df_tot


# TODO See if the prefix should be removed from the arguments
def translate_pid_cuts_to_branch_cuts(prefix: str, pid_cuts: List[str]) -> List[str]:
    """Return a list cuts that can be applied to the calib. DataFrame.

    The PID cuts are usually specified in a simplified notation such as 'DLLK
    < 4'. However, the calibration sample doesn't contain a 'DLLK' column,
    but rather 'probe_PIDK'. This function translates cuts in the simplified
    notation to the actual variable names in the calibration datasets.

    Args:
        prefix: A string to be prepended to each branch name, except nTracks.
        pid_cuts: Simplified user-level cut list, e.g., ["DLLK < 4"].
        bin_vars: Variables used for the binning (will be excluded from the
            translation).
    """
    branch_cuts = []
    for pid_cut in pid_cuts:
        pid_cut_var, pid_cut_string = pid_cut_to_branch_name_and_cut(prefix, pid_cut)
        branch_cuts.append(pid_cut_var + pid_cut_string)
    return branch_cuts


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
        axis_list.append(bh.axis.Variable(binning.binnings[particle][bin_var]))
        vals = df[bin_var].values
        vals_list.append(vals)

    # Create boost-histogram with the desired axes, and fill with sWeight applied
    hist = bh.Histogram(*axis_list)
    hist.fill(*vals_list, weight=df["sw"])

    return hist


def pid_cut_to_branch_name_and_cut(prefix: str, pid_cut: str) -> Tuple[str, str]:
    """Translate a PID cut in the simplified notation to a branch name.

    Args:
        prefix: A prefix to be prepended to the branch name, e.g., "probe".
        pid_cut: PID cut in the simplified notation, e.g., "DLLK > 4".

    Returns:
        A tuple of (translated branch name, cut string), e.g., ("probe_PIDK", "<4").
    """
    branch_names = create_branch_names(prefix)
    # Remove whitespace in the cut string
    whitespace = re.compile(r"\s+")
    pid_cut = re.sub(whitespace, "", pid_cut)
    pid_cut_var, delimiter, cut_string = re.split(r"(<|>)", pid_cut)
    # The cut variable must be in the branch_names - we can't cut on something
    # that is not present in the DataFrame
    try:
        return branch_names[pid_cut_var], delimiter + cut_string
    except KeyError:
        log.error(
            (
                f"The PID cut variable {pid_cut_var} is not among "
                f"known branch names: {branch_names.keys()}"
            )
        )
        raise


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
        log.info(
            f"{n_passing}/{n_total} ({n_passing / n_total * 100:.2f}%) events passed the cut"
        )
        hist_passing = make_hist(df_passing, particle, bin_vars)
        log.debug(f"Created 'passing' histogram")

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
