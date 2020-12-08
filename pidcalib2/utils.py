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
) -> List[str]:
    """Return a list of branch names relevant to the PID cuts and binning vars.

    Args:
        prefix: A prefix for the branch names, i.e., "probe".
        pid_cuts: Simplified user-level cut list, e.g., ["DLLK < 4"].
        bin_vars: Variables used for the binning.
    """
    branch_names = create_branch_names(prefix)
    relevant_branch_names = []
    # Add sWeight if calib sample
    if prefix == "probe":
        relevant_branch_names.append(branch_names["sw"])
    for bin_var in bin_vars:
        relevant_branch_names.append(branch_names[bin_var])
    for pid_cut in pid_cuts:
        for branch_name in branch_names:
            if (
                branch_name in pid_cut
                and branch_names[branch_name] not in relevant_branch_names
            ):
                relevant_branch_names.append(branch_names[branch_name])
    return relevant_branch_names


def get_eos_paths(year: int, magnet: str) -> List[str]:
    """Get EOS paths of calibration files for a given year and magnet."""
    eos_url = "root://eoslhcb.cern.ch/"
    calib_subpath = pidcalib_sample_dir(year, magnet)
    calib_path = f"/eos/lhcb/grid/prod/lhcb/LHCb/{calib_subpath}/0000/"

    eos_fs = xrdclient.FileSystem(eos_url)
    status, listing = eos_fs.dirlist(calib_path)  # type: ignore
    if not status.ok:  # type: ignore
        raise Exception(status)

    return [f"{eos_url}{calib_path}{xrdpath.name}" for xrdpath in listing]


def extract_branches_to_dataframe(
    paths: List[str], particle: str, branches: List[str]
) -> pd.DataFrame:
    """Read ROOT files via XRootD, extract branches, and save to a Pandas DF.

    Args:
        paths: List of XRootD URLs to the files to be read.
        particle: Particle whose decay tree branches are to be extracted.
        branches: Names of the branches to include in the DataFrame.

    Returns:
        Pandas DataFrame with the requested branches from a decay tree of a
        single particle.
    """
    df_tot = pd.DataFrame()
    log.info(f"Reading files...")
    # TODO Read *all* files
    for path in tqdm(paths[:50], leave=False):
        # log.info(f"Loading {path}")
        for mother in mothers[particle]:
            for charge in charges[particle]:
                tree_path = f"{mother}_{particle.capitalize()}{charge}Tuple/DecayTree"
                tree = uproot4.open(path)[tree_path]
                # Convert to a Pandas DF with a subset of branches
                df = tree.arrays(branches, library="pd")  # type: ignore
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

    Returns:
        bh.Histogram: [description]
    """
    axis_list = []
    vals_list = []

    # The first branch should always be the sWeight
    sweights = df.iloc[:, [0]]
    binning_branches = list(df.columns)[1 : 1 + len(bin_vars)]

    # Loop over bin dimensions and define the axes
    for i, bin_var in enumerate(bin_vars):
        axis_list.append(bh.axis.Variable(binning.binnings[particle][bin_var]))
        vals = df[binning_branches[i]].values
        vals_list.append(vals)

    # Create boost-histogram with the desired axes, and fill with sWeight applied
    hist = bh.Histogram(*axis_list)
    hist.fill(*vals_list, weight=sweights)

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


def create_histograms(
    df_total: pd.DataFrame, particle: str, pid_cuts: List[str], bin_vars: List[str]
) -> Dict[str, bh.Histogram]:
    hists = {}

    hists["total"] = make_hist(df_total, particle, bin_vars)
    zero_bins = np.count_nonzero(hists["total"].view(flow=False) == 0)
    if zero_bins:
        log.warning(f"There are {zero_bins} empty bins in the total histogram!")
        print(hists["total"].view(flow=False))

    n_total = len(df_total.index)
    for i, pid_cut in enumerate(pid_cuts):
        log.debug(f"Processing '{pid_cuts[i]}' cut")
        df_passing = df_total.query(pid_cut)
        n_passing = len(df_passing.index)
        log.debug(
            f"{n_passing}/{n_total} ({n_passing / n_total * 100:.2f}%) events passed the cut"
        )
        hists[f"pass_{pid_cut}"] = make_hist(df_passing, particle, bin_vars)
        log.debug(f"Created 'pass_{pid_cut}' histogram")

        hists[f"eff_{pid_cut}"] = hists[f"pass_{pid_cut}"].copy()
        hists[f"eff_{pid_cut}"][...] = hists[f"pass_{pid_cut}"].view(
            flow=False
        ) / hists["total"].view(flow=False)
        log.debug(f"Created 'eff_{pid_cut}' histogram")

    return hists
