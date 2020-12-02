import re
from typing import Dict, List

from tqdm import tqdm
import uproot4
import pandas as pd
from logzero import logger as log
from XRootD import client as xrdclient

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


# Make a list of variables to put into calib DataFrame, based on user-supplied
# PID and binning variables
# TODO: Write doc string
# TODO: Add type hints
def get_relevant_branch_names(name: str, pid_cuts: str, bin_vars: str) -> List[str]:
    """Return a list of variables to """
    branch_names = create_branch_names(name)
    vars = []
    # Add sWeight if calib sample
    if name == "probe":
        vars.append(branch_names["sw"])
    for bin_var in bin_vars:
        vars.append(branch_names[bin_var])
    for pid_cut in pid_cuts:
        for branch_name in branch_names:
            if branch_name in pid_cut and branch_names[branch_name] not in vars:
                vars.append(branch_names[branch_name])
    return vars


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

    log.info(f"Read {len(paths)} files")
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
    branch_names = create_branch_names(prefix)
    whitespace = re.compile(r"\s+")
    for pid_cut in pid_cuts:
        # Remove whitespace in the cut string
        pid_cut = re.sub(whitespace, "", pid_cut)
        pid_cut_var, _ = re.split(r"<|>", pid_cut)
        if pid_cut_var in branch_names:
            branch_cuts.append(pid_cut.replace(pid_cut_var, branch_names[pid_cut_var]))
        else:
            branch_cuts.append(pid_cut)
    return branch_cuts
