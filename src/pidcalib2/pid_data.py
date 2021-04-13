###############################################################################
# (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################

import json
import os
import pickle
import re
from pathlib import Path
from typing import Dict, List

import boost_histogram as bh
import pandas as pd
import uproot
import uproot3
from logzero import logger as log
from tqdm import tqdm

from . import utils

# Dict of mothers for each particle type
mothers = {"Pi": ["DSt"], "K": ["DSt"], "Mu": ["Jpsi"]}

# Dict of charges for each particle type
charges = {"Pi": ["P", "M"], "K": ["P", "M"], "Mu": ["P", "M"], "P": ["", "bar"]}

run1_samples = [
    "13b",
    "15",
    "17",
    "20",
    "20_MCTuneV2",
    "20_MCTunev3",
    "20r1",
    "20r1_MCTuneV2",
    "21",
    "21_MCTuneV4",
    "21r1",
    "21r1_MCTuneV4",
    "22",
    "23",
    "23Val",
    "23_MCTuneV1",
    "26",
    "5TeV",
]


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
        "Brunel_P": f"{prefix}_Brunel_P",
        "ETA": f"{prefix}_ETA",
        "Brunel_ETA": f"{prefix}_Brunel_ETA",
        "nTracks": "nTracks",
        "nTracks_Brunel": "nTracks_Brunel",
        "nSPDhits": "nSPDhits",
        "nSPDhits_Brunel": "nSPDhits_Brunel",
        "sWeight": f"{prefix}_sWeight",
        "TRCHI2NDOF": f"{prefix}_TRCHI2NDOF",
    }
    return branch_names


def is_run1(sample: str) -> bool:
    return sample in run1_samples


def get_relevant_branch_names(
    prefix: str, pid_cuts: List[str], bin_vars: List[str], cuts: List[str] = None
) -> Dict[str, str]:
    """Return a list of branch names relevant to the cuts and binning vars.

    Args:
        prefix: A prefix for the branch names, i.e., "probe".
        pid_cuts: Simplified user-level cut list, e.g., ["DLLK < 4"].
        bin_vars: Variables used for the binning.
        cuts: Arbitrary cut list, e.g., ["Dst_IPCHI2 < 10.0"].
    """
    branch_names = create_branch_names(prefix)

    # Remove sWeight if not a calib sample
    if prefix != "probe":
        del branch_names["sWeight"]

    # Create a list of vars in the PID cuts
    pid_cuts_vars = []
    whitespace = re.compile(r"\s+")
    for pid_cut in pid_cuts:
        pid_cut = re.sub(whitespace, "", pid_cut)
        pid_cut_var, _, _ = re.split(r"(<|>|==|!=)", pid_cut)
        pid_cuts_vars.append(pid_cut_var)

    # Remove all vars that are not used for binning or PID cuts
    for branch in tuple(branch_names):
        if branch not in [*pid_cuts_vars, *bin_vars, "sWeight"]:
            del branch_names[branch]

    # Add vars in the arbitrary cuts
    if cuts:
        for cut in cuts:
            cut = re.sub(whitespace, "", cut)
            cut_var, _, _ = re.split(r"(<|>|==|!=)", cut)
            branch_names[cut_var] = cut_var

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


def get_file_list(
    sample: str,
    magnet: str,
    particle: str,
    samples_file: str,
    max_files: int = None,
) -> List[str]:
    """Return a list of calibration files.

    Args:
        sample: Data sample name (Turbo18, etc.)
        magnet: Magnet polarity (up, down)
        particle: Particle type (K, pi, etc.)
        samples_file: File in which to look up the file list.
        max_files: Optional. The maximum number of files to get. Defaults to
            None (= all files).
    """
    magnet = "Mag" + magnet.capitalize()
    if samples_file is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        samples_file = str(Path(current_dir, "data/samples.json"))

    sample_name = "-".join([sample, magnet, particle])

    log.debug(f"Reading file lists from '{samples_file}'")
    with open(samples_file) as f:
        samples = json.load(f)

    try:
        sample = samples[sample_name]
    except KeyError:
        log.error(f"Sample '{sample_name}' not found in {samples_file}")
        raise

    file_list = []
    if isinstance(sample, dict):
        link = sample["link"]
        try:
            file_list = samples[link]
        except KeyError:
            log.error(f"Linked sample '{link}' not found in {samples_file}")
            raise
    else:
        file_list = samples[sample_name]

    if max_files:
        file_list = file_list[:max_files]
    return file_list


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
    paths: List[str], tree_paths: List[str], branches: Dict[str, str]
) -> pd.DataFrame:
    """Read ROOT files via XRootD, extract branches, and save to a Pandas DF.

    DataFrame columns are not named after the branches in the ROOT trees, but
    rather by their associated 'simple user-level' analogues. E.g., 'DLLK'
    instead of 'probe_PIDK'.

    Args:
        paths: Paths to ROOT files; either file system paths or URLs, e.g.
            ["root:///eos/lhcb/file.root"].
        tree_paths: Internal ROOT file paths to the relevant trees
        branches: Dict of the branches {simple_name: branch_name} to include
           in the DataFrame.

    Returns:
        Pandas DataFrame with the requested branches from a decay tree of a
        single particle.
    """
    df_tot = pd.DataFrame()
    for path in tqdm(paths, leave=False, desc="Reading files"):
        for tree_path in tree_paths:
            df = root_to_dataframe(path, tree_path, list(branches.values()))
            df_tot = df_tot.append(df)

    # Rename colums of the dataset from branch names to simple user-level
    # names, e.g., probe_PIDK -> DLLK.
    inverse_branch_dict = {val: key for key, val in branches.items()}
    df_tot = df_tot.rename(columns=inverse_branch_dict)  # type: ignore

    log.info(f"Read {len(paths)} files with a total of {len(df_tot.index)} events")
    return df_tot


def get_tree_paths(particle: str, sample: str) -> List[str]:
    """Return a list of internal ROOT paths to relevant trees in the files

    Args:
        particle: Particle type (K, pi, etc.)
        sample: Data sample name (Turbo18, etc.)
    """
    tree_paths = []
    if is_run1(sample):
        # Run 1 files have a simple structure with a single tree
        tree_paths.append("DecayTree")
    else:
        # Run 2 ROOT file structure with multiple trees
        for mother in mothers[particle]:
            for charge in charges[particle]:
                tree_paths.append(
                    f"{mother}_{particle.capitalize()}{charge}Tuple/DecayTree"
                )

    return tree_paths


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


def get_calib_hists(
    hist_dir: str,
    sample: str,
    magnet: str,
    ref_pars: Dict[str, List[str]],
    bin_vars: Dict[str, str],
) -> Dict[str, Dict[str, bh.Histogram]]:
    """Get calibration efficiency histograms from all necessary files.

    Args:
        hist_dir: Directory where to look for the required files.
        sample: Data sample name (Turbo18, etc.).
        magnet: Magnet polarity (up, down).
        ref_pars: Reference particle prefixes with a particle type and PID cut.
        bin_vars: Binning variables ({standard name: reference sample branch name}).

    Returns:
        Dictionary with an efficiency histogram for each reference particle.
        The reference particle prefixes are the dictionary keys.
    """
    hists: Dict[str, Dict[str, bh.Histogram]] = {}
    for ref_par in ref_pars:
        particle = ref_pars[ref_par][0]

        pid_cut = ref_pars[ref_par][1]
        whitespace = re.compile(r"\s+")
        pid_cut = re.sub(whitespace, "", pid_cut)

        bin_str = ""
        for bin_var in bin_vars:
            bin_str += f"_{bin_var}"
        calib_name = Path(
            hist_dir,
            utils.create_hist_filename(
                sample, magnet, particle, pid_cut, list(bin_vars)
            ),
        )

        log.debug(f"Loading efficiency histograms from '{calib_name}'")

        hists[ref_par] = {}
        with open(calib_name, "rb") as f:
            hists[ref_par]["eff"] = pickle.load(f)
            hists[ref_par]["passing"] = pickle.load(f)
            hists[ref_par]["total"] = pickle.load(f)
            hists[ref_par]["passing_sumw2"] = pickle.load(f)
            hists[ref_par]["total_sumw2"] = pickle.load(f)
    return hists


def save_dataframe_as_root(
    df: pd.DataFrame, name: str, filename: str, columns: List[str] = None
):
    """Save a DataFrame as a TTree in a ROOT file.

    NaN entries are changed to -999 because ROOT TTrees don't support NaNs.

    Args:
        df: DataFrame to be saved.
        name: Name of the new TTree.
        filename: Name of the file to which to save the TTree.
        columns: Optional. Names of the columns which are to be saved. If
            'None', all the columns will be saved.
    """
    df_wo_nan = df.fillna(-999)
    if columns is None:
        columns = list(df_wo_nan.keys())
    branches_w_types = {branch: df_wo_nan[branch].dtype for branch in columns}
    with uproot3.recreate(filename) as f:
        log.debug(f"Creating a TTree with the following branches: {branches_w_types}")
        f[name] = uproot3.newtree(branches_w_types)
        branch_dict = {branch: df_wo_nan[branch] for branch in branches_w_types}
        f[name].extend(branch_dict)
    log.info(f"Efficiency tree saved to {filename}")
