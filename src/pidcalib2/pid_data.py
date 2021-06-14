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
from typing import Any, Dict, List

import boost_histogram as bh
import pandas as pd
import uproot
import uproot3
from logzero import logger as log

from . import utils

tuple_names = {
    "Pi": ["DSt_PiMTuple", "DSt_PiPTuple"],
    "K": ["DSt_KPTuple", "DSt_KMTuple"],
    "Mu": ["Jpsi_MuPTuple", "Jpsi_MuMTuple"],
    "K_DsPhi": ["DsPhi_KPTuple", "DsPhi_KMTuple"],
    "K_DSt3Pi": ["DSt3Pi_KPTuple", "DSt3Pi_KMTuple"],
    "Pi_DSt3Pi": ["DSt3Pi_PiMTuple", "DSt3Pi_PiPTuple"],
    "K_Phi": ["Phi_KMTuple", "Phi_KPTuple"],
    "K_Omega": ["OmegaL_KPTuple", "OmegaL_KMTuple"],
    "K_DD": ["OmegaDDD_KPTuple", "OmegaDDD_KMTuple"],
    "Pi_KS": ["KSLL_PiPTuple", "KSLL_PiMTuple"],
    "Pi_DD": ["KSDD_PiPTuple", "KSDD_PiMTuple"],
    "P": [
        "Lam0LL_PTuple",
        "Lam0LL_PbarTuple",
        "Lam0LL_HPT_PTuple",
        "Lam0LL_HPT_PbarTuple",
        "Lam0LL_VHPT_PTuple",
        "Lam0LL_VHPT_PbarTuple",
    ],
    "P_DD": [
        "Lam0DD_PTuple",
        "Lam0DD_PbarTuple",
        "Lam0DD_HPT_PTuple",
        "Lam0DD_HPT_PbarTuple",
        "Lam0DD_VHPT_PTuple",
        "Lam0DD_VHPT_PbarTuple",
    ],
    "P_Sigmac0": ["Sigmac0_PbarTuple", "Sigmac0_PTuple"],
    "P_Sigmacpp": ["Sigmacpp_PbarTuple", "Sigmacpp_PTuple"],
    "P_B_Jpsi": ["B_Jpsi_PTuple", "B_Jpsi_PbarTuple"],
    "e_B_Jpsi": ["B_Jpsi_EPTuple", "B_Jpsi_EMTuple"],
    "Mu_B_Jpsi": ["B_Jpsi_MuPTuple", "B_Jpsi_MuMTuple"],
    "e": ["Jpsi_EPTuple", "Jpsi_EMTuple"],
    "P_LcfB": ["LbLcMu_PTuple", "LbLcMu_PbarTuple"],
    "P_IncLc": ["Lc_PTuple", "Lc_PbarTuple"],
    "Mu_nopt": ["Jpsinopt_MuPTuple", "Jpsinopt_MuMTuple"],
}

simple_samples = [
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
    "Electron15",
    "Electron16",
    "Electron17",
    "Electron18",
]

aliases = {
    "DLLK": "probe_PIDK",
    "DLLp": "probe_PIDp",
    "DLLmu": "probe_PIDmu",
    "DLLe": "probe_PIDe",
    "DLLd": "probe_PIDd",
    "MC12TuneV2_ProbNNe": "probe_MC12TuneV2_ProbNNe",
    "MC12TuneV2_ProbNNghost": "probe_MC12TuneV2_ProbNNghost",
    "MC12TuneV2_ProbNNk": "probe_MC12TuneV2_ProbNNk",
    "MC12TuneV2_ProbNNp": "probe_MC12TuneV2_ProbNNp",
    "MC12TuneV2_ProbNNpi": "probe_MC12TuneV2_ProbNNpi",
    "MC12TuneV3_ProbNNe": "probe_MC12TuneV3_ProbNNe",
    "MC12TuneV3_ProbNNghost": "probe_MC12TuneV3_ProbNNghost",
    "MC12TuneV3_ProbNNk": "probe_MC12TuneV3_ProbNNk",
    "MC12TuneV3_ProbNNp": "probe_MC12TuneV3_ProbNNp",
    "MC12TuneV3_ProbNNpi": "probe_MC12TuneV3_ProbNNpi",
    "MC12TuneV4_ProbNNe": "probe_MC12TuneV4_ProbNNe",
    "MC12TuneV4_ProbNNghost": "probe_MC12TuneV4_ProbNNghost",
    "MC12TuneV4_ProbNNk": "probe_MC12TuneV4_ProbNNk",
    "MC12TuneV4_ProbNNp": "probe_MC12TuneV4_ProbNNp",
    "MC12TuneV4_ProbNNpi": "probe_MC12TuneV4_ProbNNpi",
    "MC15TuneV1_ProbNNe": "probe_MC15TuneV1_ProbNNe",
    "MC15TuneV1_ProbNNghost": "probe_MC15TuneV1_ProbNNghost",
    "MC15TuneV1_ProbNNk": "probe_MC15TuneV1_ProbNNk",
    "MC15TuneV1_ProbNNp": "probe_MC15TuneV1_ProbNNp",
    "MC15TuneV1_ProbNNpi": "probe_MC15TuneV1_ProbNNpi",
    "Brunel_MC12TuneV2_ProbNNe": "probe_Brunel_MC12TuneV2_ProbNNe",
    "Brunel_MC12TuneV2_ProbNNghost": "probe_Brunel_MC12TuneV2_ProbNNghost",
    "Brunel_MC12TuneV2_ProbNNk": "probe_Brunel_MC12TuneV2_ProbNNk",
    "Brunel_MC12TuneV2_ProbNNp": "probe_Brunel_MC12TuneV2_ProbNNp",
    "Brunel_MC12TuneV2_ProbNNpi": "probe_Brunel_MC12TuneV2_ProbNNpi",
    "Brunel_MC12TuneV3_ProbNNe": "probe_Brunel_MC12TuneV3_ProbNNe",
    "Brunel_MC12TuneV3_ProbNNghost": "probe_Brunel_MC12TuneV3_ProbNNghost",
    "Brunel_MC12TuneV3_ProbNNk": "probe_Brunel_MC12TuneV3_ProbNNk",
    "Brunel_MC12TuneV3_ProbNNp": "probe_Brunel_MC12TuneV3_ProbNNp",
    "Brunel_MC12TuneV3_ProbNNpi": "probe_Brunel_MC12TuneV3_ProbNNpi",
    "Brunel_MC12TuneV4_ProbNNe": "probe_Brunel_MC12TuneV4_ProbNNe",
    "Brunel_MC12TuneV4_ProbNNghost": "probe_Brunel_MC12TuneV4_ProbNNghost",
    "Brunel_MC12TuneV4_ProbNNk": "probe_Brunel_MC12TuneV4_ProbNNk",
    "Brunel_MC12TuneV4_ProbNNp": "probe_Brunel_MC12TuneV4_ProbNNp",
    "Brunel_MC12TuneV4_ProbNNpi": "probe_Brunel_MC12TuneV4_ProbNNpi",
    "Brunel_MC15TuneV1_ProbNNe": "probe_Brunel_MC15TuneV1_ProbNNe",
    "Brunel_MC15TuneV1_ProbNNghost": "probe_Brunel_MC15TuneV1_ProbNNghost",
    "Brunel_MC15TuneV1_ProbNNk": "probe_Brunel_MC15TuneV1_ProbNNk",
    "Brunel_MC15TuneV1_ProbNNp": "probe_Brunel_MC15TuneV1_ProbNNp",
    "Brunel_MC15TuneV1_ProbNNpi": "probe_Brunel_MC15TuneV1_ProbNNpi",
    "P": "probe_P",
    "Brunel_P": "probe_Brunel_P",
    "ETA": "probe_ETA",
    "Brunel_ETA": "probe_Brunel_ETA",
    "nTracks": "nTracks",
    "nTracks_Brunel": "nTracks_Brunel",
    "nSPDhits": "nSPDhits",
    "nSPDhits_Brunel": "nSPDhits_Brunel",
    "TRCHI2NDOF": "probe_TRCHI2NDOF",
    "Brunel_TRCHI2NDOF": "probe_Brunel_TRCHI2NDOF",
    "IsMuon": "probe_isMuon",
    "Brunel_IsMuon": "probe_Brunel_isMuon",
    "HasRich": "probe_hasRich",
    "Brunel_HasRich": "probe_Brunel_hasRich",
}


def is_simple(sample: str) -> bool:
    """Return whether a sample has a simple directory structure.

    All Run 1 samples and Run 2 Electron samples have a single DecayTree inside.
    Standard Run 2 samples have multiple directories, each of which has a
    DecayTree inside. This function checks if the sample is in the list of
    simple samples.

    Args:
        sample: Sample name, e.g., Turbo15.
    """
    return sample in simple_samples


def get_relevant_branch_names(
    pid_cuts: List[str], bin_vars: List[str], cuts: List[str] = None
) -> Dict[str, str]:
    """Return a list of branch names relevant to the cuts and binning vars.

    Args:
        pid_cuts: Simplified user-level cut list, e.g., ["DLLK < 4"].
        bin_vars: Variables used for the binning.
        cuts: Arbitrary cut list, e.g., ["Dst_IPCHI2 < 10.0"].
    """
    branch_names = {"sWeight": "probe_sWeight"}

    whitespace = re.compile(r"\s+")
    for pid_cut in pid_cuts:
        pid_cut = re.sub(whitespace, "", pid_cut)
        pid_cut_vars = utils.extract_variable_names(pid_cut)

        for pid_cut_var in pid_cut_vars:
            if pid_cut_var not in aliases:
                log.warning(
                    (
                        f"PID cut variable '{pid_cut_var}' is not a known alias, "
                        "using raw variable"
                    )
                )
                branch_names[pid_cut_var] = pid_cut_var
            else:
                branch_names[pid_cut_var] = aliases[pid_cut_var]

    for bin_var in bin_vars:
        if bin_var not in aliases:
            log.warning(
                f"'Binning variable {bin_var}' is not a known alias, using raw variable"
            )
            branch_names[bin_var] = bin_var
        else:
            branch_names[bin_var] = aliases[bin_var]

    # Add vars in the arbitrary cuts
    if cuts:
        for cut in cuts:
            cut = re.sub(whitespace, "", cut)
            cut_vars = utils.extract_variable_names(cut)
            for cut_var in cut_vars:
                if cut_var not in aliases:
                    log.warning(
                        (
                            f"Cut variable '{cut_var}' is not a known alias, "
                            "using raw variable"
                        )
                    )
                    branch_names[cut_var] = cut_var
                else:
                    branch_names[cut_var] = aliases[cut_var]

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


def get_calibration_samples(samples_file: str = None):
    """Return a dictionary of all files for all calibration samples.

    Args:
        samples_file: JSON file with the calibration file lists.
    """
    if samples_file is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        samples_file = str(Path(current_dir, "data/samples.json"))

    log.debug(f"Reading file lists from '{samples_file}'")
    with open(samples_file) as f:
        samples_dict = json.load(f)
    return samples_dict


def get_calibration_sample(
    sample: str,
    magnet: str,
    particle: str,
    samples_file: str,
    max_files: int = None,
) -> Dict[str, Any]:
    """Return a list of calibration files.

    Args:
        sample: Data sample name (Turbo18, etc.)
        magnet: Magnet polarity (up, down)
        particle: Particle type (K, pi, etc.)
        samples_file: File in which to look up the file list.
        max_files: Optional. The maximum number of files to get. Defaults to
            None (= all files).
    """
    samples_dict = get_calibration_samples(samples_file)

    magnet = "Mag" + magnet.capitalize()
    sample_name = "-".join([sample, magnet, particle])

    try:
        sample_dict = samples_dict[sample_name]
    except KeyError:
        log.error(f"Sample '{sample_name}' not found in {samples_file}")
        raise

    calibration_sample = sample_dict.copy()
    if "link" in sample_dict:
        del calibration_sample["link"]
        link = sample_dict["link"]
        try:
            calibration_sample["files"] = samples_dict[link]["files"]
        except KeyError:
            log.error(f"Linked sample '{link}' not found in {samples_file}")
            raise

    if max_files:
        calibration_sample["files"] = calibration_sample["files"][:max_files]
    return calibration_sample


def root_to_dataframe(
    path: str, tree_names: List[str], branches: List[str]
) -> pd.DataFrame:
    """Return DataFrame with requested branches from tree in ROOT file.

    Args:
        path: Path to the ROOT file; either file system path or URL, e.g.
            root:///eos/lhcb/file.root.
        tree_names: Names of trees inside the ROOT file to read.
        branches: Branches to put in the DataFrame.
    """
    tree = None
    # EOS sometimes fails with a message saying the operation expired. It is
    # intermittent and hard to replicate. See this related issue:
    # https://github.com/scikit-hep/uproot4/issues/351. To avoid PIDCalib2
    # completely failing in these cases, we skip the file with a warning
    # message if this happens.
    try:
        root_file = uproot.open(path)
    except OSError as err:
        if "Operation expired" in err.args[0]:
            log.error(
                f"Failed to open '{path}' because an XRootD operation expired; skipping"
            )
            print(err)
            return None  # type: ignore
        else:
            raise

    dfs = []
    for tree_name in tree_names:
        try:
            tree = root_file[tree_name]
            dfs.append(tree.arrays(branches, library="pd"))  # type: ignore
        except uproot.exceptions.KeyInFileError as exc:  # type: ignore
            similar_keys = utils.find_similar_strings(exc.key, list(aliases), 0.80)
            similar_keys += utils.find_similar_strings(
                exc.key, tree.keys(), 0.80  # type: ignore
            )
            similar_keys += utils.find_similar_strings(
                "probe_" + exc.key, tree.keys(), 0.80  # type: ignore
            )
            # Remove duplicates while preserving ordering
            similar_keys = list(dict.fromkeys(similar_keys))
            log.error(
                (
                    f"Branch '{exc.key}' not found; similar aliases and/or branches "
                    f"that exist: {similar_keys}"
                )
            )
            raise
        except OSError as err:
            if "Operation expired" in err.args[0]:
                log.error(
                    (
                        f"Failed to open '{path}' because an XRootD operation "
                        "expired; skipping"
                    )
                )
                print(err)
                return None  # type: ignore
            else:
                raise

    return pd.concat(dfs, ignore_index=True)  # type: ignore


def get_tree_paths(particle: str, sample: str) -> List[str]:
    """Return a list of internal ROOT paths to relevant trees in the files

    Args:
        particle: Particle type (K, pi, etc.)
        sample: Data sample name (Turbo18, etc.)
    """
    tree_paths = []
    if is_simple(sample):
        # Run 1 (and Run 2 Electron) files have a simple structure with a single
        # tree
        tree_paths.append("DecayTree")
    else:
        # Run 2 ROOT file structure with multiple trees
        for tuple_name in tuple_names[particle]:
            tree_paths.append(f"{tuple_name}/DecayTree")

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
