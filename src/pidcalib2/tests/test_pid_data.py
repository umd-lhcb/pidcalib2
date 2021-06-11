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

import math
import os
from pathlib import Path

import pytest

from pidcalib2 import pid_data

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


@pytest.mark.xrootd
@pytest.mark.slow
def test_root_to_dataframe():
    assert (
        pid_data.root_to_dataframe(
            "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision15/PIDCALIB.ROOT/00064787/0000/00064787_00000037_1.pidcalib.root",  # noqa: E501
            ["DSt_PiMTuple/DecayTree", "DSt_PiPTuple/DecayTree"],
            ["probe_PIDK"],
        ).shape[0]
        == 18400
    )


def test_get_tree_paths():
    assert pid_data.get_tree_paths("Pi", "26") == ["DecayTree"]
    assert pid_data.get_tree_paths("Pi", "Turbo15") == [
        "DSt_PiMTuple/DecayTree",
        "DSt_PiPTuple/DecayTree",
    ]
    assert pid_data.get_tree_paths("Mu", "Turbo15-MagUp-Mu") == [
        "Jpsi_MuPTuple/DecayTree",
        "Jpsi_MuMTuple/DecayTree",
    ]


def test_get_relevant_branch_names():
    assert pid_data.get_relevant_branch_names(["DLLK < 4"], ["P"]) == {
        "sWeight": "probe_sWeight",
        "P": "probe_P",
        "DLLK": "probe_PIDK",
    }

    assert pid_data.get_relevant_branch_names(["DLLp > 4"], ["P", "ETA"]) == {
        "sWeight": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert pid_data.get_relevant_branch_names(["DLLp == 4"], ["P", "ETA"]) == {
        "sWeight": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert pid_data.get_relevant_branch_names(["DLLp != 4"], ["P", "ETA"]) == {
        "sWeight": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert pid_data.get_relevant_branch_names(["test_DLL != 4"], ["special_var"]) == {
        "sWeight": "probe_sWeight",
        "test_DLL": "test_DLL",
        "special_var": "special_var",
    }


def test_dataframe_from_local_file():
    df = pid_data.dataframe_from_local_file(
        str(Path(THIS_DIR, "test_data/cal_test_data.csv")), ["sWeight"]
    )
    assert df.shape[0] == 99
    assert df["sWeight"][0] == pytest.approx(1.1081801082842266)

    df_pkl = pid_data.dataframe_from_local_file(
        str(Path(THIS_DIR, "test_data/cal_test_data.pkl")), ["sw"]
    )
    assert df_pkl.shape[0] == 99
    assert df_pkl["sw"][0] == pytest.approx(1.1081801082842266)

    with pytest.raises(KeyError):
        pid_data.dataframe_from_local_file(
            str(Path(THIS_DIR, "test_data/cal_test_data.csv")),
            ["this key doesn't exist"],
        )

    with pytest.raises(Exception):
        pid_data.dataframe_from_local_file(
            str(Path(THIS_DIR, "test_data/file.root")), ["sWeight"]
        )


def test_get_reference_branch_names():
    ref_pars = {"Bach": ["K", "DLLK > 4"]}
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    branch_names = pid_data.get_reference_branch_names(ref_pars, bin_vars)
    assert branch_names == ["Bach_P", "Bach_ETA", "nTracks"]


def test_get_reference_branch_name():
    assert pid_data.get_reference_branch_name("Bach", "foo", "bar") == "Bach_bar"
    assert pid_data.get_reference_branch_name("Bach", "nTracks", "bar") == "bar"
    assert pid_data.get_reference_branch_name("Bach", "nSPDhits", "bar") == "bar"
    assert (
        pid_data.get_reference_branch_name("Bach", "foo", "nTracks") == "Bach_nTracks"
    )


def test_get_calib_hists():
    ref_pars = {"Bach": ["K", "DLLK > 4"]}
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    eff_histos = pid_data.get_calib_hists(
        str(Path(THIS_DIR, "test_data")), "Turbo18", "up", ref_pars, bin_vars
    )
    assert math.isnan(eff_histos["Bach"]["eff"].sum())  # type: ignore
    assert eff_histos["Bach"]["eff"][4, 2, 1] == pytest.approx(0.9672206239134171)
    assert eff_histos["Bach"]["passing"].sum() == pytest.approx(103762743.76872134)
    assert eff_histos["Bach"]["total"].sum() == pytest.approx(114809531.81162772)
    assert eff_histos["Bach"]["passing_sumw2"].sum() == pytest.approx(
        130492602.73047148
    )
    assert eff_histos["Bach"]["total_sumw2"].sum() == pytest.approx(149177697.03949496)

    with pytest.raises(FileNotFoundError):
        pid_data.get_calib_hists(
            str(Path(THIS_DIR, "test_data")),
            "Turbo34",
            "up",
            ref_pars,
            bin_vars,
        )


def test_get_calibration_sample():
    assert (
        len(
            pid_data.get_calibration_sample(
                "Turbo18", "up", "Pi", str(Path(THIS_DIR, "../data/samples.json"))
            )["files"]
        )
        == 428
    )
    assert (
        len(
            pid_data.get_calibration_sample(
                "Turbo18",
                "up",
                "Pi",
                str(Path(THIS_DIR, "../data/samples.json")),
                3,
            )["files"]
        )
        == 3
    )
    assert (
        len(
            pid_data.get_calibration_sample(
                "Electron18",
                "down",
                "e_B_Jpsi",
                str(Path(THIS_DIR, "../data/samples.json")),
            )["files"]
        )
        == 1
    )
    assert (
        len(
            pid_data.get_calibration_sample(
                "20", "down", "K", str(Path(THIS_DIR, "../data/samples.json"))
            )["files"]
        )
        == 72
    )

    with pytest.raises(KeyError):
        pid_data.get_calibration_sample(
            "Turbo34", "up", "Pi", str(Path(THIS_DIR, "../data/samples.json"))
        )
