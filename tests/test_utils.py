import pickle

import pandas as pd
import pytest
from pidcalib2 import utils


def test_pid_calib_sample_dir():
    assert utils.pidcalib_sample_dir(2018, "up") == "Collision18/PIDCALIB.ROOT/00082947"

    with pytest.raises(AssertionError):
        utils.pidcalib_sample_dir(2018, "middle")

    with pytest.raises(AssertionError):
        utils.pidcalib_sample_dir(3018, "up")


@pytest.mark.xrootd
def test_get_eos_paths():
    assert len(utils.get_eos_paths(2018, "up")) == 419
    assert len(utils.get_eos_paths(2018, "down")) == 423
    assert len(utils.get_eos_paths(2017, "up")) == 371
    assert len(utils.get_eos_paths(2017, "down")) == 477
    assert len(utils.get_eos_paths(2016, "up")) == 146
    assert len(utils.get_eos_paths(2016, "down")) == 154
    assert len(utils.get_eos_paths(2015, "up")) == 43
    assert len(utils.get_eos_paths(2015, "down")) == 78
    assert (
        utils.get_eos_paths(2018, "up")[0]
        == "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision18/PIDCALIB.ROOT/00082947/0000/00082947_00000001_1.pidcalib.root"  # noqa: E501
    )


@pytest.mark.xrootd
@pytest.mark.slow
def test_calib_root_to_dataframe():
    assert (
        utils.calib_root_to_dataframe(
            [
                "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision15/PIDCALIB.ROOT/00064787/0000/00064787_00000037_1.pidcalib.root"  # noqa: E501
            ],
            "pi",
            {"DLLK": "probe_PIDK"},
        ).shape[0]
        == 18400
    )


def test_make_hist():
    df = pd.read_csv("tests/data/cal_test_data.csv", index_col=0)
    hist = utils.make_hist(df, "pi", ["P"])
    assert hist.size == 20
    assert hist.sum() == pytest.approx(71.55106080517815)
    assert hist[3] == pytest.approx(13.581349537355582)


def test_get_relevant_branch_names():
    assert utils.get_relevant_branch_names("probe", ["DLLK < 4"], ["P"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "DLLK": "probe_PIDK",
    }

    assert utils.get_relevant_branch_names("probe", ["DLLp > 4"], ["P", "ETA"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert utils.get_relevant_branch_names("probe", ["DLLp == 4"], ["P", "ETA"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert utils.get_relevant_branch_names("probe", ["DLLp != 4"], ["P", "ETA"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    with pytest.raises(ValueError):
        utils.get_relevant_branch_names("probe", ["DLLp = 4"], ["P", "ETA"])


def test_create_eff_histograms():
    df = pd.read_csv("tests/data/cal_test_data.csv", index_col=0)
    hists = utils.create_eff_histograms(df, "pi", ["DLLK>4"], ["P"])
    assert hists["eff_DLLK>4"].sum() / hists["eff_DLLK>4"].size == pytest.approx(
        0.18751578358705173 / 20
    )


def test_dataframe_from_local_file():
    df = utils.dataframe_from_local_file("tests/data/cal_test_data.csv", ["sw"])
    assert df.shape[0] == 99
    assert df["sw"][0] == pytest.approx(1.1081801082842266)

    with pytest.raises(KeyError):
        utils.dataframe_from_local_file(
            "tests/data/cal_test_data.csv", ["this key doesn't exist"]
        )


def test_get_per_event_effs():
    df_ref = pd.read_csv("tests/data/ref_test_data.csv", index_col=0)
    ref_pars = {"Bach": ["K", "DLLK > 4"]}
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    hists = {}
    with open("tests/data/effhist_2018_up_K_DLLK>4_P_ETA_nTracks.pkl", "rb") as f:
        hists["Bach"] = pickle.load(f)
    df_ref = utils.add_bin_indexes(df_ref, ref_pars, bin_vars, hists)
    df_ref = utils.add_efficiencies(df_ref, ref_pars, hists)
    assert df_ref.eff.mean() == pytest.approx(0.8951140908087826)
