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
def test_extract_branches_to_dataframe():
    assert (
        utils.extract_branches_to_dataframe(
            [
                "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision15/PIDCALIB.ROOT/00064787/0000/00064787_00000037_1.pidcalib.root"  # noqa: E501
            ],
            "pi",
            ["probe_PIDK"],
        ).shape[0]
        == 18400
    )


def test_translate_pid_cuts_to_branch_cuts():
    assert utils.translate_pid_cuts_to_branch_cuts(
        "probe", ["DLLK<4", "ProbNNpi>3"]
    ) == ["probe_PIDK<4", "probe_MC15TuneV1_ProbNNpi>3"]

    assert utils.translate_pid_cuts_to_branch_cuts(
        "probe", ["DLLK < 4", "ProbNNpi > 3"]
    ) == ["probe_PIDK<4", "probe_MC15TuneV1_ProbNNpi>3"]


def test_make_hist():
    df = pd.read_csv("tests/data/test_data.csv", index_col=0)
    hist = utils.make_hist(df, "pi", ["P"])
    assert hist.size == 20
    assert hist.sum() == pytest.approx(71.55106080517815)
    assert hist[3] == pytest.approx(13.581349537355582)


def test_get_relevant_branch_names():
    assert utils.get_relevant_branch_names("probe", ["DLLK < 4"], ["P"]) == [
        "probe_sWeight",
        "probe_P",
        "probe_PIDK",
    ]

    assert utils.get_relevant_branch_names("probe", ["DLLp > 4"], ["P", "ETA"]) == [
        "probe_sWeight",
        "probe_P",
        "probe_ETA",
        "probe_PIDp",
    ]


def test_pid_cut_to_branch_name_and_cut():
    assert utils.pid_cut_to_branch_name_and_cut("probe", "DLLK > 4") == (
        "probe_PIDK",
        ">4",
    )

    with pytest.raises(KeyError):
        utils.pid_cut_to_branch_name_and_cut("probe", "DLLX > 4")


def test_create_eff_histograms():
    df = pd.read_csv("tests/data/test_data.csv", index_col=0)
    hists = utils.create_eff_histograms(df, "pi", ["probe_PIDK>4"], ["P"])
    assert hists["eff_probe_PIDK>4"].sum() / hists[
        "eff_probe_PIDK>4"
    ].size == pytest.approx(0.18751578358705173 / 20)


def test_dataframe_from_local_file():
    df = utils.dataframe_from_local_file("tests/data/test_data.csv")
    assert df.shape[0] == 99
    assert df["probe_sWeight"][0] == pytest.approx(1.1081801082842266)
