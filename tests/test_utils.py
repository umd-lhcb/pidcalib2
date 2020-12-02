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
        == "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision18/PIDCALIB.ROOT/00082947/0000/00082947_00000001_1.pidcalib.root"
    )


@pytest.mark.xrootd
@pytest.mark.slow
def test_extract_branches_to_dataframe():
    assert (
        utils.extract_branches_to_dataframe(
            [
                "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision15/PIDCALIB.ROOT/00064787/0000/00064787_00000037_1.pidcalib.root"
            ],
            "pi",
            ["probe_PIDK"],
        ).shape[0]
        == 18400
    )


def test_translate_pid_cuts_to_branch_cuts():
    assert utils.translate_pid_cuts_to_branch_cuts(
        "probe", ["DLLK<4", "ProbNNPi>3"], ["P", "ETA"]
    ) == ["probe_PIDK<4", "ProbNNPi>3"]

    assert utils.translate_pid_cuts_to_branch_cuts(
        "probe", ["DLLK<4", "ProbNNPi>3"], ["ETA"]
    ) == ["probe_PIDK<4", "ProbNNPi>3"]

    assert utils.translate_pid_cuts_to_branch_cuts(
        "probe", ["DLLK < 4", "ProbNNPi > 3"], ["ETA"]
    ) == ["probe_PIDK<4", "ProbNNPi>3"]
