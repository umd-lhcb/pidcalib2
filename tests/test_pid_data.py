import pytest
from pidcalib2 import pid_data


def test_pid_calib_sample_dir():
    assert (
        pid_data.pidcalib_sample_dir(2018, "up") == "Collision18/PIDCALIB.ROOT/00109276"
    )

    with pytest.raises(AssertionError):
        pid_data.pidcalib_sample_dir(2018, "middle")

    with pytest.raises(AssertionError):
        pid_data.pidcalib_sample_dir(3018, "up")


@pytest.mark.xrootd
def test_get_eos_paths():
    assert len(pid_data.get_eos_paths(2018, "up")) == 428
    assert len(pid_data.get_eos_paths(2018, "down")) == 401
    assert len(pid_data.get_eos_paths(2017, "up")) == 312
    assert len(pid_data.get_eos_paths(2017, "down")) == 371
    assert len(pid_data.get_eos_paths(2016, "up")) == 238
    assert len(pid_data.get_eos_paths(2016, "down")) == 258
    assert len(pid_data.get_eos_paths(2015, "up")) == 43
    assert len(pid_data.get_eos_paths(2015, "down")) == 78
    assert len(pid_data.get_eos_paths(2015, "down", 10)) == 10
    assert (
        pid_data.get_eos_paths(2018, "up")[0]
        == "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision18/PIDCALIB.ROOT/00109276/0000/00109276_00000001_1.pidcalib.root"  # noqa: E501
    )


@pytest.mark.xrootd
@pytest.mark.slow
def test_calib_root_to_dataframe():
    assert (
        pid_data.calib_root_to_dataframe(
            [
                "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision15/PIDCALIB.ROOT/00064787/0000/00064787_00000037_1.pidcalib.root"  # noqa: E501
            ],
            ["DSt_PiMTuple/DecayTree", "DSt_PiPTuple/DecayTree"],
            {"DLLK": "probe_PIDK"},
        ).shape[0]
        == 18400
    )


def test_get_tree_paths():
    assert pid_data.get_tree_paths("pi", 2014) == ["DecayTree"]
    assert pid_data.get_tree_paths("pi", 2015) == [
        "DSt_PiPTuple/DecayTree",
        "DSt_PiMTuple/DecayTree",
    ]
    assert pid_data.get_tree_paths("mu", 2015) == [
        "Jpsi_MuPTuple/DecayTree",
        "Jpsi_MuMTuple/DecayTree",
    ]


def test_get_relevant_branch_names():
    assert pid_data.get_relevant_branch_names("probe", ["DLLK < 4"], ["P"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "DLLK": "probe_PIDK",
    }

    assert pid_data.get_relevant_branch_names("probe", ["DLLp > 4"], ["P", "ETA"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert pid_data.get_relevant_branch_names("probe", ["DLLp == 4"], ["P", "ETA"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert pid_data.get_relevant_branch_names("probe", ["DLLp != 4"], ["P", "ETA"]) == {
        "sw": "probe_sWeight",
        "P": "probe_P",
        "ETA": "probe_ETA",
        "DLLp": "probe_PIDp",
    }

    assert pid_data.get_relevant_branch_names(
        "notprobe", ["DLLp != 4"], ["P", "ETA"]
    ) == {"P": "notprobe_P", "ETA": "notprobe_ETA", "DLLp": "notprobe_PIDp"}

    with pytest.raises(ValueError):
        pid_data.get_relevant_branch_names("probe", ["DLLp = 4"], ["P", "ETA"])


def test_dataframe_from_local_file():
    df = pid_data.dataframe_from_local_file("tests/data/cal_test_data.csv", ["sw"])
    assert df.shape[0] == 99
    assert df["sw"][0] == pytest.approx(1.1081801082842266)

    with pytest.raises(KeyError):
        pid_data.dataframe_from_local_file(
            "tests/data/cal_test_data.csv", ["this key doesn't exist"]
        )

    with pytest.raises(Exception):
        pid_data.dataframe_from_local_file("tests/data/file.root", ["sw"])


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
    eff_histos = pid_data.get_calib_hists("tests/data", 2018, "up", ref_pars, bin_vars)
    assert eff_histos["Bach"].sum() == pytest.approx(199.01361598888047)

    with pytest.raises(FileNotFoundError):
        pid_data.get_calib_hists("tests/data", 18, "up", ref_pars, bin_vars)
