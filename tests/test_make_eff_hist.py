import math
import shutil

import pandas as pd
import pytest

from pidcalib2 import make_eff_hists


def test_make_eff_hists():
    config = {
        "sample": "Turbo18",
        "magnet": "up",
        "particle": "pi",
        "pid_cuts": ["DLLK < 4", "DLLK<3"],
        "cuts": None,
        "bin_vars": ["P"],
        "output_dir": "tests/test_output",
        "local_dataframe": "tests/data/cal_test_data.csv",
        "verbose": False,
    }
    make_eff_hists.make_eff_hists(config)
    eff_histo = pd.read_pickle("tests/test_output/effhists_Turbo18_up_pi_DLLK<4_P.pkl")

    # These asserts might come in handy when some detail in the boost_histogram
    # implementation changes, thus failing the reference histogram comparison.
    assert eff_histo.sum(flow=False) == pytest.approx(17.812484216412948)
    assert eff_histo.sum(flow=True) == pytest.approx(23.57053953265155)
    assert eff_histo[1] == 0.8978157343833105

    eff_histo_reference = pd.read_pickle(
        "tests/data/effhists_Turbo18_up_pi_DLLK<4_P.pkl"
    )
    assert eff_histo == eff_histo_reference

    shutil.rmtree("tests/test_output")


def test_make_eff_hists_with_cuts():
    config = {
        "sample": "Turbo18",
        "magnet": "up",
        "particle": "pi",
        "pid_cuts": ["DLLK < 4", "DLLK<3"],
        "cuts": ["Dst_IPCHI2 < 100", "probe_TRACK_GHOSTPROB < 0.05"],
        "bin_vars": ["P"],
        "output_dir": "tests/test_output",
        "local_dataframe": "tests/data/cal_test_data.csv",
        "verbose": True,
    }
    make_eff_hists.make_eff_hists(config)
    eff_histo = pd.read_pickle("tests/test_output/effhists_Turbo18_up_pi_DLLK<4_P.pkl")

    # These asserts might come in handy when some detail in the boost_histogram
    # implementation changes, thus failing the reference histogram comparison.
    assert eff_histo.sum(flow=False) == pytest.approx(17.812484216412948)
    assert eff_histo.sum(flow=True) == pytest.approx(23.57053953265155)
    assert eff_histo[1] == 0.8978157343833105

    eff_histo_reference = pd.read_pickle(
        "tests/data/effhists_Turbo18_up_pi_DLLK<4_P.pkl"
    )
    assert eff_histo == eff_histo_reference

    shutil.rmtree("tests/test_output")

    # Test that stricter cuts have an effect
    config["cuts"] = ["Dst_IPCHI2 < 9", "probe_TRACK_GHOSTPROB < 0.01"]
    with pytest.warns(RuntimeWarning):
        make_eff_hists.make_eff_hists(config)
    eff_histo = pd.read_pickle("tests/test_output/effhists_Turbo18_up_pi_DLLK<4_P.pkl")

    # With stricter cuts, some bins in the efficiency histogram become empty
    assert math.isnan(eff_histo.sum(flow=False))
    assert math.isnan(eff_histo.sum(flow=True))
    assert eff_histo[1] == 0.806456926487749

    assert eff_histo != eff_histo_reference

    shutil.rmtree("tests/test_output")


@pytest.mark.xrootd
@pytest.mark.slow
def test_make_eff_hists_local_file_list():
    config = {
        "sample": "Turbo18",
        "magnet": "up",
        "particle": "pi",
        "pid_cuts": ["DLLK < 4", "DLLK<3"],
        "cuts": None,
        "bin_vars": ["P"],
        "local_dataframe": None,
        "output_dir": "tests/test_output",
        "file_list": "tests/data/Turbo18-MagUp-Pi",
        "verbose": False,
    }
    make_eff_hists.make_eff_hists(config)
    eff_histo = pd.read_pickle("tests/test_output/effhists_Turbo18_up_pi_DLLK<4_P.pkl")

    # These asserts might come in handy when some detail in the boost_histogram
    # implementation changes, thus failing the reference histogram comparison.
    assert eff_histo.sum(flow=False) == pytest.approx(17.378390835392615)
    assert eff_histo.sum(flow=True) == pytest.approx(833.8480366879824)
    assert eff_histo[1] == 0.9797018332573063

    shutil.rmtree("tests/test_output")


@pytest.mark.xrootd
@pytest.mark.slow
def test_make_eff_hists_file_list():
    config = {
        "sample": "Turbo18",
        "magnet": "up",
        "particle": "pi",
        "pid_cuts": ["DLLK < 4", "DLLK<3"],
        "cuts": None,
        "bin_vars": ["P"],
        "local_dataframe": None,
        "output_dir": "tests/test_output",
        "file_list": None,
        "file_list_dir": "tests/data",
        "max_files": 1,
        "verbose": False,
    }
    make_eff_hists.make_eff_hists(config)
    eff_histo = pd.read_pickle("tests/test_output/effhists_Turbo18_up_pi_DLLK<4_P.pkl")

    # These asserts might come in handy when some detail in the boost_histogram
    # implementation changes, thus failing the reference histogram comparison.
    assert eff_histo.sum(flow=False) == pytest.approx(17.378390835392615)
    assert eff_histo.sum(flow=True) == pytest.approx(833.8480366879824)
    assert eff_histo[1] == 0.9797018332573063

    shutil.rmtree("tests/test_output")


def test_decode_arguments():
    assert (
        make_eff_hists.decode_arguments(
            [
                "--sample=Turbo18",
                "-m=up",
                "--particle=K",
                "--pid-cut='DLLK>5'",
                "--bin-var=P",
            ]
        ).verbose
        is False
    )

    with pytest.raises(SystemExit):
        make_eff_hists.decode_arguments(
            [
                "--sample=Turbo18",
                "-m=up",
                "--particle=nu",
                "--pid-cut='DLLK>5'",
                "--bin-var=P",
            ]
        )

    with pytest.raises(SystemExit):
        make_eff_hists.decode_arguments(
            ["-m=up", "--particle=K", "--pid-cut='DLLK>5'", "--bin-var=P"]
        )

    with pytest.raises(SystemExit):
        make_eff_hists.decode_arguments(
            ["--sample=Turbo18", "-m=up", "--particle=K", "--bin-var=P"]
        )
