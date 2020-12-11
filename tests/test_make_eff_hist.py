import shutil

import pandas as pd
import pytest
from pidcalib2 import make_eff_hists


def test_make_eff_hists():
    config = {
        "year": 2018,
        "magnet": "up",
        "particle": "pi",
        "pid_cuts": ["DLLK < 4", "DLLK<3"],
        "bin_vars": ["P"],
        "output_dir": "tests/test_output",
        "local_dataframe": "tests/data/test_data.csv",
    }
    make_eff_hists.make_eff_hists(config)
    eff_histo = pd.read_pickle("tests/test_output/effhist_2018_up_pi_DLLK<4_P.pkl")

    # These asserts might come in handy when some detail in the boost_histogram
    # implementation changes, thus failing the reference histogram comparison.
    assert eff_histo.sum(flow=False) == pytest.approx(17.812484216412948)
    assert eff_histo.sum(flow=True) == pytest.approx(23.57053953265155)
    assert eff_histo[1] == 0.8978157343833105

    eff_histo_reference = pd.read_pickle("tests/data/effhist_2018_up_pi_probe_PIDK<4.pkl")
    assert eff_histo == eff_histo_reference

    shutil.rmtree("tests/test_output")
