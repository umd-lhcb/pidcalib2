import pytest

from pidcalib2 import ref_calib


def test_ref_calib():
    config = {
        "year": 2018,
        "magnet": "up",
        "bin_vars": '{"P": "P", "ETA": "ETA", "nTracks": "nTracks"}',
        "output_dir": "tests/data",
        "dry_run": True,
        "ref_file": "tests/data/ref_test_data.root",
        "ref_tree": "DecayTree",
        "ref_pars": '{"Bach": ["K", "DLLK > 4"]}',
    }
    assert ref_calib.ref_calib(config) == pytest.approx(0.8951140908087826)
