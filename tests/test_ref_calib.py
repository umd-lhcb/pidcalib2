import os
import pathlib

import pytest

from pidcalib2 import ref_calib


def test_ref_calib():
    config = {
        "sample": "Turbo18",
        "magnet": "up",
        "bin_vars": '{"P": "P", "ETA": "ETA", "nTracks": "nTracks"}',
        "output_dir": "tests/data",
        "dry_run": True,
        "ref_file": "tests/data/ref_test_data.root",
        "ref_tree": "DecayTree",
        "ref_pars": '{"Bach": ["K", "DLLK > 4"]}',
    }
    assert ref_calib.ref_calib(config) == pytest.approx(0.8951140908087826)
    assert pathlib.Path("tests/data/ref_test_data_PID_eff.root").exists()
    os.remove("tests/data/ref_test_data_PID_eff.root")

    with pytest.raises(SyntaxError):
        config = {
            "sample": "Turbo18",
            "magnet": "up",
            "bin_vars": '{"P", "ETA", "nTracks"}',
            "output_dir": "tests/data",
            "dry_run": True,
            "ref_file": "tests/data/ref_test_data.root",
            "ref_tree": "DecayTree",
            "ref_pars": '{"Bach": ["K", "DLLK > 4"]}',
        }
        ref_calib.ref_calib(config)

    with pytest.raises(SyntaxError):
        config = {
            "sample": "Turbo18",
            "magnet": "up",
            "bin_vars": '{"P": "P", "ETA": "ETA", "nTracks": "nTracks"}',
            "output_dir": "tests/data",
            "dry_run": True,
            "ref_file": "tests/data/ref_test_data.root",
            "ref_tree": "DecayTree",
            "ref_pars": '{"Bach"}',
        }
        ref_calib.ref_calib(config)


def test_decode_arguments():
    assert (
        ref_calib.decode_arguments(
            [
                "--sample=Turbo18",
                "-m=up",
                "--ref-file=ref_file.root",
                "--ref-pars={'Ref': ['pi', 'DLLK<5']}",
            ]
        ).magnet
        == "up"
    )

    with pytest.raises(SystemExit):
        ref_calib.decode_arguments(
            [
                "--sample=Turbo18",
                "-m=up",
                "--ref-file=ref_file.root",
                "--ref-pars={'Ref': ['pi', 'DLLK<5']}",
                "--particle=K",
            ]
        )

    with pytest.raises(SystemExit):
        ref_calib.decode_arguments(
            [
                "--sample=Turbo18",
                "-m=MagUp",
                "--ref-file=ref_file.root",
                "--ref-pars={'Ref': ['pi', 'DLLK<5']}",
            ]
        )

    with pytest.raises(SystemExit):
        ref_calib.decode_arguments(
            [
                "--sample=Turbo18",
                "-m=up",
                "--ref-file=ref_file.root",
            ]
        )
