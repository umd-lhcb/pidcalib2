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

import os
from pathlib import Path

import pytest

from pidcalib2 import ref_calib

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_ref_calib():
    config = {
        "sample": "Turbo18",
        "magnet": "up",
        "bin_vars": '{"P": "P", "ETA": "ETA", "nTracks": "nTracks"}',
        "output_dir": str(Path(THIS_DIR, "test_data")),
        "dry_run": True,
        "ref_file": str(Path(THIS_DIR, "test_data/ref_test_data.root")),
        "ref_tree": "DecayTree",
        "ref_pars": '{"Bach": ["K", "DLLK > 4"]}',
    }
    assert ref_calib.ref_calib(config) == pytest.approx(0.8951140908087826)
    assert Path(THIS_DIR, "test_data/ref_test_data_PID_eff.root").exists()
    os.remove(Path(THIS_DIR, "test_data/ref_test_data_PID_eff.root"))

    with pytest.raises(SyntaxError):
        config = {
            "sample": "Turbo18",
            "magnet": "up",
            "bin_vars": '{"P", "ETA", "nTracks"}',
            "output_dir": str(Path(THIS_DIR, "test_data")),
            "dry_run": True,
            "ref_file": str(Path(THIS_DIR, "test_data/ref_test_data.root")),
            "ref_tree": "DecayTree",
            "ref_pars": '{"Bach": ["K", "DLLK > 4"]}',
        }
        ref_calib.ref_calib(config)

    with pytest.raises(SyntaxError):
        config = {
            "sample": "Turbo18",
            "magnet": "up",
            "bin_vars": '{"P": "P", "ETA": "ETA", "nTracks": "nTracks"}',
            "output_dir": str(Path(THIS_DIR, "test_data")),
            "dry_run": True,
            "ref_file": str(Path(THIS_DIR, "test_data/ref_test_data.root")),
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