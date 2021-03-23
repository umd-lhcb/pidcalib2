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
import shutil

import pytest
import ROOT

from pidcalib2 import merge_trees


@pytest.mark.pyroot
def test_copy_tree():
    shutil.copy("tests/data/ref_test_data.root", "tests/data/ref_test_data_copy.root")

    merge_trees.copy_tree_and_set_as_friend(
        "tests/data/ref_PID_eff.root",
        "PID_eff_tree",
        "tests/data/ref_test_data_copy.root",
        "DecayTree",
    )

    new_file = ROOT.TFile("tests/data/ref_test_data_copy.root")
    new_tree = new_file.Get("DecayTree")

    assert new_tree.GetEntry(0) == 68
    assert new_tree.PID_eff == pytest.approx(0.9933816230604905)

    os.remove("tests/data/ref_test_data_copy.root")


@pytest.mark.pyroot
def test_copy_tree_bad_filenames():
    with pytest.raises(SystemExit):
        merge_trees.copy_tree_and_set_as_friend(
            "tests/data/x.root",
            "PID_eff_tree",
            "tests/data/ref_test_data_copy.root",
            "DecayTree",
        )

    with pytest.raises(SystemExit):
        merge_trees.copy_tree_and_set_as_friend(
            "tests/data/ref_PID_eff.root",
            "PID_eff_tree",
            "tests/data/x.root",
            "DecayTree",
        )


@pytest.mark.pyroot
def test_copy_tree_bad_treenames():
    shutil.copy("tests/data/ref_test_data.root", "tests/data/ref_test_data_copy.root")

    with pytest.raises(SystemExit):
        merge_trees.copy_tree_and_set_as_friend(
            "tests/data/ref_PID_eff.root",
            "x",
            "tests/data/ref_test_data_copy.root",
            "DecayTree",
        )

    with pytest.raises(SystemExit):
        merge_trees.copy_tree_and_set_as_friend(
            "tests/data/ref_PID_eff.root",
            "PID_eff_tree",
            "tests/data/ref_test_data_copy.root",
            "x",
        )

    os.remove("tests/data/ref_test_data_copy.root")
