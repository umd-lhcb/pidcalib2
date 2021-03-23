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

import pickle

import pandas as pd
import pytest

from pidcalib2 import utils


def test_make_hist():
    df = pd.read_csv("tests/data/cal_test_data.csv", index_col=0)
    hist = utils.make_hist(df, "Pi", ["P"])
    assert hist.size == 20
    assert hist.sum() == pytest.approx(71.55106080517815)
    assert hist[3] == pytest.approx(13.581349537355582)


def test_create_eff_histograms():
    df = pd.read_csv("tests/data/cal_test_data.csv", index_col=0)
    hists = utils.create_eff_histograms(df, "Pi", ["DLLK>4"], ["P"])
    assert hists["eff_DLLK>4"].sum() / hists["eff_DLLK>4"].size == pytest.approx(
        0.18751578358705173 / 20
    )


def test_get_per_event_effs():
    df_ref = pd.read_csv("tests/data/ref_test_data.csv", index_col=0)
    prefixes = ["Bach"]
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    hists = {}
    with open("tests/data/effhists_Turbo18_up_K_DLLK>4_P-ETA-nTracks.pkl", "rb") as f:
        hists["Bach"] = pickle.load(f)
    df_ref = utils.add_bin_indices(df_ref, prefixes, bin_vars, hists)
    df_ref = utils.add_efficiencies(df_ref, prefixes, hists)
    assert df_ref.PID_eff.mean() == pytest.approx(0.8951140908087826)


def test_add_bin_indices():
    df_ref = pd.read_csv("tests/data/ref_test_data.csv", index_col=0)
    prefixes = ["Bach"]
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    hists = {}
    with open("tests/data/effhists_Turbo18_up_K_DLLK>4_P-ETA-nTracks.pkl", "rb") as f:
        hists["Bach"] = pickle.load(f)
    df_ref = utils.add_bin_indices(df_ref, prefixes, bin_vars, hists)
    assert df_ref["Bach_P_PID_index"].sum() == 623
    assert df_ref["Bach_ETA_PID_index"].sum() == 120
    assert df_ref["Bach_PID_index"].sum() == 10438
