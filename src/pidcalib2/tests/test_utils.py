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
import pickle
from pathlib import Path

import pandas as pd
import pytest

from pidcalib2 import utils

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def test_make_hist():
    df = pd.read_csv(Path(THIS_DIR, "test_data/cal_test_data.csv"), index_col=0)
    hist = utils.make_hist(df, "Pi", ["P"])
    assert hist.size == 20
    assert hist.sum() == pytest.approx(71.55106080517815)
    assert hist[3] == pytest.approx(13.581349537355582)


def test_create_eff_histograms():
    df = pd.read_csv(Path(THIS_DIR, "test_data/cal_test_data.csv"), index_col=0)

    particle = "Pi"
    pid_cut = "DLLK>4"
    bin_vars = ["P"]

    hists = {}
    hists["total"] = utils.make_hist(df, particle, bin_vars)
    hists["total_sumw2"] = utils.make_hist(df, particle, bin_vars, True)

    df_passing = df.query(pid_cut)
    hists[f"passing_{pid_cut}"] = utils.make_hist(df_passing, particle, bin_vars)
    hists[f"passing_sumw2_{pid_cut}"] = utils.make_hist(
        df_passing, particle, bin_vars, True
    )

    eff_hists = utils.create_eff_histograms(hists)
    assert eff_hists["eff_DLLK>4"].sum() == pytest.approx(0.18751578358705173)
    assert eff_hists["eff_DLLK>4"].size == 20


def test_get_per_event_effs():
    df_ref = pd.read_csv(Path(THIS_DIR, "test_data/ref_test_data.csv"), index_col=0)
    prefixes = ["Bach"]
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    hists = {}
    hists["Bach"] = {}
    with open(
        Path(THIS_DIR, "test_data/effhists-Turbo18-up-K-DLLK>4-P.ETA.nTracks.pkl"), "rb"
    ) as f:
        hists["Bach"]["eff"] = pickle.load(f)
        hists["Bach"]["passing"] = pickle.load(f)
        hists["Bach"]["total"] = pickle.load(f)
        hists["Bach"]["passing_sumw2"] = pickle.load(f)
        hists["Bach"]["total_sumw2"] = pickle.load(f)
    df_ref = utils.add_bin_indices(df_ref, prefixes, bin_vars, hists)
    with pytest.warns(RuntimeWarning):
        df_ref = utils.add_efficiencies(df_ref, prefixes, hists)
    assert df_ref.PIDCalibEff.mean() == pytest.approx(0.8745793984247746)
    assert df_ref.PIDCalibErr.mean() == pytest.approx(0.010109849956086676)


def test_get_multiparticle_per_event_effs():
    df_ref = pd.read_csv(Path(THIS_DIR, "test_data/ref_test_data.csv"), index_col=0)
    prefixes = ["h1", "h2"]
    bin_vars = {"P": "P"}
    hists = {}
    hists["h1"] = {}
    with open(
        Path(THIS_DIR, "test_data/effhists-Turbo18-up-K-DLLK>0-P.pkl"), "rb"
    ) as f:
        hists["h1"]["eff"] = pickle.load(f)
        hists["h1"]["passing"] = pickle.load(f)
        hists["h1"]["total"] = pickle.load(f)
        hists["h1"]["passing_sumw2"] = pickle.load(f)
        hists["h1"]["total_sumw2"] = pickle.load(f)

    hists["h2"] = {}
    with open(
        Path(THIS_DIR, "test_data/effhists-Turbo18-up-Pi-DLLK<0-P.pkl"), "rb"
    ) as f:
        hists["h2"]["eff"] = pickle.load(f)
        hists["h2"]["passing"] = pickle.load(f)
        hists["h2"]["total"] = pickle.load(f)
        hists["h2"]["passing_sumw2"] = pickle.load(f)
        hists["h2"]["total_sumw2"] = pickle.load(f)

    df_ref = utils.add_bin_indices(df_ref, prefixes, bin_vars, hists)
    df_ref = utils.add_efficiencies(df_ref, prefixes, hists)
    assert df_ref.PIDCalibEff.mean() == pytest.approx(0.928205513381871)
    assert df_ref.PIDCalibErr.mean() == pytest.approx(0.004805856952295592)


def test_add_bin_indices():
    df_ref = pd.read_csv(Path(THIS_DIR, "test_data/ref_test_data.csv"), index_col=0)
    prefixes = ["Bach"]
    bin_vars = {"P": "P", "ETA": "ETA", "nTracks": "nTracks"}
    hists = {}
    hists["Bach"] = {}
    with open(
        Path(THIS_DIR, "test_data/effhists-Turbo18-up-K-DLLK>4-P.ETA.nTracks.pkl"), "rb"
    ) as f:
        hists["Bach"]["eff"] = pickle.load(f)
        hists["Bach"]["passing"] = pickle.load(f)
        hists["Bach"]["total"] = pickle.load(f)
        hists["Bach"]["passing_sumw2"] = pickle.load(f)
        hists["Bach"]["total_sumw2"] = pickle.load(f)
    df_ref = utils.add_bin_indices(df_ref, prefixes, bin_vars, hists)
    assert df_ref["Bach_P_PIDCalibBin"].sum() == 623
    assert df_ref["Bach_ETA_PIDCalibBin"].sum() == 120
    assert df_ref["Bach_PIDCalibBin"].sum() == 10438
