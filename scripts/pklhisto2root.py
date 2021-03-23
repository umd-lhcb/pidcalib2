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

import itertools
import pathlib
import pickle
import sys

import boost_histogram as bh
import ROOT


def convert_to_root_histo(bh_histo):
    histo = None
    if len(bh_histo.axes) == 1:
        histo = ROOT.TH1D("eff_histo", "Efficiency Histogram", 3, 0, 1)  # type: ignore
        histo.SetBins(bh_histo.axes[0].size, bh_histo.axes[0].edges)
        histo.GetXaxis().SetTitle(bh_histo.axes[0].metadata["name"])
    elif len(bh_histo.axes) == 2:
        histo = ROOT.TH2D(  # type: ignore
            "eff_histo", "Efficiency Histogram", 3, 0, 1, 3, 0, 1
        )
        histo.SetBins(
            bh_histo.axes[0].size,
            bh_histo.axes[0].edges,
            bh_histo.axes[1].size,
            bh_histo.axes[1].edges,
        )
        histo.GetXaxis().SetTitle(bh_histo.axes[0].metadata["name"])
        histo.GetYaxis().SetTitle(bh_histo.axes[1].metadata["name"])
    elif len(bh_histo.axes) == 3:
        histo = ROOT.TH3D(  # type: ignore
            "eff_histo", "Efficiency Histogram", 3, 0, 1, 3, 0, 1, 3, 0, 1
        )
        histo.SetBins(
            bh_histo.axes[0].size,
            bh_histo.axes[0].edges,
            bh_histo.axes[1].size,
            bh_histo.axes[1].edges,
            bh_histo.axes[2].size,
            bh_histo.axes[2].edges,
        )
        histo.GetXaxis().SetTitle(bh_histo.axes[0].metadata["name"])
        histo.GetYaxis().SetTitle(bh_histo.axes[1].metadata["name"])
        histo.GetZaxis().SetTitle(bh_histo.axes[2].metadata["name"])
    else:
        raise Exception(f"{len(bh_histo.axes)}D histograms not supported by ROOT")

    indices_ranges = [list(range(n)) for n in bh_histo.axes.size]
    for indices_tuple in itertools.product(*indices_ranges):
        root_indices = [index + 1 for index in indices_tuple]  # type: ignore
        histo.SetBinContent(histo.GetBin(*root_indices), bh_histo[indices_tuple])

    return histo


def main():
    pkl_path = pathlib.Path(sys.argv[1])
    with open(pkl_path, "rb") as f:
        bh_histo = pickle.load(f)
        if isinstance(bh_histo, bh.Histogram):
            histo = convert_to_root_histo(bh_histo)
            root_path = pkl_path.with_suffix(".root")
            root_file = ROOT.TFile(str(root_path), "RECREATE")  # type: ignore
            histo.Write()
            root_file.Close()


if __name__ == "__main__":
    main()
