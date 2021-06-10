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

"""Convert pickled PIDCalib2 histograms to TH*D & save them in a ROOT file.

Only 1D, 2D, and 3D histograms are supported by ROOT. Attempting to convert
higher-dimensional histograms will result in an exception.
"""

import itertools
import pathlib
import pickle
import sys

import boost_histogram as bh
import ROOT

from . import utils


def convert_to_root_histo(
    name: str, bh_histo: bh.Histogram, bh_error_histo: bh.Histogram = None
):
    """Convert boost_histogram histogram to a ROOT histogram.

    Only 1D, 2D, and 3D histograms are supported by ROOT. Attempting to convert
    higher-dimensional histograms will result in an exception.

    Args:
        name: Name of the new ROOT histogram.
        bh_histo: The histogram to convert.
        bh_error_histo: Optional. Histogram from which to bin read errors. If
            none is supplied, ROOT will calculate Poisson errors.

    Returns:
        The converted ROOT histogram. Type depends on dimensionality.
    """
    histo = None
    if len(bh_histo.axes) == 1:
        histo = ROOT.TH1D(name, name, 3, 0, 1)
        histo.SetBins(bh_histo.axes[0].size, bh_histo.axes[0].edges)
        histo.GetXaxis().SetTitle(bh_histo.axes[0].metadata["name"])
    elif len(bh_histo.axes) == 2:
        histo = ROOT.TH2D(name, name, 3, 0, 1, 3, 0, 1)
        histo.SetBins(
            bh_histo.axes[0].size,
            bh_histo.axes[0].edges,
            bh_histo.axes[1].size,
            bh_histo.axes[1].edges,
        )
        histo.GetXaxis().SetTitle(bh_histo.axes[0].metadata["name"])
        histo.GetYaxis().SetTitle(bh_histo.axes[1].metadata["name"])
    elif len(bh_histo.axes) == 3:
        histo = ROOT.TH3D(name, name, 3, 0, 1, 3, 0, 1, 3, 0, 1)
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
        root_indices = [index + 1 for index in indices_tuple]
        histo.SetBinContent(histo.GetBin(*root_indices), bh_histo[indices_tuple])
        if bh_error_histo is not None:
            histo.SetBinError(
                histo.GetBin(*root_indices), bh_error_histo[indices_tuple]
            )

    return histo


def main():
    pkl_path = pathlib.Path(sys.argv[1])
    eff_histos = {}
    with open(pkl_path, "rb") as f:
        eff_histos["eff"] = pickle.load(f)
        eff_histos["passing"] = pickle.load(f)
        eff_histos["total"] = pickle.load(f)
        eff_histos["passing_sumw2"] = pickle.load(f)
        eff_histos["total_sumw2"] = pickle.load(f)
        eff_histos["error"] = utils.create_error_histogram(eff_histos)

        for item in eff_histos.values():
            assert isinstance(item, bh.Histogram)

        root_path = pkl_path.with_suffix(".root")
        root_file = ROOT.TFile(str(root_path), "RECREATE")

        eff_histo = convert_to_root_histo(
            "eff_histo", eff_histos["eff"], eff_histos["error"]
        )
        eff_histo.Write()

        passing_histo = convert_to_root_histo("passing_histo", eff_histos["passing"])
        passing_histo.Write()

        total_histo = convert_to_root_histo("total_histo", eff_histos["total"])
        total_histo.Write()

        root_file.Close()


if __name__ == "__main__":
    main()
