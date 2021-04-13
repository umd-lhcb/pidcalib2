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

from typing import List

import numpy as np
from logzero import logger as log

valid_particles = ["Pi", "K", "P", "Mu"]


def p_binning(particle: str, low: float = 3000, high: float = 100000) -> List[float]:
    """Return a binning for the momentum.

    Args:
        particle: Particle type ["Pi", "K", ...]
        low: Optional. Lowest momentum [MeV]. Defaults to 3000.
        high: Optional. Highest momentum [MeV]. Defaults to 100000.
    """
    if particle not in valid_particles:
        log.error(f"'{particle}' is not a valid particle for P binning")
        raise KeyError

    bins = []
    if particle in ["Pi", "K", "P"]:
        bins.append(low)
        bins.append(9300)  # R1 kaon threshold
        bins.append(15600)  # R2 kaon threshold
        # Uniform bin boundaries
        uniform_bins = np.linspace(19000, high, 16).tolist()  # type:ignore
        bins.extend(uniform_bins)
    elif particle == "Mu":
        bins = [
            low,
            6000,
            8000,
            10000,
            12000,
            14500,
            17500,
            21500,
            27000,
            32000,
            40000,
            60000,
            70000,
            high,
        ]
    return bins


def eta_binning(particle, low: float = 1.5, high: float = 5.0) -> List[float]:
    bins = np.linspace(low, high, 5).tolist()  # type:ignore
    return bins


def ntracks_binning(particle, low: float = 0, high: float = 500) -> List[float]:
    bins = [low, 50, 200, 300, high]
    return bins


def nspdhits_binning(particle, low: float = 0, high: float = 1000) -> List[float]:
    bins = [low, 200, 400, 600, 800, high]
    return bins


def trchi2_binning(particle, low: float = 0.0, high: float = 3.0) -> List[float]:
    bins = np.linspace(low, high, 4).tolist()  # type:ignore
    return bins


# Dict of binnings for each track type and variable
binnings = {}

binnings["Pi"] = {
    "P": p_binning("Pi"),
    "Brunel_P": p_binning("Pi"),
    "ETA": eta_binning("Pi"),
    "Brunel_ETA": eta_binning("Pi"),
    "nTracks": ntracks_binning("Pi"),
    "nTracks_Brunel": ntracks_binning("Pi"),
    "nSPDhits": nspdhits_binning("Pi"),
    "nSPDhits_Brunel": nspdhits_binning("Pi"),
    "TRCHI2NDOF": trchi2_binning("Pi"),
}

binnings["K"] = {
    "P": p_binning("K"),
    "Brunel_P": p_binning("K"),
    "ETA": eta_binning("K"),
    "Brunel_ETA": eta_binning("K"),
    "nTracks": ntracks_binning("K"),
    "nTracks_Brunel": ntracks_binning("K"),
    "nSPDhits": nspdhits_binning("K"),
    "nSPDhits_Brunel": nspdhits_binning("K"),
    "TRCHI2NDOF": trchi2_binning("K"),
}

binnings["Mu"] = {
    "P": p_binning("Mu"),
    "ETA": eta_binning("Mu"),
    "nTracks": ntracks_binning("Mu"),
    "nTracks_Brunel": ntracks_binning("Mu"),
    "nSPDhits": nspdhits_binning("Mu"),
    "nSPDhits_Brunel": nspdhits_binning("Mu"),
    "TRCHI2NDOF": trchi2_binning("Mu"),
}
