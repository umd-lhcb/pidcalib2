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

import json
from typing import List

import numpy as np
from logzero import logger as log

valid_particles = ["Pi", "K", "P", "Mu", "e"]


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
    if particle in ["Pi", "K", "P", "e"]:
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
    "P": {"bin_edges": p_binning("Pi")},
    "Brunel_P": {"bin_edges": p_binning("Pi")},
    "ETA": {"bin_edges": eta_binning("Pi")},
    "Brunel_ETA": {"bin_edges": eta_binning("Pi")},
    "nTracks": {"bin_edges": ntracks_binning("Pi")},
    "nTracks_Brunel": {"bin_edges": ntracks_binning("Pi")},
    "nSPDhits": {"bin_edges": nspdhits_binning("Pi")},
    "nSPDhits_Brunel": {"bin_edges": nspdhits_binning("Pi")},
    "TRCHI2NDOF": {"bin_edges": trchi2_binning("Pi")},
}

binnings["K"] = {
    "P": {"bin_edges": p_binning("K")},
    "Brunel_P": {"bin_edges": p_binning("K")},
    "ETA": {"bin_edges": eta_binning("K")},
    "Brunel_ETA": {"bin_edges": eta_binning("K")},
    "nTracks": {"bin_edges": ntracks_binning("K")},
    "nTracks_Brunel": {"bin_edges": ntracks_binning("K")},
    "nSPDhits": {"bin_edges": nspdhits_binning("K")},
    "nSPDhits_Brunel": {"bin_edges": nspdhits_binning("K")},
    "TRCHI2NDOF": {"bin_edges": trchi2_binning("K")},
}

binnings["Mu"] = {
    "P": {"bin_edges": p_binning("Mu")},
    "Brunel_P": {"bin_edges": p_binning("Mu")},
    "ETA": {"bin_edges": eta_binning("Mu")},
    "Brunel_ETA": {"bin_edges": eta_binning("Mu")},
    "nTracks": {"bin_edges": ntracks_binning("Mu")},
    "nTracks_Brunel": {"bin_edges": ntracks_binning("Mu")},
    "nSPDhits": {"bin_edges": nspdhits_binning("Mu")},
    "nSPDhits_Brunel": {"bin_edges": nspdhits_binning("Mu")},
    "TRCHI2NDOF": {"bin_edges": trchi2_binning("Mu")},
}

binnings["P"] = {
    "P": {"bin_edges": p_binning("P")},
    "Brunel_P": {"bin_edges": p_binning("P")},
    "ETA": {"bin_edges": eta_binning("P")},
    "Brunel_ETA": {"bin_edges": eta_binning("P")},
    "nTracks": {"bin_edges": ntracks_binning("P")},
    "nTracks_Brunel": {"bin_edges": ntracks_binning("P")},
    "nSPDhits": {"bin_edges": nspdhits_binning("P")},
    "nSPDhits_Brunel": {"bin_edges": nspdhits_binning("P")},
    "TRCHI2NDOF": {"bin_edges": trchi2_binning("P")},
}

binnings["e"] = {
    "P": {"bin_edges": p_binning("e")},
    "Brunel_P": {"bin_edges": p_binning("e")},
    "ETA": {"bin_edges": eta_binning("e")},
    "Brunel_ETA": {"bin_edges": eta_binning("e")},
    "nTracks": {"bin_edges": ntracks_binning("e")},
    "nTracks_Brunel": {"bin_edges": ntracks_binning("e")},
    "nSPDhits": {"bin_edges": nspdhits_binning("e")},
    "nSPDhits_Brunel": {"bin_edges": nspdhits_binning("e")},
    "TRCHI2NDOF": {"bin_edges": trchi2_binning("e")},
}


def set_binning(particle: str, variable: str, bin_edges: List[float]) -> None:
    """Set a new binning for a variable of a particle.

    Either a binning for a new particle/variable is added or the existing
    binning is rewritten.

    Args:
        particle: Particle name.
        variable: Variable name, e.g., "P" or "Brunel_ETA"
        bin_edges: A list of all bin edges.
    """
    if not isinstance(bin_edges, list):
        log.error("bin_edges parameter is not a list.")
        raise TypeError

    if particle not in binnings:
        binnings[particle] = {}

    binnings[particle][variable] = {"bin_edges": bin_edges}


def get_binning(particle: str, variable: str, verbose: bool = False) -> List[float]:
    """Return a suitable binning for a particle and variable.

    Args:
        particle: Particle name.
        variable: Variable name, e.g., "P" or "Brunel_ETA"
        verbose: Optional. Print message when alternative binning is used.
            Defaults to False.
    """
    if particle not in binnings or variable not in binnings[particle]:
        # Remove particle suffix, e.g., 'DsPhi' in 'K_DsPhi'
        pure_particle = particle.split("_", 1)[0]
        if pure_particle not in binnings or variable not in binnings[pure_particle]:
            log.error(f"No '{variable}' binning defined for particle {particle}")
            raise KeyError
        else:
            if verbose:
                log.info(
                    (
                        f"No '{variable}' binning defined for particle '{particle}'. "
                        f"Falling back to particle '{pure_particle}' binning."
                    )
                )
            return binnings[pure_particle][variable]["bin_edges"]

    else:
        return binnings[particle][variable]["bin_edges"]


def load_binnings(path: str) -> None:
    """Load binnings from a JSON file."""
    new_binnings = {}
    log.info(f"Loading binnings from {path}")
    with open(path) as f:
        new_binnings = json.load(f)
    for particle, variables in new_binnings.items():
        for variable, binning in variables.items():
            set_binning(particle, variable, binning)
