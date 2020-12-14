from typing import List

import numpy as np
from logzero import logger as log

valid_particles = ["pi", "K", "p", "mu"]


def p_binning(particle: str, low: float = 3000, high: float = 100000) -> List[float]:
    """Return a binning for the momentum.

    Args:
        particle (str): Particle type ["pi", "K", "p"]
        low (float, optional): [description]. Defaults to 3000.
        high (float, optional): [description]. Defaults to 100000.

    Returns:
        List[float]: [description]
    """
    if particle not in valid_particles:
        log.error(f"'{particle}' is not a valid particle for P binning")
        raise KeyError

    bins = []
    if particle in ["pi", "K", "p"]:
        bins.append(low)
        bins.append(9300)  # R1 kaon threshold
        bins.append(15600)  # R2 kaon threshold
        # Uniform bin boundaries
        uniform_bins = np.linspace(19000, high, 16).tolist()
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
    bins = np.linspace(low, high, 5).tolist()
    return bins


def ntracks_binning(particle, low: float = 0, high: float = 500) -> List[float]:
    bins = [low, 50, 200, 300, high]
    return bins


def trchi2_binning(particle, low: float = 0.0, high: float = 3.0) -> List[float]:
    bins = np.linspace(low, high, 4).tolist()
    return bins


# Dict of binnings for each track type and variable
binnings = {}

binnings["pi"] = {
    "P": p_binning("pi"),
    "ETA": eta_binning("pi"),
    "nTracks": ntracks_binning("pi"),
    "TRCHI2NDOF": trchi2_binning("pi"),
}

binnings["K"] = {
    "P": p_binning("K"),
    "ETA": eta_binning("K"),
    "nTracks": ntracks_binning("K"),
    "TRCHI2NDOF": trchi2_binning("K"),
}
