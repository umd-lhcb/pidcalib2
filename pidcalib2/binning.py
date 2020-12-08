from typing import List

import numpy as np

valid_particles = ["pi", "K", "p", "mu"]


def P_binning(par: str, low: float = 3000, high: float = 100000) -> List[float]:
    assert par in valid_particles
    bins = []
    if par in ["pi", "K", "p"]:
        bins.append(low)
        bins.append(9300)  # R1 kaon threshold
        bins.append(15600)  # R2 kaon threshold
        # Uniform bin boundaries
        uniform_bins = np.linspace(19000, high, 16).tolist()
        for b in uniform_bins:
            bins.append(b)
    elif par == "Mu":
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


def ETA_binning(par, low=1.5, high=5.0):
    bins = np.linspace(low, high, 5).tolist()
    return bins


def nTracks_binning(par, low=0, high=500):
    bins = [low, 50, 200, 300, high]
    return bins


def TRCHI2_binning(par, low=0.0, high=3.0):
    bins = np.linspace(low, high, 4).tolist()
    return bins


# Dict of binnings for each track type and variable
binnings = {}

binnings["pi"] = {
    "P": P_binning("pi"),
    "ETA": ETA_binning("pi"),
    "nTracks": nTracks_binning("pi"),
    "TRCHI2NDOF": TRCHI2_binning("pi"),
}

binnings["K"] = {
    "P": P_binning("K"),
    "ETA": ETA_binning("K"),
    "nTracks": nTracks_binning("K"),
    "TRCHI2NDOF": TRCHI2_binning("K"),
}
