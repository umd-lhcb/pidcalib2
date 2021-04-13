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

import re
from typing import Dict, List, Union

import boost_histogram as bh
import numpy as np
import pandas as pd
from logzero import logger as log

from . import binning, pid_data


def make_hist(
    df: pd.DataFrame, particle: str, bin_vars: List[str], square_weights: bool = False
) -> bh.Histogram:
    """Create a histogram of sWeighted events with appropriate binning

    Args:
        df: DataFrame from which to histogram events.
        particle: Particle type (K, Pi, etc.).
        bin_vars: Binning variables in the user-convention, e.g., ["P", "ETA"].
        square_weights: Use square of sWeights instead of sWeights as the weights.
    """
    axis_list = []
    vals_list = []

    # Loop over bin dimensions and define the axes
    for bin_var in bin_vars:
        if (
            particle not in binning.binnings
            or bin_var not in binning.binnings[particle]
        ):
            log.error(f"No binning defined for particle {particle} variable {bin_var}")

        axis_list.append(
            bh.axis.Variable(
                binning.binnings[particle][bin_var], metadata={"name": bin_var}
            )
        )
        vals = df[bin_var].values
        vals_list.append(vals)

    # Create boost-histogram with the desired axes, and fill with sWeight applied
    hist = bh.Histogram(*axis_list)
    if square_weights:
        hist.fill(*vals_list, weight=np.square(df["sWeight"]))  # type: ignore
    else:
        hist.fill(*vals_list, weight=df["sWeight"])

    return hist


def create_eff_histograms(
    df_total: pd.DataFrame, particle: str, pid_cuts: List[str], bin_vars: List[str]
) -> Dict[str, bh.Histogram]:
    """Create efficiency histograms for all supplied PID cuts.

    Args:
        df_total: DataFrame with all events.
        particle: Particle type (K, Pi, etc.).
        pid_cuts: Simple user-level cut list, e.g., ["DLLK<4"].
        bin_vars: Variables used for binning.

    Returns:
        A dictionary with all the efficiency histograms, with the PID cuts as
        keys.
    """

    hists = {}
    hists["total"] = make_hist(df_total, particle, bin_vars)
    hists["total_sumw2"] = make_hist(df_total, particle, bin_vars, True)

    zero_bins = np.count_nonzero(hists["total"].view(flow=False) == 0)
    if zero_bins:
        log.warning(f"There are {zero_bins} empty bins in the total histogram!")
        log.warning(hists["total"].view(flow=False))
        hist_total_nan = hists["total"].view()
        hist_total_nan[hist_total_nan == 0] = np.nan
        hists["total"][...] = hist_total_nan

    n_total = len(df_total.index)
    for i, pid_cut in enumerate(pid_cuts):
        log.info(f"Processing '{pid_cuts[i]}' cut")
        df_passing = df_total.query(pid_cut)
        n_passing = len(df_passing.index)
        percent_passing = n_passing / n_total * 100
        log.info(
            f"{n_passing}/{n_total} ({percent_passing:.2f}%) events passed the cut"
        )
        hists[f"passing_{pid_cut}"] = make_hist(df_passing, particle, bin_vars)
        hists[f"passing_sumw2_{pid_cut}"] = make_hist(
            df_passing, particle, bin_vars, True
        )
        log.debug("Created 'passing' histogram")

        hists[f"eff_{pid_cut}"] = hists[f"passing_{pid_cut}"].copy()
        hists[f"eff_{pid_cut}"][...] = hists[f"passing_{pid_cut}"].view(
            flow=False
        ) / hists["total"].view(flow=False)
        log.debug(f"Created 'eff_{pid_cut}' histogram")

    return hists


def log_config(config: dict) -> None:
    """Pretty-print a config/dict."""
    longest_key = len(max(config, key=len))
    log.info("=" * longest_key)
    for entry in config:
        if config[entry] is not None:
            log.info(f"{entry:{longest_key}}: {config[entry]}")
    log.info("=" * longest_key)


def add_bin_indices(
    df: pd.DataFrame,
    prefixes: List[str],
    bin_vars: Dict[str, str],
    eff_hists: Dict[str, Dict[str, bh.Histogram]],
) -> pd.DataFrame:
    """Return a DataFrame with added indices of bins for each event.

    The binnings of binning variables are taken from efficiency histograms.
    Each event falls into a certain bin in each binning variable. This bin's
    index is added to the DataFrame. The same procedure is repeated for each
    variable. Finally, a global index of the N-dimensional bin of the
    efficiency histogram where the event belongs is added. Multiple
    efficiency histograms can be specified since the total PID efficiency for
    the event can depend on multiple particles.

    Args:
        df: Input dataframe.
        prefixes: Branch prefixes of the particles in the reference sample.
        bin_vars: Variables used for binning.
        eff_hists: Efficiency histograms for each prefix/particle.
    """
    df_new = df.copy()
    for prefix in prefixes:
        eff_histo = eff_hists[prefix]["eff"]
        axes = [
            pid_data.get_reference_branch_name(
                prefix, axis.metadata["name"], bin_vars[axis.metadata["name"]]
            )
            for axis in eff_histo.axes
        ]
        for bin_var, branch_name in bin_vars.items():
            ref_branch_name = pid_data.get_reference_branch_name(
                prefix, bin_var, branch_name
            )
            bins = []
            for axis in eff_histo.axes:
                if axis.metadata["name"] == bin_var:
                    bins = axis.edges
            df_new[f"{ref_branch_name}_PIDCalibBin"] = pd.cut(
                df_new[ref_branch_name],
                bins,
                labels=False,
                include_lowest=True,
                right=False,
                precision=0,
            )

        df_nan = df_new[df_new.isna().any(axis=1)]  # type: ignore
        df_new.dropna(inplace=True)
        index_names = [f"{axis}_PIDCalibBin" for axis in axes]
        indices = np.ravel_multi_index(
            df_new[index_names].transpose().to_numpy().astype(int),  # type: ignore
            eff_histo.axes.size,
        )
        df_new[f"{prefix}_PIDCalibBin"] = indices
        df_new = pd.concat([df_new, df_nan]).sort_index()  # type: ignore
    log.debug("Bin indices assigned")
    return df_new  # type: ignore


def add_efficiencies(
    df: pd.DataFrame, prefixes: List[str], eff_hists: Dict[str, Dict[str, bh.Histogram]]
) -> pd.DataFrame:
    """Return a DataFrame with added efficiencies for each event.

    Each particle correspondig to a prefix is assigned an efficiency. The
    total event efficiency is the product of individual efficiencies.

    Args:
        df: Input dataframe.
        prefixes: Branch prefixes of the particles in the reference sample.
        eff_hists: Efficiency histograms for each prefix/particle.
    """
    df_new = df.copy()

    # We separate the dataframe into two parts: one where all the events have
    # PID indices (events inside the PID binning) and those that don't.
    # Efficiency is added only for events inside the PID binning.
    df_nan = df_new[df_new.isna().any(axis=1)]
    df_new.dropna(inplace=True)

    df_new["PIDCalibEff"] = 1
    df_new["PIDCalibErr2"] = 0

    for prefix in prefixes:
        efficiency_table = eff_hists[prefix]["eff"].view().flatten()  # type: ignore
        error_table = (
            create_error_histogram(eff_hists[prefix]).view().flatten()  # type: ignore
        )
        np.nan_to_num(efficiency_table, False)  # Replicate old PIDCalib's behavior

        # Assign efficiencies by taking the efficiency value from the relevant bin
        df_new[f"{prefix}_PIDCalibEff"] = np.take(
            efficiency_table, df_new[f"{prefix}_PIDCalibBin"]
        )
        # Assign errors by taking the error value from the relevant bin
        df_new[f"{prefix}_PIDCalibErr"] = np.take(
            error_table, df_new[f"{prefix}_PIDCalibBin"]
        )
        df_new["PIDCalibEff"] = df_new["PIDCalibEff"] * df_new[f"{prefix}_PIDCalibEff"]
        df_new["PIDCalibErr2"] += df_new[f"{prefix}_PIDCalibErr"] ** 2

    df_new["PIDCalibErr"] = np.sqrt(df_new["PIDCalibErr2"])  # type: ignore
    df_new.drop(columns=["PIDCalibErr2"], inplace=True)

    df_new = pd.concat([df_new, df_nan]).sort_index()
    log.debug("Particle efficiencies assigned")

    num_outside_range = len(df_nan.index)
    num_outside_range_frac = len(df_nan.index) / len(df_new.index)
    log.debug(
        f"Events out of range: {num_outside_range} ({num_outside_range_frac:.2%})"
    )
    return df_new


def create_hist_filename(
    sample: str, magnet: str, particle: str, pid_cut: str, bin_vars: List[str]
) -> str:
    """Return effhists filename corresponding to parameters.

    Args:
        sample: Data sample name (Turbo18, etc.).
        magnet: Magnet polarity (up, down).
        particle: Particle type (K, Pi, etc.).
        pid_cut: Simplified user-level cut, e.g., "DLLK < 4".
        bin_vars: Variables used for binning.
    """
    whitespace = re.compile(r"\s+")
    cut = re.sub(whitespace, "", pid_cut)

    return f"effhists_{sample}_{magnet}_{particle}_{cut}_{'-'.join(bin_vars)}.pkl"


def binomial_uncertainty(
    num_pass: Union[float, np.ndarray],
    num_total: Union[float, np.ndarray],
    err_pass_sq: Union[float, np.ndarray],
    err_tot_sq: Union[float, np.ndarray],
) -> float:
    """Return the uncertainty of binomial experiments.

    The parameters can be either floats or numpy arrays.

    The uncertainty is calculated the way ROOT does it in TH1::Divide() when
    binomial errors are specified. This approach has known problems when
    num_pass == num_total or 0. We use this approach to ensure compatibility
    with the original PIDCalib, and because these edge-cases are unlikely to
    occur.

    Args:
        num_pass: Number of passing events.
        num_total: Total number of events.
        err_pass_sq: Squared uncertainty on the number of passing events (sum
            of squares of the event weights).
        err_tot_sq: Squared uncertainty on the number of total events (sum of
            squares of the event weights).
    """
    prob = num_pass / num_total
    prob_sq = prob ** 2  # type: ignore
    num_total_sq = num_total ** 2  # type: ignore
    return np.sqrt(  # type: ignore
        abs(((1 - 2 * prob) * err_pass_sq + err_tot_sq * prob_sq) / num_total_sq)
    )


def create_error_histogram(eff_hists: Dict[str, bh.Histogram]) -> bh.Histogram:
    uncertainty = binomial_uncertainty(
        eff_hists["passing"].view(flow=False),  # type: ignore
        eff_hists["total"].view(flow=False),  # type: ignore
        eff_hists["passing_sumw2"].view(flow=False),  # type: ignore
        eff_hists["total_sumw2"].view(flow=False),  # type: ignore
    )

    err_histo = eff_hists["passing"].copy()
    err_histo[...] = uncertainty
    return err_histo
