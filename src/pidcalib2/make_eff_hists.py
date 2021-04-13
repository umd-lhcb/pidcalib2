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

"""Module to make LHCb PID efficiency histograms.

This module creates histograms that can be used to estimate the PID
efficiency of a user's sample.

Examples:
    Create a single efficiency histogram for a single PID cut::

        $ python -m pidcalib2.make_eff_hists --sample=Turbo18 --magnet=up \
            --particle=Pi --pid-cut="DLLK > 4" --bin-var=P --bin-var=ETA \
            --bin-var=nSPDHits --output-dir=pidcalib_output

    Create multiple histograms in one run (most of the time is spent reading
    in data, so specifying multiple cuts is much faster than running
    make_eff_hists sequentially)::

        $ python -m pidcalib2.make_eff_hists --sample=Turbo16 --magnet=up \
            --particle=Pi --pid-cut="DLLK > 0" --pid-cut="DLLK > 4" \
            --pid-cut="DLLK > 6" --bin-var=P --bin-var=ETA \
            --bin-var=nSPDHits --output-dir=pidcalib_output
"""

import argparse
import logging
import pathlib
import pickle
import re
import sys

import logzero
from logzero import logger as log

from . import pid_data, utils

try:
    from .version import version  # type: ignore
except ImportError:
    version = "N/A"


def decode_arguments(args):
    """Decode CLI arguments."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--sample",
        help="calibration sample (Turbo18, Electron16, ...)",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-m",
        "--magnet",
        help="magnet polarity",
        required=True,
        choices=["up", "down"],
    )
    parser.add_argument(
        "-p",
        "--particle",
        help="particle type",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--pid-cut",
        help=(
            "PID cut string, e.g., 'DLLK < 4.0' (-i can be used multiple times for "
            "multiple cuts). Supported operators are <, >, ==, and !="
        ),
        action="append",
        dest="pid_cuts",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--cut",
        help=(
            "arbitrary cut string, e.g., 'Dst_IPCHI2 < 10.0' (-c can be used multiple "
            "times for multiple cuts). Supported operators are <, >, ==, and !="
        ),
        action="append",
        dest="cuts",
    )
    parser.add_argument(
        "-b",
        "--bin-var",
        help="binning variable (-b can be used multiple times for multiple variables)",
        action="append",
        dest="bin_vars",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="pidcalib_output",
        help="directory where to save output files",
    )
    parser.add_argument(
        "-l",
        "--local-dataframe",
        help="(debug) read a calibration DataFrame from file",
    )
    parser.add_argument(
        "-f",
        "--file-list",
        help="(debug) read calibration file paths from a text file",
    )
    parser.add_argument(
        "-a",
        "--samples-file",
        help="(debug) read calibration sample lists from a custom file",
    )
    parser.add_argument(
        "-n",
        "--max-files",
        type=int,
        help="(debug) a max number of files to read",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="(debug) increase verbosity",
    )
    parser.add_argument("-V", "--version", action="version", version=version)
    parsed_args = parser.parse_args(args)
    return parsed_args


def make_eff_hists(config: dict) -> None:
    """Create sWeighted PID calibration histograms and save them to disk.

    Calibration samples from EOS are read and relevant branches extracted to
    a DataFrame. Each PID cut is applied to the DataFrame in turn and the
    results are histogrammed (each event with its associated sWeight).
    Particle type and binning variables are used to select an appropriate
    predefined binning. The efficiency histograms are saved to a requested
    output directory.

    Args:
        config: A configuration dictionary. See decode_arguments(args) for
            details.
    """
    if config["verbose"]:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)

    pattern = re.compile(r"\s+")
    config["pid_cuts"] = [
        re.sub(pattern, "", pid_cut) for pid_cut in config["pid_cuts"]
    ]

    log.info("Running PIDCalib2 make_eff_hists with the following config:")
    utils.log_config(config)

    output_dir = pathlib.Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    calibration_prefix = "probe"
    branch_names = pid_data.get_relevant_branch_names(
        calibration_prefix, config["pid_cuts"], config["bin_vars"], config["cuts"]
    )
    log.info(f"Branches to be read: {branch_names}")

    df_total = None
    if config["local_dataframe"]:
        df_total = pid_data.dataframe_from_local_file(
            config["local_dataframe"], list(branch_names)
        )
    else:
        if config["file_list"]:
            with open(config["file_list"]) as f_list:
                eos_paths = f_list.read().splitlines()
        else:
            eos_paths = pid_data.get_file_list(
                config["sample"],
                config["magnet"],
                config["particle"],
                config["samples_file"],
                config["max_files"],
            )
        tree_paths = pid_data.get_tree_paths(config["particle"], config["sample"])

        log.info(f"{len(eos_paths)} calibration files from EOS will be processed")
        for path in eos_paths:
            log.debug(f"  {path}")
        df_total = pid_data.calib_root_to_dataframe(eos_paths, tree_paths, branch_names)
        # df_total.to_pickle("local_dataframe.pkl")
        # df_total.to_csv("local_dataframe.csv")

    if config["cuts"]:
        cut_string = " & ".join(config["cuts"])
        log.debug(f"Applying generic cuts: '{cut_string}'")
        num_before = df_total.shape[0]
        df_total.query(cut_string, inplace=True)
        num_after = df_total.shape[0]
        log.debug(
            f"{num_after}/{num_before} ({num_after/num_before:.1%}) events "
            "passed generic cuts"
        )

    eff_hists = utils.create_eff_histograms(
        df_total, config["particle"], config["pid_cuts"], config["bin_vars"]
    )

    for name in eff_hists:
        if name.startswith("eff_"):
            cut = name.replace("eff_", "")
            hist_filename = utils.create_hist_filename(
                config["sample"],
                config["magnet"],
                config["particle"],
                cut,
                config["bin_vars"],
            )
            eff_hist_path = output_dir / hist_filename
            with open(eff_hist_path, "wb") as f:
                pickle.dump(eff_hists[f"eff_{cut}"], f)
                pickle.dump(eff_hists[f"passing_{cut}"], f)
                pickle.dump(eff_hists["total"], f)
                pickle.dump(eff_hists[f"passing_sumw2_{cut}"], f)
                pickle.dump(eff_hists["total_sumw2"], f)

            log.info(f"Efficiency histograms saved to '{eff_hist_path}'")


def main():
    config = vars(decode_arguments(sys.argv[1:]))
    make_eff_hists(config)


if __name__ == "__main__":
    main()
