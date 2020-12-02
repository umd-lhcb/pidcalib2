import argparse
import pathlib
import re

import pandas
from logzero import logger as log

from . import utils


def decode_arguments():
    """Decode CLI arguments."""
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "-y",
        "--year",
        help="year of data-taking",
        required=True,
        type=int,
        choices=[2018, 2017, 2016, 2015],
    )
    parser.add_argument(
        "-m", "--magnet", help="magnet polarity", required=True, choices=["up", "down"],
    )
    parser.add_argument(
        "-p",
        "--particle",
        help="particle type",
        required=True,
        choices=["K", "pi", "p", "mu", "e"],
    )
    parser.add_argument(
        "-c",
        "--pid-cut",
        help="PID cut string, e.g., 'DLLK<4.0' (-c can be used multiple times for multiple cuts)",
        action="append",
        dest="pid_cuts",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--bin-var",
        help="binning variable (-v can be used multiple times for multiple variables)",
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
    args = parser.parse_args()
    return args


def log_config(config: dict) -> None:
    log.info("Running PIDCalib2 make_perf_hists with the following config:")
    longest_key = len(max(config, key=len))
    log.info("=" * longest_key)
    for entry in config:
        log.info(f"{entry:{longest_key}}: {config[entry]}")
    log.info("=" * longest_key)


def make_perf_hists(config: dict) -> None:
    """Create sWeighted PID calibration histograms and save them to disk.

    Bla bla.
    #TODO Finish docstring
    """
    # log_format = '%(color)s[%(levelname)1.1s %(module)s:%(lineno)d]%(end_color)s %(message)s'
    # formatter = logzero.LogFormatter(fmt=log_format)
    # logzero.setup_default_logger(formatter=formatter)

    # TODO Setup log file
    # pattern = re.compile(r"\s+")
    # config["pid_cuts"] = [re.sub(pattern, '', pid_cut) for pid_cut in config["pid_cuts"]]
    log_config(config)

    pathlib.Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

    # TODO Change name/review if necessary
    calib_par_name = "probe"
    branch_names = utils.get_relevant_branch_names(
        calib_par_name, config["pid_cuts"], config["bin_vars"]
    )
    log.info(f"Branches to be read: {branch_names}")

    eos_paths = utils.get_eos_paths(config["year"], config["magnet"])
    log.info(f"{len(eos_paths)} calibration files from EOS will be processed")

    df_tot = None
    if config["local_dataframe"] is not None:
        df_tot = pandas.read_pickle(config["local_dataframe"])
    else:
        df_tot = utils.extract_branches_to_dataframe(
            eos_paths, config["particle"], branch_names
        )
    # df_tot.to_pickle("local_dataframe.pkl")

    pid_cuts = utils.translate_pid_cuts_to_branch_cuts(
        calib_par_name, config["pid_cuts"]
    )


if __name__ == "__main__":
    config = vars(decode_arguments())
    make_perf_hists(config)
