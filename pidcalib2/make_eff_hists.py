import argparse
import pathlib
import re

import pickle
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
        choices=[2018, 2017, 2016, 2015, 2014, 2013],
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
        help=(
            "PID cut string, e.g., 'DLLK<4.0' (-c can be used multiple times for "
            "multiple cuts)"
        ),
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
    parser.add_argument(
        "-f",
        "--file-list",
        help="(debug) read calibration file paths from a text file",
    )
    parser.add_argument(
        "-n", "--max-files", type=int, help="(debug) a max number of files to read",
    )
    args = parser.parse_args()
    return args


def make_eff_hists(config: dict) -> None:
    """Create sWeighted PID calibration histograms and save them to disk.

    Calibration samples from EOS are read and relevant branches extracted to
    a DataFrame. Each PID cut is applied to the DataFrame in turn and the
    results are histogrammed (each event with its associated sWeight).
    Particle type and binning variables are used to select an appropriate
    predefined binning. The efficiency histograms are saved to a requested
    output directory.
    """
    # log_format = '%(color)s[%(levelname)1.1s %(module)s:%(lineno)d]%(end_color)s %(message)s'  # noqa
    # formatter = logzero.LogFormatter(fmt=log_format)
    # logzero.setup_default_logger(formatter=formatter)

    # TODO Setup log file
    # Remove whitespace from PID cuts to avoid DLLK < 4 != DLLK<4
    pattern = re.compile(r"\s+")
    config["pid_cuts"] = [
        re.sub(pattern, "", pid_cut) for pid_cut in config["pid_cuts"]
    ]

    log.info("Running PIDCalib2 make_perf_hists with the following config:")
    utils.log_config(config)

    output_dir = pathlib.Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    calibration_prefix = "probe"
    branch_names = utils.get_relevant_branch_names(
        calibration_prefix, config["pid_cuts"], config["bin_vars"]
    )
    log.info(f"Branches to be read: {branch_names}")

    df_total = None
    if config["local_dataframe"]:
        df_total = utils.dataframe_from_local_file(
            config["local_dataframe"], list(branch_names)
        )
    else:
        if config["file_list"]:
            with open(config["file_list"]) as f:
                eos_paths = f.read().splitlines()
                print(eos_paths)
        else:
            eos_paths = utils.get_eos_paths(
                config["year"], config["magnet"], config["max_files"]
            )
        tree_paths = utils.get_tree_paths(config["particle"], config["year"])

        log.info(f"{len(eos_paths)} calibration files from EOS will be processed")
        df_total = utils.calib_root_to_dataframe(eos_paths, tree_paths, branch_names)
        # df_total.to_pickle("local_dataframe.pkl")
        # df_total.to_csv("local_dataframe.csv")

    eff_hists = utils.create_eff_histograms(
        df_total, config["particle"], config["pid_cuts"], config["bin_vars"]
    )

    for name in eff_hists:
        if name.startswith("eff_"):
            cut = name.replace("eff_", "")
            hist_filename = (
                f"effhists_"
                f"{config['year']}_"
                f"{config['magnet']}_"
                f"{config['particle']}_"
                f"{cut}_"
                f"{'_'.join(config['bin_vars'])}"
                ".pkl"
            )
            eff_hist_path = output_dir / hist_filename
            with open(eff_hist_path, "wb") as f:
                pickle.dump(eff_hists[f"eff_{cut}"], f)
                pickle.dump(eff_hists[f"passing_{cut}"], f)
                pickle.dump(eff_hists["total"], f)

            log.info(f"Efficiency histograms saved to '{eff_hist_path}'")


if __name__ == "__main__":
    config = vars(decode_arguments())
    make_eff_hists(config)
