import argparse
import ast
from typing import Dict
import time

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
        "-v",
        "--bin-vars",
        help=(
            "dictionary of binning variables (keys) and their associated names in "
            "reference sample (values) e.g. \"{'P': 'P', 'ETA' : 'Eta'}\""
        ),
        default="{'P' : 'P', 'ETA' : 'ETA', 'nTracks' : 'nTracks'}",
        dest="bin_vars",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="pidcalib_output",
        help="directory where to save output files",
    )
    parser.add_argument(
        "-f", "--ref-file", help="reference sample file", required=True,
    )
    parser.add_argument(
        "-t", "--ref-tree", help="reference sample tree name", default="DecayTree",
    )
    parser.add_argument(
        "-p",
        "--ref-pars",
        help=(
            "JSON dictionary of particles from reference sample to apply "
            "cuts to, where the keys represent the particle branch names, "
            "and the values passed are a list containing particle type and "
            "PID cut e.g. \"{'D0_K' : ['K','DLLK>4.0'], 'D0_Pi' : ['Pi','DLLK<4.0']}\""
        ),
        required=True,
    )
    args = parser.parse_args()
    return args


def ref_calib(config: Dict) -> float:
    """TODO: Write docstring

    Args:
        config ([type]): [description]
    """
    log.info("Running PIDCalib2 ref_calib with the following config:")
    utils.log_config(config)

    try:
        bin_vars = ast.literal_eval(config["bin_vars"])
    except SyntaxError:
        log.error("The --bin-vars string is not a valid Python dict")
        raise

    try:
        ref_pars = ast.literal_eval(config["ref_pars"])
    except SyntaxError:
        log.error("The --ref-pars string is not valid Python dict")
        raise

    ref_branches = utils.get_reference_branch_names(ref_pars, bin_vars)

    log.info(f"Loading reference sample '{config['ref_file']}' ...")
    df_ref = utils.root_to_dataframe(
        config["ref_file"], config["ref_tree"], ref_branches
    )
    log.debug(
        f"Reference sample '{config['ref_file']}' with {len(df_ref.index)} events loaded"  # noqa
    )

    # TODO: Rename output_dir to hist_dir (?)
    eff_histos = utils.get_calib_hists(
        config["output_dir"], config["year"], config["magnet"], ref_pars, bin_vars
    )

    start = time.perf_counter()
    df_ref = utils.add_bin_indices(df_ref, ref_pars, bin_vars, eff_histos)
    df_ref = utils.add_efficiencies(df_ref, ref_pars, eff_histos)
    end = time.perf_counter()
    log.debug(f"Efficiency calculation took {end-start:.2f}s")

    # Calculate average of the per-event effs
    # Use only data with valid eff values (those events falling inside the
    # calibration hist)
    avg_eff = df_ref["eff"].dropna().mean()
    log.info(f"Average per-event PID efficiency: {avg_eff:.2%}")
    return avg_eff  # type: ignore


if __name__ == "__main__":
    config = vars(decode_arguments())
    ref_calib(config)
