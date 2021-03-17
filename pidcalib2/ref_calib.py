import argparse
import ast
import pathlib
import time
from typing import Dict

from logzero import logger as log

from . import merge_trees, pid_data, utils


def decode_arguments():
    """Decode CLI arguments."""
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument(
        "-s",
        "--sample",
        help="calibration sample (Turbo18, Electron16, ...)",
        required=True,
        type=str,
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
    parser.add_argument(
        "-d", "--dry-run", help="do not update the reference file", action="store_true"
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

    ref_branches = pid_data.get_reference_branch_names(ref_pars, bin_vars)

    log.info(f"Loading reference sample '{config['ref_file']}' ...")
    df_ref = pid_data.root_to_dataframe(
        config["ref_file"], config["ref_tree"], ref_branches
    )
    log.debug(
        f"Reference sample '{config['ref_file']}' with {len(df_ref.index)} events loaded"  # noqa
    )

    # TODO: Rename output_dir to hist_dir (?)
    eff_histos = pid_data.get_calib_hists(
        config["output_dir"], config["sample"], config["magnet"], ref_pars, bin_vars
    )

    start = time.perf_counter()
    df_ref = utils.add_bin_indices(df_ref, ref_pars, bin_vars, eff_histos)
    df_ref = utils.add_efficiencies(df_ref, ref_pars, eff_histos)
    end = time.perf_counter()
    log.debug(f"Efficiency calculation took {end-start:.2f}s")

    # Calculate average of the per-event effs
    # Use only data with valid eff values (those events falling inside the
    # calibration hist)
    avg_eff = df_ref["PID_eff"].dropna().mean()
    log.info(f"Average per-event PID efficiency: {avg_eff:.2%}")

    output_path = pathlib.Path(config["output_dir"])
    ref_path = pathlib.Path(config["ref_file"])

    eff_path = output_path / ref_path.name.replace(".root", "_PID_eff.root")

    pid_data.save_dataframe_as_root(
        df_ref[[key for key in df_ref.keys() if key not in ref_branches]],
        "PID_eff_tree",
        str(eff_path),
    )

    if not config["dry_run"]:
        merge_trees.copy_tree_and_set_as_friend(
            str(eff_path), "PID_eff_tree", config["ref_file"], config["ref_tree"]
        )
    else:
        log.warning("This is a dry run, the reference file was not updated")

    return avg_eff  # type: ignore


if __name__ == "__main__":
    config = vars(decode_arguments())
    ref_calib(config)
