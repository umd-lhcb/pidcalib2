# -*- coding: utf-8 -*-
"""Copy tree from a file and set it as friend of a tree in another file.

PIDCalib2.ref_calib creates a file with a tree where the efficiency for each
track and event from the user's reference sample is stored. This PID
efficiency information has to be added to the existing event information.
This script is used to copy the efficiency tree to the user's reference file
and set it as a friend of the reference tree. This way one can use the
reference tree as if it had the branches from the efficiency tree.

Example:
    $ python -m pidcalib2.merge_trees eff.root eff_tree ref.root DecayTree
"""

import sys

import ROOT
from logzero import logger as log


def copy_tree_and_set_as_friend(
    source_file_name: str,
    source_tree_name: str,
    dest_file_name: str,
    dest_tree_name: str,
) -> None:
    """Copy tree from a file and set it as friend of a tree in another file.

    Args:
        source_file_name: File from which to copy the tree.
        source_tree_name: Name of the tree to be copied.
        dest_file_name: File to which to copy the tree.
        dest_tree_name: Name of the tree on which to call AddFriend.
    """
    source_file = ROOT.TFile(source_file_name, "read")  # type: ignore
    source_tree = source_file.Get(source_tree_name)
    log.info(f"Reading {source_tree_name} from {source_file_name}")
    if not source_tree:
        print(f"ERROR: '{source_tree_name}' not found in '{source_file_name}'.")
        sys.exit(1)

    dest_file = ROOT.TFile(dest_file_name, "update")  # type: ignore
    dest_file.cd()

    log.info(f"Reading {dest_tree_name} from {dest_file_name}")
    dest_tree = dest_file.Get(dest_tree_name)
    if not dest_tree:
        print(f"ERROR: '{dest_tree_name}' not found in '{dest_file_name}'.")
        sys.exit(2)

    dest_tree.AddFriend(source_tree)

    log.info(f"Copying {source_tree_name} to {dest_file_name}")
    source_tree.CloneTree().Write()
    # Avoid ROOT adding a new tree with an incremented cycle number
    dest_tree.Write("", ROOT.TObject.kOverwrite)  # type: ignore

    dest_file.Close()


def main():
    """The main function."""
    efficiency_file = sys.argv[1]
    efficiency_tree = sys.argv[2]
    reference_file = sys.argv[3]
    reference_tree = sys.argv[4]

    copy_tree_and_set_as_friend(
        efficiency_file, efficiency_tree, reference_file, reference_tree
    )


if __name__ == "__main__":
    main()
