""" pymol pairwise alignment script

Usage: 
    python pymol_pairwise_align.py [-f <pdb_file1.pdb> <pdb_file2.pdb> ...] 
                                   [-i <pdb_id1> <pdb_id2> ...] 
                                   [-o <outdir>]
Inputs:
    [-f --pdb_files] : pdb filepaths to align.
    [-i --pdb_ids] : PDB IDs to fetch from the online database.
    [-o --outdir] : output directory

Outputs:
    Creates an output directory <outdir> storing results of each pairwise 
    comparison.

"""

import sys, os
import argparse
import pymol
from pymol import cmd


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--pdb_files", type=str, nargs="*", 
                        help="pdb files to load")
    parser.add_argument("-i", "--pdb_ids", type=str, nargs="*", 
                        help="pdb IDs to fetch from PDB database")
    parser.add_argument("-o", "--outdir", type=str, 
                        help="output directory")
    return parser.parse_args(args)


def run_pymol_align_pairwise(item_list):
    nitems = len(item_list)
    alignments = [[None for _ in range(nitems)] for _ in range(nitems)]
    
    pymol.finish_launching(["pymol", "-cq"])  # quiet + no GUI

    for i in range(nitems):
        item1 = item_list[i]
        for j in range(nitems):
            if i == j:
                continue
            item2 = item_list[j]

            alignment = pymol_align_structures(item1, item2)
            alignments[i][j] = alignment

    # Finish PyMOL session
    cmd.quit()
    
    return alignments


class Item:

    def __init__(self, value, kind):
        self.value = value
        self.kind = kind


def load_pymol_item(ref, item):
    if item.kind == "file":
        return cmd.load(item.value, ref)
    elif item.kind == "id":
        return cmd.fetch(item.value, ref)
    else:
        raise RuntimeError(f"Unknown kind {item.kind}")


def pymol_align_structures(
        item1, item2, 
        extra_cmds1=[], 
        extra_cmds2=[],
        verbosity=1,
):
    
    if verbosity:
        print(f"Aligning items:")
        print(f"\t{item1.value} ({item1.kind})")
        print(f"\t{item2.value} ({item2.kind})")

    # Load the two structures
    for i, item in enumerate([item1, item2]):
        load_pymol_item(f"struct{i+1}", item)
    
    # Align struct1 to struct2 and get RMSD
    alignment_result = cmd.align('struct1', 'struct2')
    rmsd = alignment_result[0]
    
    if verbosity:
        print(f"\tRMSD: {rmsd:.3f}")

    # Delete the references
    cmd.delete("struct1")
    cmd.delete("struct2")

    return rmsd


def save_alignments(alignments, outdir):
    print(f"Saving alignment results to {outdir}")
    # TODO: save outputs
    return


def main(args):
    print(f"Handled args: {args}")
    pdb_files = args.pdb_files
    pdb_ids = args.pdb_ids
    outdir = args.outdir

    print("Loading files:", pdb_files)
    print("Fetching ids:", pdb_ids)
    print("Saving output to", outdir)

    loaded_structures = []
    
    for fpath in pdb_files:
        struct = Item(fpath, "file")
        loaded_structures.append(struct)
    
    for id in pdb_ids:
        struct = Item(id, "id")
        loaded_structures.append(struct)

    alignments = run_pymol_align_pairwise(loaded_structures)

    # Save alignments
    os.makedirs(outdir, exist_ok=True)
    save_alignments(alignments, outdir=outdir)

    return


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)


# Load arguments
# load $1, protein1
# load $2, protein2
# # Align protein2 to protein1
# align protein2, protein1

# # Get alignment score (RMSD)
# set rmsd = cmd.align("protein2", "protein1")[0]

# # Print to terminal
# print("RMSD: ", rmsd)
