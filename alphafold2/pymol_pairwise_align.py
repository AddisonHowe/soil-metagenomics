"""Pymol pairwise alignment script

Usage: 
    python pymol_pairwise_align.py [-f <pdb_file1.pdb> <pdb_file2.pdb> ...] 
                                   [-n <name1> <name2> ...]
                                   [-i <pdb_id1> <pdb_id2> ...] 
                                   [-o <outdir>]
Inputs:
    [-f --pdb_files] : pdb filepaths to align.
    [-n --names] : names corresponding to each file
    [-i --pdb_ids] : PDB IDs to fetch from the online database.
    [-o --outdir] : output directory

Outputs:
    Creates an output directory <outdir> storing results of each pairwise 
    alignment comparison.

"""

import sys, os
import argparse
import tqdm
import pymol
from pymol import cmd

from pymol_helpers import Item, load_pymol_item, get_printv


FETCHED_PDBID_CACHE = "./out/pymol/fetched_pdbids"


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--pdb_files", type=str, nargs="*", 
                        help="pdb files to load")
    parser.add_argument("-n", "--names", type=str, nargs="*", 
                        default=None, help="names")
    parser.add_argument("-i", "--pdb_ids", type=str, nargs="*", 
                        help="pdb IDs to fetch from PDB database")
    parser.add_argument("-o", "--outdir", type=str, 
                        help="output directory")
    parser.add_argument("-v", "--verbosity", type=int, default=1)
    parser.add_argument("--pbar", action="store_true")
    return parser.parse_args(args)


def run_pymol_align_pairwise(
        item_list, 
        names,
        verbosity=1,
        use_pbar=False,
):
    printv = get_printv(verbosity)
    nitems = len(item_list)
    os.makedirs(FETCHED_PDBID_CACHE, exist_ok=True)
    
    pymol.finish_launching(["pymol", "-cq"])  # quiet + no GUI

    alignments = {}
    for i in tqdm.trange(nitems-1, desc="Item i", disable=not use_pbar):
        item1 = item_list[i]
        name1 = names[i]
        for j in tqdm.trange(i, nitems, desc="Item j", disable=not use_pbar, leave=False):
            if i == j:
                continue
            item2 = item_list[j]
            name2 = names[j]

            printv(f"Aligning items:", pbar=use_pbar)
            printv(f"\t{item1.value} ({item1.kind})", pbar=use_pbar)
            printv(f"\t{item2.value} ({item2.kind})", pbar=use_pbar)

            alignment = pymol_align_structures(
                item1, item2, 
            )
            alignments[(name1, name2)] = alignment
            printv((f"\tRMSD: {alignment:.3f}"), pbar=use_pbar)

    # Finish PyMOL session
    cmd.quit()
    
    return alignments



def pymol_align_structures(
        item1, item2, 
        extra_cmds1=[], 
        extra_cmds2=[],
):
    # Load the two structures
    for i, item in enumerate([item1, item2]):
        load_pymol_item(f"struct{i+1}", item, fetch_path=FETCHED_PDBID_CACHE)
    # Align struct1 to struct2 and get RMSD
    alignment_result = cmd.align('struct1', 'struct2')
    rmsd = alignment_result[0]
    # Delete the references
    cmd.delete("struct1")
    cmd.delete("struct2")
    return rmsd


def save_alignment_data(alignments, outdir):
    print(f"Saving alignment results to {outdir}")
    with open(f"{outdir}/pairwise_alignments.tsv", "w") as f:
        for pair, alignment in alignments.items():
            f.write(f"{pair[0]}\t{pair[1]}\t{alignment}\n")   
    return


def main(args):
    print(f"Handled args: {args}")
    pdb_files = args.pdb_files if args.pdb_files else []
    pdb_ids = args.pdb_ids if args.pdb_ids else []
    names = args.names if args.names else pdb_files
    outdir = args.outdir
    verbosity = args.verbosity
    use_pbar = args.pbar

    print("Loading files:", pdb_files)
    print("Fetching ids:", pdb_ids)
    print("Saving output to", outdir)

    if len(pdb_files) != len(names):
        raise RuntimeError("Length of names should match length of files")

    loaded_structures = []
    structure_names = []
    
    for i, fpath in enumerate(pdb_files):
        struct = Item(fpath, "file")
        loaded_structures.append(struct)
        structure_names.append(names[i])
    
    for id in pdb_ids:
        struct = Item(id, "id")
        loaded_structures.append(struct)
        structure_names.append(id)

    alignments = run_pymol_align_pairwise(
        loaded_structures,
        names=structure_names,
        verbosity=verbosity,
        use_pbar=use_pbar,
    )

    # Save alignments
    os.makedirs(outdir, exist_ok=True)
    save_alignment_data(alignments, outdir=outdir)

    return


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
