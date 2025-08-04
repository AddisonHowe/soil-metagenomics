"""Pairwise alignment script

Usage: 
    python pairwise_align.py [-f <pdb_file1.pdb> <pdb_file2.pdb> ...] 
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

from Bio.Align import PairwiseAligner

from pymol_helpers import get_printv

from mgsa.io import load_pdb_structure
from mgsa.structures import get_ca_atoms, align_atoms_in_structures
from mgsa.structures import align_atoms_in_structures, extract_residue_sequence


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


def run_align_pairwise(
        structure_list, 
        names,
        verbosity=1,
        use_pbar=False,
):
    printv = get_printv(verbosity)
    nitems = len(structure_list)
    os.makedirs(FETCHED_PDBID_CACHE, exist_ok=True)
    
    printv("Subsetting C-alpha atoms", pbar=use_pbar)
    ca_atoms = []
    for i in tqdm.trange(nitems, desc="Item", disable=not use_pbar):
        struct = structure_list[i]
        ca_atoms.append(get_ca_atoms(struct))

    # seq_alignments = {}

    aligner = PairwiseAligner(scoring="blastp")
    aligner.mode = "global"

    struct_alignments = {}
    for i in tqdm.trange(nitems-1, desc="Item i", disable=not use_pbar):
        struct1 = structure_list[i]
        name1 = names[i]
        residues1, seq1 = extract_residue_sequence(struct1)
        for j in tqdm.trange(i, nitems, desc="Item j", disable=not use_pbar, leave=False):
            if i == j:
                continue
            
            struct2 = structure_list[j]
            name2 = names[j]
            residues2, seq2 = extract_residue_sequence(struct2)

            # Perform sequence alignment
            printv(f"Aligning sequences {i}, {j}", pbar=use_pbar, v=2)
            printv(f"\tSeq {i}: {seq1[0:20]}...", pbar=use_pbar, v=3)
            printv(f"\tSeq {j}: {seq2[0:20]}...", pbar=use_pbar, v=3)
            alignments = align_sequences(seq1, seq2, aligner=aligner)
            alignment = alignments[0]
            aligned_seq1 = alignment.aligned[0]
            aligned_seq2 = alignment.aligned[1]
            printv(str(alignment[:,0:20]).replace("\n", "...\n"), pbar=use_pbar, v=2)
            printv(f"aligned seq i:\n{aligned_seq1}", pbar=use_pbar, v=3)
            printv(f"aligned seq j:\n{aligned_seq2}", pbar=use_pbar, v=3)

            # Align the C-alpha atoms
            atoms1, atoms2 = get_aligned_atoms(
                aligned_seq1, aligned_seq2, residues1, residues2
            )
            
            # Perform structural alignment
            struct_alignments[(name1, name2)] = align_atoms_in_structures(
                struct1, struct2, atoms1, atoms2
            )
    
    return struct_alignments


def get_aligned_atoms(aligned_seq1, aligned_seq2, residues1, residues2):
    atoms1 = []
    atoms2 = []
    i, j = 0, 0
    for a1, a2 in zip(aligned_seq1, aligned_seq2):
        if a1 != '-' and a2 != '-':
            if "CA" in residues1[i] and "CA" in residues2[j]:
                atoms1.append(residues1[i]["CA"])
                atoms2.append(residues2[j]["CA"])
        if a1 != '-':
            i += 1
        if a2 != '-':
            j += 1
    return atoms1, atoms2


def align_sequences(seq1, seq2, aligner=None):
    if aligner is None:
        aligner = PairwiseAligner(scoring="blastp")
    alignment = aligner.align(seq1, seq2)
    return alignment


def save_alignment_data(alignments, outdir):
    print(f"Saving alignment results to {outdir}")
    print(alignments)
    with open(f"{outdir}/pairwise_alignments.tsv", "w") as f:
        for pair, alignment in alignments.items():
            f.write(f"{pair[0]}\t{pair[1]}\t{alignment}\n")   
    return


def main(args):
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
        struct = load_pdb_structure(fpath, f"{names[i]}")
        loaded_structures.append(struct)
        structure_names.append(names[i])
    
    # for id in pdb_ids:
    #     struct = Item(id, "id")
    #     loaded_structures.append(struct)
    #     structure_names.append(id)

    alignments = run_align_pairwise(
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
