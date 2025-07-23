"""Structure analysis script

Usage: 
    python structure_analysis.py [-f <pdb_file1.pdb> <pdb_file2.pdb> ...] 
                                 [-n <name1> <name2> ...]
                                 [-o <outdir>] 
                                 [--outfilename <outfilename>]
Inputs:
    [-f --pdb_files] : pdb filepaths
    [-n --names] : names corresponding to each file
    [-o --outdir] : output directory
    [--outfname] : name of results file (Default is structure_analysis.tsv)

Outputs:
    Creates a tsv results file in the output directory, each row corresponding 
    to one of the input files.

"""

import sys, os
import argparse
import tqdm
import pymol
from pymol import cmd

from pymol_helpers import Item, load_pymol_item, get_printv

from mgsa.io import load_pdb_structure
import mgsa.structures as mgsastruct


FETCHED_PDBID_CACHE = "./out/pymol/fetched_pdbids"

KEY_ORDER = [
    "id",
    "rg",
    "helix",
    "sheet",
    "coil",
    "sasa",
    "instability",
    "isoelectric_point",
    "length",
]

KEYS = {k: k for k in KEY_ORDER}  # Default header
# Custom header names
KEYS["rg"] = "radius of gyration"
KEYS["isoelectric_point"] = "isoelectric point"
KEYS["helix"] = "prop helix"
KEYS["sheet"] = "prop sheet"
KEYS["coil"] = "prop coil"


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--pdb_files", type=str, nargs="*", 
                        help="pdb files to load")
    parser.add_argument("-n", "--names", type=str, nargs="*", 
                        default=None, help="names")
    parser.add_argument("-o", "--outdir", type=str, 
                        help="output directory")
    parser.add_argument("--outfname", type=str, default="structure_analysis",
                        help="name of output file")
    parser.add_argument("-v", "--verbosity", type=int, default=1)
    parser.add_argument("--pbar", action="store_true")
    return parser.parse_args(args)


def analyze_protein(struct) -> dict:
    results = {}
    # ID value
    results['id'] = struct.id
    # Quantitative structural measures
    results['rg'] = mgsastruct.rgyrate(struct)
    ss_props = mgsastruct.get_secondary_structure_proportions(struct)
    results["helix"] = ss_props["helix"]
    results["sheet"] = ss_props["sheet"]
    results["coil"] = ss_props["coil"]
    results["sasa"] = mgsastruct.compute_sasa(struct, level="S")
    # Direct sequence analysis
    analyzed_seq = mgsastruct.analyze_sequence(struct)
    results['instability'] = analyzed_seq.instability_index()
    results['isoelectric_point'] = analyzed_seq.isoelectric_point()
    results['length'] = analyzed_seq.length
    return results


def write_result_line(result, f):
    vals = [str(result[k]) for k in KEY_ORDER]
    line = "\t".join(vals)
    f.write(line + "\n")
    return line


def write_results(results, outfpath):
    with open(outfpath, "w") as f:
        f.write("\t".join([KEYS[k] for k in KEY_ORDER]) + "\n")
    for res in results:
        with open(outfpath, "a") as f:
            write_result_line(res, f)
    return


def run_structure_analysis(
        item_list, 
        outdir,
        align_item=None,
        fname="structure_analysis.tsv",
        verbosity=1,
        use_pbar=False,
):
    printv = get_printv(verbosity)
    nitems = len(item_list)
    os.makedirs(FETCHED_PDBID_CACHE, exist_ok=True)

    pymol.finish_launching(["pymol", "-cq"])  # quiet + no GUI

    if align_item:
        load_pymol_item("struct0", align_item, fetch_path=FETCHED_PDBID_CACHE)

    results = []
    for i in tqdm.trange(nitems, desc="Item", disable=not use_pbar):
        item = item_list[i]
        printv(f"Loading {item.value}", pbar=use_pbar)
        load_pymol_item(f"struct", item, fetch_path=FETCHED_PDBID_CACHE)
        protein_struct = load_pdb_structure(item.value, id=item.id)
        results.append(analyze_protein(protein_struct))
        cmd.delete("struct")

    write_results(results, outfpath=f"{outdir}/{fname}")

    # Finish PyMOL session
    cmd.quit()
    return


def main(args):
    pdb_files = args.pdb_files if args.pdb_files else []
    names = args.names if args.names else []
    outdir = args.outdir
    outfname = args.outfname
    verbosity = args.verbosity
    use_pbar = args.pbar
    
    if len(pdb_files) != len(names):
        raise RuntimeError("Length of names should match length of files")
    
    if verbosity:
        print("Loading files:", pdb_files)
        print("Saving output to", outdir)

    loaded_structures = []
    for i, fpath in enumerate(pdb_files):
        loaded_structures.append(Item(fpath, "file", id=names[i]))

    os.makedirs(outdir, exist_ok=True)
    run_structure_analysis(
        loaded_structures, 
        outdir=outdir,
        fname=outfname,
        verbosity=verbosity,
        use_pbar=use_pbar,
    )

    return


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
