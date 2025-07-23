"""Pymol image generation script

Usage: 
    python pymol_generate_images.py [-f <pdb_file1.pdb> <pdb_file2.pdb> ...] 
                                    [-i <pdb_id1> <pdb_id2> ...] 
                                    [-af <align_to_file>]
                                    [-ai <align_to_id>]
                                    [-sf <savefile1> <savefile2> ...]
                                    [-si <saveid1> <saveid2> ...]
                                    [-o <outdir>]
Inputs:
    [-f --pdb_files] : pdb filepaths to view.
    [-i --pdb_ids] : PDB IDs to fetch from the online database.
    [-af --align_to_file] : pdb filepath to align to
    [-ai --align_to_id] : PDB ID to align to
    [-sf --saveas_files] : saveas names for files
    [-si --saveas_ids] : saveas names for ids
    [-o --outdir] : output directory

Outputs:
    Creates an output directory <outdir> storing results of each structure.

"""

import sys, os
import argparse
import tqdm
import pymol
from pymol import cmd

from pymol_helpers import Item, load_pymol_item, get_printv


FETCHED_PDBID_CACHE = "./out/pymol/fetched_pdbids"
ALIGN_TO_COLOR = 'cyan'
COLOR_CMD = "spectrum b, red blue,"


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--pdb_files", type=str, nargs="*", 
                        help="pdb files to load")
    parser.add_argument("-i", "--pdb_ids", type=str, nargs="*", 
                        help="pdb IDs to fetch from PDB database")
    parser.add_argument("-af", "--align_to_file", type=str, default=None, 
                        help="pdb filepath to align to")
    parser.add_argument("-ai", "--align_to_id", type=str, default=None, 
                        help="PDB ID to align to")
    parser.add_argument("-sf", "--saveas_files", type=str, nargs="*", 
                        default=None, help="saveas names")
    parser.add_argument("-si", "--saveas_ids", type=str, nargs="*", 
                        default=None, help="saveas names")
    parser.add_argument("-o", "--outdir", type=str, 
                        help="output directory")
    parser.add_argument("-v", "--verbosity", type=int, default=1)
    parser.add_argument("--pbar", action="store_true")
    return parser.parse_args(args)


def run_pymol_generate_images(
        item_list, 
        outdir,
        align_item=None,
        saveas=None,
        verbosity=1,
        use_pbar=False,
):
    printv = get_printv(verbosity)
    nitems = len(item_list)
    os.makedirs(FETCHED_PDBID_CACHE, exist_ok=True)

    pymol.finish_launching(["pymol", "-cq"])  # quiet + no GUI

    if align_item:
        load_pymol_item("struct0", align_item, fetch_path=FETCHED_PDBID_CACHE)

    for i in tqdm.trange(nitems, desc="Item", disable=not use_pbar):
        item = item_list[i]
        
        printv(f"Loading {item.value}", pbar=use_pbar)
        load_pymol_item(f"struct", item, fetch_path=FETCHED_PDBID_CACHE)
        
        cmd.spectrum("b", "red blue", "struct")
        
        if align_item:
            cmd.align("struct", "struct0")
            cmd.color(ALIGN_TO_COLOR, "struct0")

        cmd.zoom("all")

        if saveas:
            fname = saveas[i]
        else:
            fname = f"structure_{i}"
        
        printv(f"Saving as {outdir}/{fname}.png", pbar=use_pbar)

        cmd.png(
            f"{outdir}/{fname}.png", 
            # width=1920, 
            # height=1080, 
            dpi=300, 
            ray=1
        )
        cmd.delete("struct")
    
    # Finish PyMOL session
    cmd.quit()
    return


def main(args):
    pdb_files = args.pdb_files if args.pdb_files else []
    pdb_ids = args.pdb_ids if args.pdb_ids else []
    saveas_files = args.saveas_files if args.saveas_files else []
    saveas_ids = args.saveas_ids if args.saveas_ids else []
    outdir = args.outdir
    verbosity = args.verbosity
    use_pbar = args.pbar

    if len(pdb_files) != len(saveas_files):
        raise RuntimeError("Length of saveas_files should match length of files")
    if len(pdb_ids) != len(saveas_ids):
        raise RuntimeError("Length of saveas_ids should match length of ids")
    
    if args.align_to_file:
        align_to = args.align_to_file
        align_item = Item(align_to, "file")
    elif args.align_to_id:
        align_to = args.align_to_id
        align_item = Item(align_to, "id")
    else:
        align_to = None
        align_item = None

    
    print(f"Loading {len(pdb_files)} files")
    print(f"Fetching {len(pdb_ids)} ids")
    if align_item:
        print(f"Aligning to {align_item.kind} {align_item.value}")
    print("Saving output to", outdir)

    loaded_structures = []
    saveas_names = []
    
    for i, fpath in enumerate(pdb_files):
        struct = Item(fpath, "file")
        loaded_structures.append(struct)
        if saveas_files:
            saveas_names.append(saveas_files[i])
    
    for i, id in enumerate(pdb_ids):
        struct = Item(id, "id")
        loaded_structures.append(struct)
        if saveas_ids:
            saveas_names.append(saveas_ids[i])

    os.makedirs(outdir, exist_ok=True)
    
    run_pymol_generate_images(
        loaded_structures, 
        outdir=outdir,
        align_item=align_item,
        saveas=saveas_names,
        verbosity=verbosity,
        use_pbar=use_pbar,
    )

    return


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
