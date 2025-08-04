"""

"""

import os, sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import csv
import gzip
import tqdm as tqdm
import time

FIX_NAMES = {
    "NC_010943.1_Stenotrophomonas_maltophilia_K279a__strain_K279a":
        "NC_010943.1_Stenotrophomonas_maltophilia_K279a,_strain_K279a",
}

INCLUDE_SCAFFOLDS = [
    "Lee_A8Q_1_Ecoli_contig_1_polypolish",
]

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--seqdir", type=str, required=True,
                        help="sequencing data directory "
                        "(e.g. data/raw_data/sequencing/)")
    parser.add_argument("--datdir", type=str, required=True,
                        help="raw data alignments directory "
                        "(e.g. data/raw_data/sequencing/raw_data_alignment).")
    parser.add_argument("--outdir", type=str, required=True,
                        help="output directory")
    parser.add_argument("--sam_header_fname", type=str, 
                        default="header.sam", help="Name of SAM header files.")
    parser.add_argument("--scaffolds", type=str, nargs="+", 
                        default=INCLUDE_SCAFFOLDS,
                        help="scaffold names to include.")
    parser.add_argument("--pbar", action="store_true")
    return parser.parse_args(args)


def main(args):    
    DISABLE_PBAR = not args.pbar
    SEQDIR = args.seqdir
    DATDIRBASE = args.datdir
    sam_header_fname = args.sam_header_fname
    INCLUDE_SCAFFOLDS = args.scaffolds
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)
    

    RAWNAME_TO_SAMPLE = {}
    SAMPLE_TO_RAWNAME = {}
    with open(f"{SEQDIR}/raw_data_name_mapping.tsv", "r") as f:
        csvreader = csv.reader(f, delimiter="\t")
        for row in csvreader:  # process each row
            RAWNAME_TO_SAMPLE[row[0]] = row[1]
            SAMPLE_TO_RAWNAME[row[1]] = row[0]


    sampleids = np.sort(list(SAMPLE_TO_RAWNAME.keys()))
    for sampleid in tqdm.tqdm(sampleids, total=len(sampleids), disable=DISABLE_PBAR):
        rawid = SAMPLE_TO_RAWNAME[sampleid]
        datdir = f"{DATDIRBASE}/{rawid}"
        if not os.path.isdir(datdir):
            continue
        basecov_fname = rawid.replace(".fastq.gz", ".sam.pileup.basecov.gz")
        basecov_fpath = f"{datdir}/{basecov_fname}"    
        # Read in lengths of contigs from SAM header
        scaffold_lengths = {}
        with open(f"{datdir}/{sam_header_fname}", "r") as f:
            csvreader = csv.reader(f, delimiter="\t")
            for row in csvreader:  # process each row
                if row[0] == "@SQ":
                    name = row[1].removeprefix("SN:")
                    name = FIX_NAMES.get(name, name)
                    n = int(row[2].removeprefix("LN:"))
                    scaffold_lengths[name] = n
        # Load data
        name_list, array_list = load_data(
            basecov_fpath, 
            scaffold_lengths,
            include_names=INCLUDE_SCAFFOLDS,
            disable_pbar=True,
            verbosity=1
        )
        named_arrays = {name: arr for name, arr in zip(name_list, array_list)}
        np.savez_compressed(f"{outdir}/coverage_arrays_{sampleid}.npz", **named_arrays)
    return


def load_data(
        basecov_fpath, 
        length_map,
        include_names="all",
        disable_pbar=False,
        verbosity=1,
):
    if include_names == "all":
        include_names = list(length_map.keys())
    prev_name = ""
    cnt = 0
    file_index = 0
    saved_files = 0
    total_files = len(length_map)
    name_list = []
    array_list = []
    time0 = time.time()
    skip_rows = 0
    with gzip.open(basecov_fpath, "rt") as f:
        csvreader = csv.reader(f, delimiter="\t",)
        next(csvreader)  # Skip header
        pbar = tqdm.tqdm(desc="Progress", disable=disable_pbar)
        while True:
            # Handle row skipping
            if skip_rows:
                # Skip specified number of rows
                for skipiter in range(skip_rows):
                    next(csvreader)
                    if skipiter % 1000 == 0:
                        pbar.update(1000)
                pbar.update(pbar.total - pbar.n)
                skip_rows = 0
            
            # Retrieve next row
            try:
                row = next(csvreader)
                name, position, nreads = row
                position = int(position)
                nreads = int(nreads)
            except StopIteration:
                break

            # Check if starting on new scaffold
            if prev_name != name:
                # Exit early if number of saved files will exceed num to save
                if saved_files == len(include_names):
                    pbar.update(pbar.total - pbar.n)
                    pbar.refresh()
                    return name_list, array_list
                # Update values to begin processing new scaffold
                prev_name = name
                file_index += 1
                cnt = 0
                nrows = length_map[name]
                # Update the progress bar
                if name in include_names:
                    desc = f"File {file_index}/{total_files}: {name}"
                else:
                    desc = f"File {file_index}/{total_files} (skip): {name}"
                pbar.set_description(desc)
                pbar.total = nrows
                pbar.refresh()
                pbar.reset()
                
                if name in include_names:
                    saved_files += 1
                else:
                    skip_rows = nrows - 1
                    pbar.update()
                    continue

                name_list.append(name)
                read_array = np.zeros(nrows, dtype=int)
                array_list.append(read_array)
                if disable_pbar and verbosity:
                    print(f"Reading {name}")

            # Update read array with read count at position 
            read_array[position] = nreads
            if cnt != position:
                raise RuntimeError("count and position mismatch")
            cnt += 1
            if cnt % 1000 == 0:
                pbar.update(1000)

        pbar.update(pbar.total - pbar.n)

    if verbosity:
        print(f"Finished loading data in {time.time() - time0:.4g} seconds")
    pbar.refresh()
    return name_list, array_list


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
