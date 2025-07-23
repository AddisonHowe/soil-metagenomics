#!/usr/bin/env bash
#=============================================================================
#
# FILE: run_pairwise_alignment.sh
#
# USAGE: run_pairwise_alignment.sh pdbdir outdir
#
# DESCRIPTION: Run pymol_pairwise_align.py script on all pdb files in a directory.
#
#   Args:
#       pdbdir: directory containing .pdb files to analyze.
#       outdir: directory to store output.
#
#
# EXAMPLE: sh run_pairwise_alignment.sh \
#               out/structure/<KO> \
#               out/structure_analysis/<KO>
#=============================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pdbdir> <outdir>"
    exit 1
fi

pdbdir=$1       # "out/structure/<KO>"
outdir=$2       # "out/structure_analysis/<KO>"

filenames=$(ls -l $pdbdir)
suffix=".pdb"

files=()
names=()
for fpath in ${pdbdir}/*.pdb; do
    [[ -f "$fpath" ]] || continue  # skip if no matching files
    filename=$(basename "${fpath%${suffix}}")
    files+=("$fpath")
    names+=("$filename")
done
    
python alphafold2/pymol_pairwise_align.py \
    -f "${files[@]}" -n "${names[@]}" \
    -o ${outdir} --pbar
