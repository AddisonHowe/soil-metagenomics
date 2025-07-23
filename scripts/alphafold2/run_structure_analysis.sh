#!/usr/bin/env bash
#=============================================================================
#
# FILE: run_structure_analysis.sh
#
# USAGE: run_structure_analysis.sh pdbdir outdir outfname
#
# DESCRIPTION: Run structure_analysis.py script on all pdb files in a directory.
#
#   Args:
#       pdbdir: directory containing .pdb files to analyze.
#       outdir: directory to store output.
#       outfname: name of output file to create within the output directory.
#
# EXAMPLE: sh run_structure_analysis.sh \
#               out/structure/<KO> \
#               out/structure_analysis/<KO> \
#               structure_metrics_<KO>.tsv
#=============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <pdbdir> <outdir> <outfname>"
    exit 1
fi

pdbdir=$1       # "out/structure/<KO>"
outdir=$2       # "out/structure_analysis/<KO>"
outfname=$3     # "structure_metrics_<KO>.tsv"

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

python alphafold2/structure_analysis.py \
    -f "${files[@]}" --names "${names[@]}" \
    -o ${outdir} --outfname ${outfname} --pbar
