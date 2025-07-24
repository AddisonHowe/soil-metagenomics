#!/usr/bin/env bash
#=============================================================================
#
# FILE: run_generate_images.sh
#
# USAGE: run_generate_images.sh pdbdir outdir alignto
#
# DESCRIPTION: Run pymol_generate_images.py script on all pdb files in a directory.
#
#   Args:
#       pdbdir: directory containing .pdb files to analyze.
#       outdir: directory to store output.
#       alignto: PDB ID of protein to align to.
#
#
# EXAMPLE: sh run_generate_images.sh \
#               out/structure/<KO> \
#               out/structure_analysis/<KO>/images \
#               <PDBID>
#=============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <pdbdir> <outdir> <alignto>"
    exit 1
fi

pdbdir=$1       # "out/structure/<KO>"
outdir=$2       # "out/structure_analysis/<KO>/images"
alignto=$3      # <PDBID>

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

python alphafold2/pymol_generate_images.py \
    -f "${files[@]}" --saveas_files "${names[@]}" \
    -o ${outdir} \
    -ai ${alignto} --pbar
