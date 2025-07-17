#!/usr/bin/env bash
#=============================================================================
#
# FILE: run_alphafold_cpu.sh
#
# USAGE: run_alphafold_cpu.sh fasta_paths outdir
#
# DESCRIPTION: Run the CPU component of AlphaFold on one or more fasta files. 
#  Should be run from the project directory, with `fasta_paths` and `outdir` 
#  arguments absolute paths. Comma separate multiple fasta files.
#
# EXAMPLE: sh run_alphafold_cpu.sh data/sequences.fasta out/alphafold
#=============================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 fpath outdir"
    exit 1
fi

fasta_paths=$1
outdir=$2

echo "Running alphafold2 (CPU) on fasta file(s) ${fasta_paths}"
echo "Saving to ${outdir}"

alphafold-monomer --fasta_paths=${fasta_paths} \
    --max_template_date=2022-01-01 \
    --model_preset=monomer \
    --db_preset=full_dbs \
    --only_msas=true \
    --use_gpu_relax=False \
    --output_dir=${outdir}
