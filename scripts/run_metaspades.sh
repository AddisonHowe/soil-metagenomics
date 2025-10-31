#!/bin/bash
#SBATCH --account=p32818
#SBATCH --partition=normal
#SBATCH --job-name=metaspades
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=metaspades_%j.out

# Load required modules (adjust for your cluster)
module load spades

metaspades.py \
    -1 ../out/reconstruction/spikein_R1.fasta \
    -2 ../out/reconstruction/spikein_R2.fasta \
    -o ../out/reconstruction/metaspades_assembly \
    -t 8 \
    -m 64 \
    --only-assembler

echo "Assembly complete! Contigs: ../out/reconstruction/metaspades_assembly/contigs.fasta"
