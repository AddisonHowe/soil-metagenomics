#!/bin/bash
#SBATCH --account=p32818
#SBATCH --partition=short
#SBATCH --job-name=prodigal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=prodigal_%j.out

# Load Prodigal module (adjust for your cluster)
module load prodigal

# Input and output files
CONTIGS="../out/reconstruction/metaspades_assembly/contigs.fasta"
OUTPUT_PREFIX="../out/reconstruction/prodigal_orfs"

echo "Running Prodigal on contigs..."
echo "Input: $CONTIGS"

# Run Prodigal in meta mode for metagenomic contigs
prodigal -i $CONTIGS \
         -o ${OUTPUT_PREFIX}.gff \
         -a ${OUTPUT_PREFIX}.faa \
         -d ${OUTPUT_PREFIX}.fna \
         -p meta \
         -f gff

echo "Prodigal completed!"
echo "Output files:"
echo "  Protein sequences: ${OUTPUT_PREFIX}.faa"
echo "  DNA sequences: ${OUTPUT_PREFIX}.fna" 
echo "  Gene annotations: ${OUTPUT_PREFIX}.gff"

# Count the number of ORFs found
NUM_ORFS=$(grep -c ">" ${OUTPUT_PREFIX}.faa)
echo "Number of ORFs predicted: $NUM_ORFS"
