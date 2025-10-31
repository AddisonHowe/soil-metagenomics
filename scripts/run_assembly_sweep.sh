#!/bin/bash
#SBATCH --account=p32818
#SBATCH --partition=normal
#SBATCH --job-name=assembly_pipeline
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=assembly_pipeline_%j.out

# Load required modules
module load spades
module load prodigal

# Define coverage and error rate combinations
coverages=(1 2 3 5 8 13 22 36 60 100)
error_rates=(0.000 0.001 0.005)

# Create output directory for ORFs only
mkdir -p ../out/reconstruction/orfs

for cov in "${coverages[@]}"; do
    for err in "${error_rates[@]}"; do
        echo "Processing coverage ${cov}, error rate ${err}"
        
        # Define file paths
        R1="../out/reconstruction/raw_reads/spikein_R1_${cov}_${err}.fasta"
        R2="../out/reconstruction/raw_reads/spikein_R2_${cov}_${err}.fasta"
        ASSEMBLY_DIR="../out/reconstruction/temp_assembly_${cov}_${err}"
        ORFS_PREFIX="../out/reconstruction/orfs/orfs_${cov}_${err}"
        
        # Step 1: Assembly with metaSPAdes
        echo "Running metaSPAdes..."
        metaspades.py \
            -1 $R1 \
            -2 $R2 \
            -o $ASSEMBLY_DIR \
            -t 8 \
            -m 64 \
            --only-assembler
        
        # Step 2: ORF prediction with Prodigal directly from assembly output
        echo "Running Prodigal..."
        prodigal -i $ASSEMBLY_DIR/contigs.fasta \
                 -a ${ORFS_PREFIX}.faa \
                 -p meta \
                 -q
        
        # Clean up temporary assembly directory (removes contigs too)
        rm -rf $ASSEMBLY_DIR
        
        echo "Completed coverage ${cov}, error rate ${err}"
        echo "ORFs: ${ORFS_PREFIX}.faa"
        echo "----------------------------------------"
    done
done

echo "All assemblies and ORF predictions completed!"
