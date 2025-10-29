#!/bin/bash
#SBATCH --account=p32818
#SBATCH --job-name=diamond_search
#SBATCH --output=diamond_%j.out
#SBATCH --error=diamond_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=short

# Load required modules
module load diamond seqtk

# Run your script
./run_search.sh
