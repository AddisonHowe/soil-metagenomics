#!/bin/bash
#SBATCH --account=p32818
#SBATCH --partition=normal
#SBATCH --job-name=aa                  # Job name
#SBATCH --ntasks=1                     # Number of tasks 
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=120G                     # Total memory
#SBATCH --time=47:00:00                # Time limit (hh:mm:ss)


DIR='/projects/p32818/metagenomic_data/scripts'

cd $DIR

# Run script
python aa_sequence.py 
