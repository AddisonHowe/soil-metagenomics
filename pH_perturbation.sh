#!/bin/bash
#SBATCH --account=p32818
#SBATCH --partition=short
#SBATCH --job-name=pH_perturbation     # Job name
#SBATCH --ntasks=1                     # Number of tasks 
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=32G                       # Total memory
#SBATCH --time=02:00:00                # Time limit (hh:mm:ss)


DIR='/projects/p32818/metagenomic_data'

cd $DIR

# Run script
python pH_perturbation.py 'Soil3'
