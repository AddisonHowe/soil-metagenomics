#!/bin/bash
#SBATCH --account=p32818
#SBATCH --partition=long
#SBATCH --job-name=mmseq               # Job name
#SBATCH --ntasks=1                     # Number of tasks 
#SBATCH --cpus-per-task=32             # Number of CPU cores per task
#SBATCH --mem=300G                     # Total memory
#SBATCH --time=120:00:00                # Time limit (hh:mm:ss)


#Load MMseq
module load mmseqs2/15-6f452-openmpi-intel-2021.4.0


DIR='/projects/p32818/metagenomic_data/scripts'

cd $DIR


# Run script
bash 03-2_MMseqs2_2nd_clustering_all_soil.sh 32
