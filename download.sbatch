#!/bin/sh
#SBATCH -A p31512
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=1gb
#SBATCH --job-name="interpro"
#SBATCH --output="logs/interpro.%A.o"
#SBATCH --error="logs/interpro.%A.e"
#SBATCH --mail-user=addi.howe@gmail.com
#SBATCH --mail-type=END,FAIL

echo Slurm job ID: $SLURM_JOB_ID

cd /projects/p31512/aeh581/soil-metagenomics

module purge
eval "$($CONDA_EXE shell.bash hook)"
conda activate env


cd data/interpro

python download_script.py
