#!/bin/bash
#SBATCH -A p32818
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mem=2gb
#SBATCH --array=0-38
#SBATCH --job-name="subset_data"
#SBATCH --output="logs/subset_data.%A_%a.o"
#SBATCH --error="logs/subset_data.%A_%a.e"
#SBATCH --mail-user=addi.howe@gmail.com
#SBATCH --mail-type=END,FAIL

echo "Slurm job ID: ${SLURM_JOB_ID}"
echo "Slurm task ID: ${SLURM_ARRAY_TASK_ID}"

cd /projects/p32818/metagenomic_data

# Read KO values from the list of all KO values to store
kolist_fpath=data/KO_subsets/ko_list.txt
mapfile -t kovals < "$kolist_fpath"

for i in "${!kovals[@]}"; do
    if (( i % SLURM_ARRAY_TASK_COUNT == SLURM_ARRAY_TASK_ID )); then
        koval="${kovals[i]}"
        echo "Handling task ${i}: $koval"
        sh data/KO_subsets/subset_ko.sh $koval
    fi
done
