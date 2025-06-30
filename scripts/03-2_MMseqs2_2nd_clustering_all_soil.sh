#!/bin/bash
set -e  # Exit script if any command fails

# Usage: ./02_MMseqs2_each_soil_separately.sh <input_faa> <output_dir> <tmp_dir> <n_threads>
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <n_threads>"
    exit 1
fi

# Get the arguments
n_threads=$1

# Set the number of threads for MMseqs2
export OMP_NUM_THREADS=$n_threads

echo "**************************************"
echo "************** MMseqs2 using $n_threads threads ********************"
echo "**************************************"


# Define input and output directories
input_faa="../data/raw_data/all.coassembly_proteins_1st_ClusterDB_repseq.fasta"
echo "Processing file: $input_faa"
output_dir="../data/raw_data"
tmp_dir="${output_dir}/tmp_mmseqs"

# Create output and temp directories if they don't exist
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"

# Extract the base filename (without extension)
base_name=$(basename "$input_faa" .fasta)

# Define MMseqs2 database names
db_name="$tmp_dir/${base_name}_DB"
cluster_db="$tmp_dir/${base_name}_2ndClusterDB"


# 1) Create an MMseqs2 database from the protein sequences
echo "********************* Creating MMseqs2 database... ***************************"
mmseqs createdb "$input_faa" "$db_name"

# 2) Perform clustering
echo "********************* Clustering sequences... *****************************"
mmseqs cluster "$db_name" "$cluster_db" "$tmp_dir" \
       --min-seq-id 0.8 -c 0.7 --split-memory-limit 80G --remove-tmp-files

# 3) cluster result1
#echo "********************* Get cluster result in a fasta-like format... **************************"
#mmseqs createseqfiledb "$db_name" "$cluster_db" "${db_name}_rep"
#mmseqs result2flat "$db_name" "$db_name" "${db_name}_rep" "$rep_seqs"

# 4) cluster result2: Extract representative sequence
echo "********************* Getting representative sequences to FASTA... **************************"
mmseqs createsubdb "$cluster_db" "$db_name" "${cluster_db}_subdb"
mmseqs convert2fasta "${cluster_db}_subdb" "${cluster_db}_repseq.fasta"
mv "${cluster_db}_repseq.fasta" "$output_dir/"

# 5) Make tsv file for cluster membership
echo "********************* Make tsv file for cluster membership... **************************"
mmseqs createtsv "$db_name" "$db_name" "$cluster_db" "${cluster_db}.tsv"

# 6) calculate the frequency table
awk '{print $1}' "${cluster_db}.tsv" | sort | uniq -c | awk '{print $1}' | sort | uniq -c | sort -nr | awk '{print $2 "\t" $1}' > "${cluster_db}_frequency.tsv"

mv "${tmp_dir}/*.tsv" "$output_dir/"

echo "Finished processing"
echo "Output saved"
echo "-----------------------------------------------------"


# Cleanup temporary files (optional)
#rm -rf "$tmp_dir"/*_DB*  # Remove database files
#rm -rf "$tmp_dir"/*_rep* # Remove representative sequence files

echo "***********************"
echo "Done."
echo "***********************"
