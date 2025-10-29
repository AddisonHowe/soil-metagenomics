#!/bin/bash

# Input files with your actual names
REFERENCE="K00370_rep.faa"
CONTIGS="JW70609_contigs.fasta"
SAMPLES="complete.fasta"  # This contains PROTEIN sequences

echo "=== DIAMOND GENE VARIANT ANALYSIS ==="
echo "Reference: $REFERENCE"
echo "Contigs: $CONTIGS" 
echo "Samples: $SAMPLES"
echo ""

echo "Step 1: Creating reference database..."
diamond makedb --in $REFERENCE -d reference_db

echo "Step 2: Searching contigs against reference..."
diamond blastx --query $CONTIGS \
               --db reference_db.dmnd \
               --out contigs_vs_reference.tsv \
               --outfmt 6 \
               --evalue 1e-20 \
               --max-target-seqs 1 \
               --sensitive \
               --threads 4

echo "Step 3: Filtering for BEST hit only..."
if [ -s contigs_vs_reference.tsv ]; then
    sort -k12,12gr contigs_vs_reference.tsv | head -n 1 > best_contig_hit.tsv
    echo "Best hit found in contigs:"
    cat best_contig_hit.tsv
else
    echo "ERROR: No hits found in contigs vs reference!"
    exit 1
fi

echo "Step 4: Extracting ONLY the best contig..."
best_contig=$(awk '{print $1}' best_contig_hit.tsv)
echo "Extracting contig: $best_contig"
echo "$best_contig" > best_contig.txt
seqtk subseq $CONTIGS best_contig.txt > best_contig_nucleotide.fasta

echo "Step 5: Translating best contig to protein..."
transeq best_contig_nucleotide.fasta best_contig_protein.fasta -frame 6

echo "Step 6: Creating protein database from best contig..."
diamond makedb --in best_contig_protein.fasta -d best_contig_db

echo "Step 7: Searching samples against BEST contig..."
# Use standard output format without scovhsp
diamond blastp --query $SAMPLES \
               --db best_contig_db.dmnd \
               --out samples_vs_best_contig.tsv \
               --outfmt 6 \
               --evalue 1e-10 \
               --max-target-seqs 1 \
               --sensitive \
               --threads 4

echo "Step 8: Checking if search produced results..."
if [ ! -s samples_vs_best_contig.tsv ]; then
    echo "WARNING: No results in samples_vs_best_contig.tsv"
    echo "Trying with more permissive e-value..."
    diamond blastp --query $SAMPLES \
                   --db best_contig_db.dmnd \
                   --out samples_vs_best_contig.tsv \
                   --outfmt 6 \
                   --evalue 1e-3 \
                   --max-target-seqs 1 \
                   --threads 4
fi

echo "Step 9: Basic filtering and sorting..."
if [ -s samples_vs_best_contig.tsv ]; then
    # Simple filter: identity > 50% and bitscore > 50
    awk '$3 > 50 && $12 > 50' samples_vs_best_contig.tsv > filtered_samples_vs_best_contig.tsv
    
    if [ -s filtered_samples_vs_best_contig.tsv ]; then
        sort -k12,12gr filtered_samples_vs_best_contig.tsv > sorted_samples_vs_best_contig.tsv
        echo "Top 10 closest matches in complete.fasta:"
        head -10 sorted_samples_vs_best_contig.tsv | awk '{printf "%-40s %-40s %5.1f%% identity %8.1f bitscore\n", $1, $2, $3, $12}'
        echo "Total filtered matches: $(wc -l < filtered_samples_vs_best_contig.tsv)"
    else
        echo "No matches passed basic filters. Showing raw results:"
        sort -k12,12gr samples_vs_best_contig.tsv > sorted_samples_vs_best_contig.tsv
        head -10 sorted_samples_vs_best_contig.tsv | awk '{printf "%-40s %-40s %5.1f%% identity %8.1f bitscore\n", $1, $2, $3, $12}'
    fi
else
    echo "ERROR: No results found even with permissive parameters!"
    echo "Check that your complete.fasta file contains valid protein sequences."
fi

echo "Analysis complete!"
