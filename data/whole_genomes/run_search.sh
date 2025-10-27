#!/bin/bash

# Input files
REFERENCE="K00370_rep.faa"
CONTIGS="JW70607_contigs.fasta"
SAMPLES="complete.fasta"

# Step 1: Create reference database
echo "Creating reference database..."
diamond makedb --in $REFERENCE -d reference_db

# Step 2: First DIAMOND search
echo "Searching contigs against reference..."
diamond blastx --query $CONTIGS \
               --db reference_db.dmnd \
               --out contigs_vs_reference.tsv \
               --outfmt 6 \
               --evalue 1e-10 \
               --max-target-seqs 1 \
               --threads 4

# Step 3: Extract positive contigs (simple approach)
echo "Extracting positive contigs..."
awk '{print $1}' contigs_vs_reference.tsv | sort -u > positive_contigs.txt
seqtk subseq $CONTIGS positive_contigs.txt > extracted_genes.fasta

# Step 4: Create database from extracted genes
echo "Creating extracted genes database..."
diamond makedb --in extracted_genes.fasta -d extracted_genes_db

# Step 5: Second DIAMOND search
echo "Searching samples against extracted genes..."
diamond blastx --query $SAMPLES \
               --db extracted_genes_db.dmnd \
               --out samples_vs_extracted.tsv \
               --outfmt 6 \
               --evalue 1e-5 \
               --threads 4

echo "Analysis complete!"
echo "Results: contigs_vs_reference.tsv and samples_vs_extracted.tsv"
