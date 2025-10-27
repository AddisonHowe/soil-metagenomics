#!/bin/bash

# Input files
REFERENCE="K00370_rep.faa"
CONTIGS="JW70607_contigs.fasta" 
SAMPLES="complete.fasta"  # This contains PROTEIN sequences

echo "Step 1: Creating reference database..."
diamond makedb --in $REFERENCE -d reference_db

echo "Step 2: Searching contigs against reference..."
diamond blastx --query $CONTIGS \
               --db reference_db.dmnd \
               --out contigs_vs_reference.tsv \
               --outfmt 6 \
               --evalue 1e-10 \
               --max-target-seqs 1 \
               --threads 4

echo "Step 3: Extracting positive contigs..."
awk '{print $1}' contigs_vs_reference.tsv | sort -u > positive_contigs.txt
seqtk subseq $CONTIGS positive_contigs.txt > extracted_genes_nucleotide.fasta

echo "Step 4: TRANSLATE nucleotide contigs to protein..."
# Method 1: Using transeq (EMBOSS) - recommended
transeq extracted_genes_nucleotide.fasta extracted_genes_protein.fasta -frame 6

# Method 2: If transeq not available, use ORFfinder or similar
# Or use this Python one-liner as alternative:
# python -c "from Bio.Seq import Seq; from Bio import SeqIO; [print(f'>{rec.id}\\n{Seq(rec.seq).translate()}') for rec in SeqIO.parse('extracted_genes_nucleotide.fasta', 'fasta')]" > extracted_genes_protein.fasta

echo "Step 5: Create protein database from translated sequences..."
diamond makedb --in extracted_genes_protein.fasta -d extracted_genes_db

echo "Step 6: Search protein samples against protein database..."
diamond blastp --query $SAMPLES \
               --db extracted_genes_db.dmnd \
               --out samples_vs_extracted.tsv \
               --outfmt 6 \
               --evalue 1e-5 \
               --max-target-seqs 10 \
               --threads 4

echo "Analysis complete!"
echo "Found $(wc -l < samples_vs_extracted.tsv 2>/dev/null || echo 0) hits in final search"
