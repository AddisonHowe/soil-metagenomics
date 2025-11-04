#!/usr/bin/env bash

# Length 1100 threshold
python scripts/preprocess_orfs.py \
    -ko K00370 \
    -i out/aaseqs/orf_sequences_K00370_90.tsv \
    -o out/aaseqs/K00370 \
    -r data/reference_seqs/K00370/references.fasta \
    -t 1100 \
    --alphafold_batchsize 6
