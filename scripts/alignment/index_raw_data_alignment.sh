#!/usr/bin/env bash
#=============================================================================
#
# FILE: index_raw_data_alignment.sh
#
# USAGE: index_raw_data_alignment.sh prefix datdir
#
# DESCRIPTION: Index a BAM file and modify the header to replace commas in the
#   @SQ sequence names with underscores.
#
#   Args:
#       prefix: prefix of the fastq.gz file to process.
#       datdir: directory containing fastq.gz files.
#   Outputs:
#       Generates the following files in <datdir>/<prefix>.fastq.gz/
#           header.sam: updated SAM header
#           53005.2.525483.[TAG].fixed.sorted.bam: updated BAM file with header
#           53005.2.525483.[TAG].fixed.sorted.bam.bai: BAM index
#
# EXAMPLE: sh index_raw_data_alignment.sh 53016.2.536117.[TAG] raw_data_alignment
#=============================================================================

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 prefix datdir"
    exit 1
fi

prefix=$1
datdir=$2

subdir=${prefix}.fastq.gz
bamfname=${prefix}.sorted.bam
bamfpath=${datdir}/${subdir}/${bamfname}

echo "Processing data in ${datdir}/${subdir}"

# Build header.sam file
header_fpath=${datdir}/${subdir}/header.sam
fixed_bam_fpath=${datdir}/${subdir}/${bamfname/.sorted/.fixed.sorted}
samtools view -H ${bamfpath} > ${header_fpath}
echo "  Generated header.sam file"

# Fix bam file by modifying the header
sed -i.bak 's/,/_/g' ${header_fpath}
samtools reheader ${header_fpath} ${bamfpath} > ${fixed_bam_fpath}
samtools index ${fixed_bam_fpath}
echo "  Generated fixed.sorted.bam file"

echo "  Done"
