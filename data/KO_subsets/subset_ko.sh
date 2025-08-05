#!/bin/bash
#=============================================================================
#
# FILE: subset_ko.sh
#
# USAGE: subset_ko.sh koval
#
# DESCRIPTION: Wrapper to subset a KO value.
#
# EXAMPLE: sh subset_ko.sh <koval>
#=============================================================================

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 koval"
    exit 1
fi

koval=$1
annotation_dir=data/raw_data/coassembly
bedfiledir=/scratch/aie7773

echo KO="${koval}"
outdir=data/KO_subsets/subset_${koval}
mkdir -p $outdir

echo Subsetting annotation data...
sh scripts/subset_annotation_data.sh $koval $annotation_dir $outdir

echo Subsetting bed files...
sh scripts/subset_bed_files.sh $koval $bedfiledir $outdir

echo Done!
