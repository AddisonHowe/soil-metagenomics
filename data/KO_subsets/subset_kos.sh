#!/bin/bash

annotation_dir=data/raw_data/coassembly
bedfiledir=/scratch/aie7773

while IFS= read -r koval; do
    echo "${koval}"
    outdir=data/KO_subsets/subset_${koval}
    echo Subsetting annotation data...
    sh scripts/subset_annotation_data.sh $koval $annotation_dir $outdir
    echo Subsetting bed files...
    sh scripts/subset_bed_files.sh $koval $bedfiledir $outdir
done < data/KO_subsets/ko_list.txt
