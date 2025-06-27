#!/bin/bash

rawdatapath=data/raw_data
koval=K00370
outdir=data/subset_${koval}

filelist=(
    T0_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil3_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil5_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil6_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil9_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil11_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil12_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil14_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil15_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil16_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    Soil17_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
)

oldsuffix="_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed"
newsuffix="_all_samples_${koval}.bed"
annotsuffix=".coassembly_annotations_${koval}.tsv"

for f in ${filelist[@]}; do 
    echo $f
    awk -v koval="$koval" '
        FNR==NR { keys[$1]; next }
        $2 in keys
        ' ${outdir}/${f/${oldsuffix}/${annotsuffix}} ${rawdatapath}/$f \
        > ${outdir}/${f/${oldsuffix}/${newsuffix}}
done
