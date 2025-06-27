#!/bin/bash

filelist=(
    T0_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil3_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil5_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil6_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil9_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil11_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil12_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil14_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil15_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil16_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
    # Soil17_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed
)

suffix=_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed

for f in ${filelist[@]}; do 
    echo $f
    cat data/shared_data/${f} | awk '$2 == "K00370" { 
            print $1, $2, $3
        }' OFS='\t' > data/subset_K00370/${f/${suffix}/_all_samples_K00370}.bed
done
