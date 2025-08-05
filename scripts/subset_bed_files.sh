#!/bin/bash
#=============================================================================
#
# FILE: subset_bed_files.sh
#
# USAGE: subset_bed_files.sh ko rawdatapath outdir
#
# DESCRIPTION: Subset bed files in the raw data directory by KO value, and 
#   store results in the output directory.
#
# EXAMPLE: sh subset_bed_files.sh <ko> <rawdatapath> <outdir>
#=============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 ko rawdatapath outdir"
    exit 1
fi

koval="$1"          # <KO>
rawdatapath="$2"    # /scratch/aie7773
outdir="$3"         # data/KO_subsets/subset_<KO>

# bed files
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
