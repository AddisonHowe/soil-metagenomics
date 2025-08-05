#!/bin/bash
#=============================================================================
#
# FILE: subset_annotation_data.sh
#
# USAGE: subset_annotation_data.sh ko annotation_dir outdir
#
# DESCRIPTION: Subset annotation data by KO value, and store results in the 
#   output directory.
#
# EXAMPLE: sh subset_annotation_data.sh <ko> <annotation_dir> <outdir>
#=============================================================================

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 ko annotation_dir outdir"
    exit 1
fi

koval="$1"          # <KO>
annotation_dir="$2" # data/raw_data/coassembly
outdir="$3"         # data/KO_subsets/subset_<KO>

# Annotation files
filelist=(
    T0.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil3.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil5.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil6.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil9.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil11.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil12.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil14.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil15.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil16.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
    Soil17.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
)

oldsuffix=.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv
newsuffix=.coassembly_annotations_${koval}.tsv

for f in ${filelist[@]}; do 
    echo $f
    cat ${annotation_dir}/${f} | awk -F'\t' -v koval="$koval" \
        '$4 == koval { 
            print $1, $2, $3, $4, $5
        }' \
        OFS='\t' > ${outdir}/${f/${oldsuffix}/${newsuffix}}
done
