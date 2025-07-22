#!/bin/bash

rawdatapath=../data/raw_data
koval=K00371
outdir=../data/subset_${koval}

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
    cat ${rawdatapath}/${f} | awk -F'\t' -v koval="$koval" \
        '$4 == koval { 
            print $1, $2, $3, $4, $5
        }' \
        OFS='\t' > ${outdir}/${f/${oldsuffix}/${newsuffix}}
done
