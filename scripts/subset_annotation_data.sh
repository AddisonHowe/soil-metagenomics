#!/bin/bash

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

suffix=.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv

for f in ${filelist[@]}; do 
    echo $f
    cat data/shared_data/${f} | awk -F'\t' '$4 == "K00370" { 
            print $1, $2, $3, $4, $5
        }' OFS='\t' > data/subset_K00370/${f/${suffix}/.coassembly_annotations_K00370}.tsv
done
