#!/usr/bin/env bash
#=============================================================================
#
# FILE: search_combined_genome.sh
#
# USAGE: search_combined_genome.sh searchfile genome outdir
#
# DESCRIPTION: Search the combined genome file `genome` for search terms listed 
#   in the input file `searchfile`. Each line of `searchfile` is a term to 
#   search.
#
# EXAMPLE: sh search_combined_genome.sh searchfile genome output
#=============================================================================


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 searchfile genome outdir"
    exit 1
fi


searchfile=$1
genomefpath=$2  #data/raw_data/sequencing/combined_genome.gff
outdir=$3

mkdir -p $outdir

while IFS= read -r search; do 
    cat ${genomefpath} | grep -i "$search" > ${outdir}/${search}
    echo $search "("$(cat "${outdir}/${search}" | wc -l) results")"
done < ${searchfile}
