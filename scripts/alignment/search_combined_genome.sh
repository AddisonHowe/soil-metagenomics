#!/usr/bin/env bash
#=============================================================================
#
# FILE: search_combined_genome.sh
#
# USAGE: search_combined_genome.sh searchfile outdir
#
# DESCRIPTION: Search the combined genome file for search terms listed in the 
#   input file `searchfile`. Each line of `searchfile` is a term to search.
#
# EXAMPLE: sh search_combined_genome.sh /path/to/searchfile /path/to/output
#=============================================================================


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 searchfile outdir"
    exit 1
fi


searchfile=$1
outdir=$2

genomefpath=data/raw_data/sequencing/combined_genome.gff

mkdir -p $outdir

for search in $(cat ${searchfile}); do 
    cat ${genomefpath} | grep -i $search > ${outdir}/${search}
    echo $search "("$(cat ${outdir}/${search} | wc -l) results")"
done
