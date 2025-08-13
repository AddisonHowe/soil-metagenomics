#!/usr/bin/env bash
#=============================================================================
#
# FILE: search_combined_genome.sh
#
# USAGE: search_combined_genome.sh searchfile genome outdir
#
# DESCRIPTION: Search the combined genome file `genome` for search terms listed 
#   in the input file `searchfile`. Each line of `searchfile` consists of two 
#   tab-separated values: a search version and a term to search.
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

while IFS= read -r line; do
    [[ "$line" =~ ^# ]] && continue  # skip comments
    IFS=$'\t' read -r searchtype search <<< "$line"

    if [ $searchtype == "basic" ]; then
        grep_start="grep -i \""
        grep_end="\""
    elif [ $searchtype == "exact" ]; then
        grep_start="grep \""
        grep_end="\""
    elif [ $searchtype == "gene" ]; then
        grep_start="grep -i -E \"gene="
        grep_end="(;|\$)\""
    elif [ $searchtype == "name" ]; then
        grep_start="grep -i -E \"Name="
        grep_end="(;|\$)\""
    else
        echo "searchtype must be one of: basic, exact, gene, name"
        exit 1
    fi

    grep_cmd=${grep_start}${search}${grep_end}
    cat ${genomefpath} | eval "$grep_cmd" > ${outdir}/${searchtype}_${search}
    echo ${searchtype}_$search "("$(cat "${outdir}/${searchtype}_${search}" | wc -l) results")"
done < ${searchfile}
