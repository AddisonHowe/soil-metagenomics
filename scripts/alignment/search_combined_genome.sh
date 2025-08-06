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


if [ "$#" -eq 3 ]; then
    version="basic"
elif [ "$#" -eq 4 ]; then
    version=$4
else
    echo "Usage: $0 searchfile genome outdir [version]"
    exit 1
fi

searchfile=$1
genomefpath=$2  #data/raw_data/sequencing/combined_genome.gff
outdir=$3

# Process version
if [ $version == "basic" ]; then
    grep_start="grep -i \""
    grep_end="\""
elif [ $version == "exact" ]; then
    grep_start="grep \""
    grep_end="\""
elif [ $version == "gene" ]; then
    grep_start="grep -i -E \"gene="
    grep_end="(;|\$)\""
else
    echo "version must be one of: basic, exact, gene"
    exit 1
fi


mkdir -p $outdir

while IFS= read -r search; do
    [[ "$search" =~ ^# ]] && continue  # skip comments
    grep_cmd=${grep_start}${search}${grep_end}
    cat ${genomefpath} | eval "$grep_cmd" > ${outdir}/${search}
    echo $search "("$(cat "${outdir}/${search}" | wc -l) results")"
done < ${searchfile}
