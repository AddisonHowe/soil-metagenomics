#!/usr/bin/env bash
#=============================================================================
#
# FILE: run_wrapper.sh
#
# USAGE: run_wrapper.sh
#
# DESCRIPTION: Run a series of analysis scripts on alphafold output.
#
# EXAMPLE: sh run_wrapper.sh
#=============================================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~ run_generate_images

sh scripts/alphafold2/run_generate_images.sh out/structure/K00370 \
    out/structure_analysis/K00370/images 1Q16

sh scripts/alphafold2/run_generate_images.sh out/structure/K02567 \
    out/structure_analysis/K02567/images 2NYA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~ run_pairwise_alignment

sh scripts/alphafold2/run_pairwise_alignment.sh \
    out/structure/K00370 out/structure_analysis/K00370

sh scripts/alphafold2/run_pairwise_alignment.sh \
    out/structure/K02567 out/structure_analysis/K02567

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~ run_structure_analysis

sh scripts/alphafold2/run_structure_analysis.sh out/structure/K00370 \
    out/structure_analysis/K00370 structure_metrics_K00370.tsv

sh scripts/alphafold2/run_structure_analysis.sh out/structure/K02567 \
    out/structure_analysis/K02567 structure_metrics_K02567.tsv
