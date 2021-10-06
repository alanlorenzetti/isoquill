#!/bin/bash

# alorenzetti 202108

# description ####
# this set of scripts will
# perform de isoseq3 regular
# pipeline in order to process
# Qbra transcriptome obtained
# using PacBio technology

# the scripts will perform

# 01. the obtention of circular consensus sequence (CCS)
for i in {1..12} ; do 
    bash scripts/01_process_isoseq.sh yes $i 12
done

# 02. the merging of chunks obtained from step 01
bash scripts/02_merge_ccs.sh

# 03. trimming of primers and correction of tx orientation
bash scripts/03_trimming.sh

# 04. removing poly-A tails and concatamers
# followed by the clustering of transcripts
bash scripts/04_refine_and_cluster.sh

# 05. removing redundancy
bash scripts/05_cogent.sh

# 06. using busco to check completeness
bash scripts/06_busco.sh

# 07. using codan to detect complete and partial CDS
bash scripts/07_codan.sh

# 99. misc stats and figures
# Rscript 99_stats_and_figures.R