#!/bin/bash

# alorenzetti 202108

# description ####
# this script will run busco
# for sequences after collapsing isoforms
# the purpose of this script is
# to assess the completeness of the transcriptome
# using the set of single copy orthologous gene
# shared among embryophyta

# requires busco v5.2.2 from bioconda
# requires perl

# setting up num threads
threads=14

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate busco_env

if [[ ! -d ./busco ]] ; then
   mkdir ./busco

    # running busco for collapsed version of transcriptome (cogent output)
    # busco does not like the slash character on fasta headers
    # so I am gonna replace it by something else
    perl -pe 's/^(>PB.*?)\|.*/\1/' \
         /home/alorenzetti/quillaja_bucket/quillaja_isoseq/cupcake/hq_transcripts.fasta.collapsed.rep.fa > ./busco/collapsed_isoforms.fa

    busco -i ./busco/collapsed_isoforms.fa \
          -o collapsed_isoforms \
          --out_path ./busco \
          --download_path ./busco/download \
          -l embryophyta_odb10 \
          -m tran \
          -c $threads > ./busco/collapsed_isoforms_log.log 2> ./busco/collapsed_isoforms_err.log
fi
