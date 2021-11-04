#!/bin/bash

# alorenzetti 20211103

# description ####
# this script will take
# the qbra transcriptome and
# will search for terpene synthases
# using search_TPS
# https://github.com/liliane-sntn/TPS
# and cytochrome p450s
# using the appropriate PFAM hmm

# requires hmmer from bioconda
# requires git

# setting up num threads
threads=14

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate hmmer_env

# creating analysis dir
if [[ ! -d find_tps_p450 ]] ; then mkdir find_tps_p450 ; fi

# downloading required software
# and unzipping required dbs
if [[ ! -d find_tps_p450/TPS ]] ; then
    cd find_tps_p450

    git clone https://github.com/liliane-sntn/TPS

    cd TPS/Tutorial

    unzip class-specific_hmms.zip
    unzip score_tables_dir.zip
    unzip PFAMs_dir.zip

    cd ../../..
fi

# searching for TPS
if [[ ! -d find_tps_p450/find_tps ]] ; then
    cd find_tps_p450

    perl TPS/search_TPS.pl -i /home/alorenzetti/quillaja_bucket/quillaja_isoseq/codan/codan_full_pred/PEP_sequences.fa \
                           -d TPS/Tutorial/class-specific_hmms \
                           -p TPS/Tutorial/PFAMs_dir \
                           -t TPS/Tutorial/score_tables_dir/all_scores.csv \
                           -o find_tps \
                           -s 1
    
    cd ..
fi

# searching for P450
if [[ ! -d find_tps_p450/find_p450 ]] ; then
    mkdir find_tps_p450/find_p450
    
    cd find_tps_p450/find_p450

    wget http://pfam.xfam.org/family/PF00067/hmm -O p450.hmm

    hmmsearch --tblout p450_table.txt \
              -o p450_output.txt \
              --cpu $threads p450.hmm /home/alorenzetti/quillaja_bucket/quillaja_isoseq/codan/codan_full_pred/PEP_sequences.fa > hmmsearch_out.log 2> hmmsearch_err.log

    cd ../..
fi