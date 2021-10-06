#!/bin/bash

# alorenzetti 202108 

# description ####
# this script will take the
# ccs file and trim the adapters
# and set up the true orientation
# of transcripts using lima

# requires lima from bioconda

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

# requires lima from bioconda
if [[ ! -d ./lima ]] ; then
    mkdir ./lima

    # observing the ccs reads I was able to detect
    # the adaps are those from NEB reported in
    # https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md
    # I will create the required files
    echo -e ">NEB_5p\nGCAATGAAGTCGCAGGGTTGGGG\n>NEB_3p\nGTACTCTGCGTTGATACCACTGCTT" > ./lima/adaps.fa

    # performing lima
    lima ./ccs/m64168e_210807_154604_ccs.bam \
         ./lima/adaps.fa \
         ./lima/m64168e_210807_154604_ccs_lima.bam \
         --isoseq --peek-guess -j $threads > ./lima/lima.log 2> ./lima/lima.err

fi
