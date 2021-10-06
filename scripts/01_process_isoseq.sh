#!/bin/bash

# alorenzetti 202108 
# description ####
# this script will perform the isoseq v3 pipeline
# according to
# https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md

# requires pbccs from bioconda

# usage 
# bash process_isoseq.sh yes 1 4

# setting vars ####
threads=15

# getting args
# use chunking
chunking=$1

# which chunk
chunk=$2

# total chunks
totchunks=$3

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

# creating directory to store results
if [[ $1 == "yes" ]] ; then
    if [[ ! -d ./ccs_${3}_parts ]] ; then mkdir ./ccs_${3}_parts ; fi

    if [[ ! -f ./ccs_${3}_parts/m64168e_210807_154604_ccs_${2}.bam ]] ; then
        ccs \
        raw/m64168e_210807_154604.subreads.bam \
        ./ccs_${3}_parts/m64168e_210807_154604_ccs_${2}.bam \
        --chunk ${2}/${3} -j $threads 
    fi

else
    if [[ ! -d ./ccs ]] ; then
        mkdir ./ccs

        ccs -j $threads \
        raw/m64168e_210807_154604.subreads.bam \
        ./ccs/m64168e_210807_154604_ccs.bam
    fi
fi

