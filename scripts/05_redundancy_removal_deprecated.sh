#!/bin/bash

# alorenzetti 202108 

# description ####
# this script will take the
# clusters from last step
# and perform redundancy
# removal by using cd-hit

# requires cd-hit from bioconda

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate cdhit_env

if [[ ! -d ./red_rem ]] ; then
    mkdir ./red_rem

    cd-hit-est -i ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster.hq.fasta.gz \
               -o ./red_rem/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster_hq_nr.fa \
               -c 0.99 -sc \
               -T $threads > ./red_rem/cd_hit_est.log 2> ./red_rem/cd_hit_est.err
fi