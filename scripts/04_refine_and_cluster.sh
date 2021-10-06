#!/bin/bash

# alorenzetti 202108 

# description ####
# this script will take the
# ccs trimmed files and remove
# poly-a and concatamers by
# using isoseq refine bringing up
# flnc tx (full-length non-concatamer)
# then output files
# will be submitted to 
# the clustering step

# requires isoseq3 from bioconda

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

if [[ ! -d ./isoseq_refine_cluster ]] ; then
    mkdir ./isoseq_refine_cluster

    # isoseq3 refine -j $threads \
    #                ./lima/m64168e_210807_154604_ccs_lima.NEB_5p--NEB_3p.bam \
    #                ./lima/adaps.fa \
    #                ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine.bam > ./isoseq_refine_cluster/refine.log 2> ./isoseq_refine_cluster/refine.err

    isoseq3 refine -j $threads --require-polya \
                   ./lima/m64168e_210807_154604_ccs_lima.NEB_5p--NEB_3p.bam \
                   ./lima/adaps.fa \
                   ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA.bam > ./isoseq_refine_cluster/refine_onlyPolyA.log 2> ./isoseq_refine_cluster/refine_onlyPolyA.err

    # isoseq3 cluster -j $threads --use-qvs --verbose \
    #                 ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine.bam \
    #                 ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_cluster.bam > ./isoseq_refine_cluster/cluster.log 2> ./isoseq_refine_cluster/cluster.err

    isoseq3 cluster -j $threads --use-qvs --verbose \
                    ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA.bam \
                    ./isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster.bam > ./isoseq_refine_cluster/cluster_onlyPolyA.log 2> ./isoseq_refine_cluster/cluster_onlyPolyA.err

fi