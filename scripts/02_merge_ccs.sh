#!/bin/bash

# alorenzetti 202108 
# description ####
# this script will merge partitioned
# ccs files generated separately

# requires pbbam from bioconda
# requires pbccs from bioconda

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate isoseqv3_env

if [[ ! -d ./ccs ]] ; then
    mkdir ./ccs

    pbmerge -o ./ccs/m64168e_210807_154604_ccs.bam ./ccs_12_parts/m64168e_210807_154604_ccs_*.bam > ./ccs/pbmerge.log 2> ./ccs/pbmerge.err
    pbindex ./ccs/m64168e_210807_154604_ccs.bam > ./ccs/pbindex.log 2> ./ccs/pbindex.err
fi