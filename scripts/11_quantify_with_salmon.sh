#!/bin/bash

# alorenzetti 20211104

# description ####
# this script will
# use salmon to generate bootstrapped
# quantification of transcriptomes

# requires salmon from bioconda
# conda install -c conda-forge boost-cpp=1.74.0
# conda install -c bioconda salmon=1.5.2

# setting up number of threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate salmon_env

# creating a new dir
if [[ ! -d salmon_quant ]] ; then
    mkdir salmon_quant
    
    cp /home/alorenzetti/quillaja_bucket/quillaja_isoseq/cupcake/hq_transcripts.fasta.collapsed.rep.fa salmon_quant/qbrasiliensis_transcriptome.fa
fi

# making index for salmon
# from ncbi refseq
if [[ ! -d salmon_quant/Qbrasiliensis ]] ; then
    # making ref transcriptome index
    salmon index -p $threads -t salmon_quant/qbrasiliensis_transcriptome.fa -i salmon_quant/Qbrasiliensis > salmon_quant/salmon_index_out.log 2> salmon_quant/salmon_index_err.log
fi

# salmon quant
if [[ ! -d salmon_quant/libs ]] ; then
    mkdir -p salmon_quant/libs

    # running salmon quant for each sample
    for i in /home/alorenzetti/quillaja_bucket/quillaja_illumina/raw/*.fq.gz ; do
        prefix=${i%%_[12]*}
        libname=${prefix##*/}

            if [[ ! -d salmon_quant/libs/${libname} ]] ; then
                salmon quant -i salmon_quant/Qbrasiliensis \
                            -l A \
                            -1 ${prefix}_1.fq.gz \
                            -2 ${prefix}_2.fq.gz \
                            -p $threads \
                            -o salmon_quant/libs/${libname} \
                            --softclip \
                            --numBootstraps 20 > salmon_quant/libs/${libname}.log 2>&1
            fi
    done

    # getting the names of all output dirs
    outdirs=`ls /home/alorenzetti/quillaja_bucket/quillaja_illumina/raw/*.fq.gz | sed 's/.*\///g;s/_.*//g;s/^/salmon_quant\/libs\//g' | sort | uniq | xargs`

    # merging results from samples
    salmon quantmerge --quants $outdirs -o salmon_quant/libs/libs_merged_tpm.sf -c tpm
    salmon quantmerge --quants $outdirs -o salmon_quant/libs/libs_merged_numreads.sf -c numreads
fi