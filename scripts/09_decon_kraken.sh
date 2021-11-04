#!/bin/bash

# alorenzetti 202111

# description ####
# this script will take
# the qbra transcriptome
# and perform decontamination
# using dbs of known microorganisms
# and viruses
# contamination is not likely
# but we wanna make sure our final
# transcriptome dataset is purely
# comprised by plant transcripts

# requires kraken2 from bioconda

# setting up num threads
threads=14

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate kraken2_env

# getting and creating dbs
if [[ ! -d decon ]] ; then mkdir decon ; fi

if [[ ! -d decon/dbs ]] ; then
    mkdir decon/dbs

    cd decon/dbs

    # downloading ready to use kraken index
    mkdir k2_pluspf_20210517
    cd k2_pluspf_20210517
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210517.tar.gz
    tar -zxvf k2_pluspf_20210517.tar.gz
    cd ../../..
fi

if [[ ! -d decon/decon_tx ]] ; then
    mkdir decon/decon_tx

    # going down to processing dir
    cd decon/decon_tx

    kraken2 --db ../dbs/k2_pluspf_20210517 --threads $threads \
            --unclassified-out hq_transcripts_fasta_collapsed_rep_unclassified.fa \
            --classified-out hq_transcripts_fasta_collapsed_rep_classified.fa \
            --output hq_transcripts_fasta_collapsed_rep_output.txt \
            --report hq_transcripts_fasta_collapsed_rep_report.txt \
            /home/alorenzetti/quillaja_bucket/quillaja_isoseq/cupcake/hq_transcripts.fasta.collapsed.rep.fa

    cd ../..
fi