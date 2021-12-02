#!/bin/bash

# alorenzetti 202109

# description ####
# this script will take the hq
# clusters and align to a reference
# genome in order to perform an exploratory
# observation of isoforms detected and
# transcript quality (checking mismatches)

# requires minimap2 v2.22 from bioconda

# setting up number of threads
threads=6

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate minimap2_env

# minimap2 parameters were taken from
# https://github.com/Magdoll/cDNA_Cupcake/wiki/Best-practice-for-aligning-Iso-Seq-to-reference-genome:-minimap2,-deSALT,-GMAP,-STAR,-BLAT

# aligning

# downloading Qsap genome from
# NCBI Assembly GenBank version
if [[ ! -f data/Qsap.fa ]] ; then
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/338/715/GCA_003338715.1_Draft.Quillaja.v1.0/GCA_003338715.1_Draft.Quillaja.v1.0_genomic.fna.gz \
         -O data/Qsap.fa.gz

    gunzip -c data/Qsap.fa.gz > data/Qsap.fa
fi

# mapping transcripts against Qsap draft genome
if [[ ! -d cogent_vs_refGen_qsap ]] ; then
    mkdir cogent_vs_refGen_qsap

    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/Qsap.fa \
             cupcake/hq_transcripts.fasta.collapsed.rep.fa > cogent_vs_refGen_qsap/hq_cogent_vs_refGen_qsap.sam 2> cogent_vs_refGen_qsap/hq_cogent_vs_refGen_qsap_err.log

    samtools sort -@ $threads -o cogent_vs_refGen_qsap/hq_cogent_vs_refGen_qsap.bam cogent_vs_refGen_qsap/hq_cogent_vs_refGen_qsap.sam
    samtools index -b cogent_vs_refGen_qsap/hq_cogent_vs_refGen_qsap.bam
fi

# mapping transcripts against Qbra draft genome
if [[ ! -d cogent_vs_refGen_qbra ]] ; then
    mkdir cogent_vs_refGen_qbra

    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/final_gapClosed.fa \
             cupcake/hq_transcripts.fasta.collapsed.rep.fa > cogent_vs_refGen_qbra/hq_cogent_vs_refGen_qbra.sam 2> cogent_vs_refGen_qbra/hq_cogent_vs_refGen_qbra_err.log

    samtools sort -@ $threads -o cogent_vs_refGen_qbra/hq_cogent_vs_refGen_qbra.bam cogent_vs_refGen_qbra/hq_cogent_vs_refGen_qbra.sam
    samtools index -b cogent_vs_refGen_qbra/hq_cogent_vs_refGen_qbra.bam
fi
