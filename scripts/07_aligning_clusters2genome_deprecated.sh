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
if [[ ! -d clusters_vs_refGen ]] ; then

    mkdir clusters_vs_refGen
    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/final_gapClosed.fa \
             isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster.hq.fasta.gz > clusters_vs_refGen/hq_clusters_vs_refGen.sam 2> clusters_vs_refGen/hq_clusters_vs_refGen_err.log

    samtools sort -@ $threads -o clusters_vs_refGen/hq_clusters_vs_refGen.bam clusters_vs_refGen/hq_clusters_vs_refGen.sam
    samtools index -b clusters_vs_refGen/hq_clusters_vs_refGen.bam

fi 

if [[ ! -d remred_vs_refGen ]] ; then
    mkdir remred_vs_refGen
    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/final_gapClosed.fa \
             rem_red/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster_hq_nr.fa > remred_vs_refGen/hq_remred_vs_refGen.sam 2> remred_vs_refGen/hq_remred_vs_refGen_err.log

    samtools sort -@ $threads -o remred_vs_refGen/hq_remred_vs_refGen.bam remred_vs_refGen/hq_remred_vs_refGen.sam
    samtools index -b remred_vs_refGen/hq_remred_vs_refGen.bam
fi

# mapping filtered transcripts against draft genome
if [[ ! -d cogent_vs_refGen ]] ; then
    mkdir cogent_vs_refGen
    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/final_gapClosed.fa \
             cupcake/hq_transcripts.fasta.no5merge.collapsed.filtered.rep.fa > cogent_vs_refGen/hq_cogentFiltRep_vs_refGen.sam 2> cogent_vs_refGen/hq_cogentFiltRep_vs_refGen_err.log

    samtools sort -@ $threads -o cogent_vs_refGen/hq_cogentFiltRep_vs_refGen.bam cogent_vs_refGen/hq_cogentFiltRep_vs_refGen.sam
    samtools index -b cogent_vs_refGen/hq_cogentFiltRep_vs_refGen.bam

    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/final_gapClosed.fa \
             cupcake/hq_transcripts.fasta.no5merge.collapsed.rep.fa > cogent_vs_refGen/hq_cogentRep_vs_refGen.sam 2> cogent_vs_refGen/hq_cogentRep_vs_refGen_err.log

    samtools sort -@ $threads -o cogent_vs_refGen/hq_cogentRep_vs_refGen.bam cogent_vs_refGen/hq_cogentRep_vs_refGen.sam
    samtools index -b cogent_vs_refGen/hq_cogentRep_vs_refGen.bam

    minimap2 -a -x splice -uf -C5 -O6,24 -B4 --secondary=no \
             -t $threads \
             data/final_gapClosed.fa \
             cupcake_v2/hq_transcripts.fasta.collapsed.rep.fa > cogent_vs_refGen/hq_cogent_vs_refGen.sam 2> cogent_vs_refGen/hq_cogent_vs_refGen_err.log

    samtools sort -@ $threads -o cogent_vs_refGen/hq_cogent_vs_refGen.bam cogent_vs_refGen/hq_cogent_vs_refGen.sam
    samtools index -b cogent_vs_refGen/hq_cogent_vs_refGen.bam
fi