#!/bin/bash

# alorenzetti 20211005

# description ####
# this script will take as input
# the cogent redundancy reduced
# Quillaja brasiliensis transcriptome
# and predict full and partial CDS using CodAn

# requires a few dependencies from conda
# conda create -n codan_env python=3.7
# conda install biopython
# conda install -c bioconda perl perl-mce blast codan

# requires emboss infoseq
# requires blast
# requires R

# setting up number of threads
threads=6

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate codan_env

if [[ ! -d codan ]] ; then
    mkdir codan

    # adjusting input fasta headers
    perl -pe 's/^(>PB.*?)\|.*/\1/' \
         cupcake/hq_transcripts.fasta.collapsed.rep.fa > codan/collapsed_isoforms.fa

    cd codan

    # writing names of all transcripts
    grep ">" collapsed_isoforms.fa | \
    sed 's/>//' > transcript_list.txt

    # downloading and preparing codan models
    if [[ ! -d codan_models ]] ; then
        mkdir codan_models

        cd codan_models

        wget https://github.com/pedronachtigall/CodAn/raw/master/models/PLANTS_full.zip
        unzip PLANTS_full.zip

        wget https://github.com/pedronachtigall/CodAn/raw/master/models/PLANTS_partial.zip
        unzip PLANTS_partial.zip

        cd ..
    fi

    # running codan using full model
    if [[ ! -d codan_full_pred ]] ; then
        codan.py -s plus -c $threads -t collapsed_isoforms.fa -o codan_full_pred -m codan_models/PLANTS_full

        # writing names of full cds transcripts
        grep ">" codan_full_pred/ORF_sequences.fasta | \
        sed 's/>//' > full_cds_transcript_list.txt
    fi

    # running codan using partial model
    # for only those that have not been
    # classified as full
    if [[ ! -d codan_partial_pred ]] ; then
        # getting names of transcripts with no full CDS prediction
        R -s -e 'x = read.table(file = "full_cds_transcript_list.txt", header = F)$V1 ;
                 y = read.table(file = "transcript_list.txt", header = F)$V1 ;
                 fin = y[!(y %in% x)] ;
                 write.table(x = fin, file = "non_full_cds_transcript_list.txt", col.names = F, row.names = F, quote = F)'

        # creating a blast database to extract sequences
        makeblastdb -in collapsed_isoforms.fa -dbtype nucl -parse_seqids

        # extracting using blastdbcmd
        blastdbcmd -db collapsed_isoforms.fa -entry_batch non_full_cds_transcript_list.txt > non_full_cds_transcript_list.fa

        # running codan
        codan.py -s plus -c $threads -t non_full_cds_transcript_list.fa -o codan_partial_pred -m codan_models/PLANTS_partial

        # writing names of partial cds transcripts
        grep ">" codan_partial_pred/ORF_sequences.fasta | \
        sed 's/>//' > partial_cds_transcript_list.txt
    fi

    # extracting transcripts without full cds predictions
    # and without partial cds predictions
    # getting a list of transcripts non full non partial
    R -s -e 'x = read.table(file = "partial_cds_transcript_list.txt", header = F)$V1 ;
                y = read.table(file = "non_full_cds_transcript_list.txt", header = F)$V1 ;
                fin = y[!(y %in% x)] ;
                write.table(x = fin, file = "non_full_non_partial_cds_transcript_list.txt", col.names = F, row.names = F, quote = F)'

    # extracting using blastdbcmd
    blastdbcmd -db collapsed_isoforms.fa -entry_batch non_full_non_partial_cds_transcript_list.txt | \
    sed 's/ $//' > non_full_non_partial_cds_transcript_list.fa
    
    cd ..
fi