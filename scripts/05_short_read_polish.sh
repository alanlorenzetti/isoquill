#!/bin/bash

# alorenzetti 20211026

# description ####
# this script will take
# the FLNC clustered reads
# and polish using a short read
# dataset obtained from similar samples

# requires trimmomatic from bioconda: trimmomatic_env
# requires lordec from bioconda: lordec_env
# requires seqtk from bioconda: seqtk_env
# requires pigz: seqtk_env

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
# activating trimmomatic env
conda activate trimmomatic_env

if [[ ! -d short_read_polish ]] ; then
    # creating adapter fasta according to
    # novogene report
    echo -e ">fiveprime_adap\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" > short_read_polish/adap.fa
    echo -e ">threeprime_adap\nGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG" >> short_read_polish/adap.fa

    mkdir short_read_polish
fi

if [[ ! -d short_read_polish/trimmed ]] ; then
    mkdir short_read_polish/trimmed

    # the first step is to trim short-read
    # dataset to include only high quality reads

    # running trimmomatic for samples
    prefixes=`ls ../quillaja_illumina/raw/*.fq.gz | sed 's/_[12].*.fq.gz//' | sort | uniq`

    for prefix in $prefixes ; do
        sample=${prefix##*/}

        trimmomatic PE \
                    -threads $threads \
                    ${prefix}_1.fq.gz ${prefix}_2.fq.gz \
                    short_read_polish/trimmed/${sample}-paired_1.fq.gz short_read_polish/trimmed/${sample}-unpaired_1.fq.gz \
                    short_read_polish/trimmed/${sample}-paired_2.fq.gz short_read_polish/trimmed/${sample}-unpaired_2.fq.gz \
                    ILLUMINACLIP:short_read_polish/adap.fa:1:30:10 \
                    SLIDINGWINDOW:4:30 \
                    MINLEN:50 > short_read_polish/trimmed/${sample}_out.log 2> short_read_polish/trimmed/${sample}_err.log
    done
fi

# activating seqtk env
conda activate seqtk_env

# subsetting short reads is necessary
# we cannot use more than 128GB RAM
# using 1/10 of the reads
prefixes=`ls short_read_polish/trimmed/*.fq.gz | sed 's/-.*_[12].*.fq.gz//' | sort | uniq`
for prefix in $prefixes ; do 
    if [[ -d short_read_polish/trimmed ]] ; then 
        if [[ ! -f ${prefix}-paired-sub_1.fq.gz || ! -f ${prefix}-paired-sub_2.fq.gz ]] ; then 
            lines=`unpigz -p $threads -c ${prefix}-paired_1.fq.gz | wc -l`
            subsample=$((lines/4/10))

            seqtk sample -s 665 ${prefix}-paired_1.fq.gz $subsample | pigz -p $threads -c > ${prefix}-paired-sub_1.fq.gz
            seqtk sample -s 665 ${prefix}-paired_2.fq.gz $subsample | pigz -p $threads -c > ${prefix}-paired-sub_2.fq.gz
        fi

        if [[ ! -f ${prefix}-unpaired-sub_1.fq.gz ]] ; then 
            lines=`unpigz -p $threads -c ${prefix}-unpaired_1.fq.gz | wc -l`
            subsample=$((lines/4/10))

            seqtk sample -s 665 ${prefix}-unpaired_1.fq.gz $subsample | pigz -p $threads -c > ${prefix}-unpaired-sub_1.fq.gz
        fi

        if [[ ! -f ${prefix}-unpaired-sub_2.fq.gz ]] ; then 
            lines=`unpigz -p $threads -c ${prefix}-unpaired_2.fq.gz | wc -l`
            subsample=$((lines/4/10))

            seqtk sample -s 665 ${prefix}-unpaired_2.fq.gz $subsample | pigz -p $threads -c > ${prefix}-unpaired-sub_2.fq.gz
        fi
    fi
done

# activating lordec env
conda activate lordec_env

# the next step is to run fmlrc to correct
# clustered FLNC reads
if [[ ! -d short_read_polish/polished ]] ; then 
    mkdir short_read_polish/polished

    cd short_read_polish/polished

    trimmedreads=`ls ../trimmed/*-sub_[12].fq.gz | xargs`
          
    lordec-correct -T $threads \
                   -2 $trimmedreads \
                   -k 15 \
                   -s 3 \
                   -i ../../quillaja_bucket/quillaja_isoseq/isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster.hq.fasta.gz \
                   -o m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster_hq_polished.fasta \
                   > lordec_out.log 2> lordec_err.log

    cd ../..
fi