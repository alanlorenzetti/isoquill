#!/bin/bash

# alorenzetti 202112

# description ####
# this script will run cogent
# in order to collapse clusters
# output by isoseq3

# requires cogent and its dependencies
# cogent should be installed following 
# https://github.com/Magdoll/Cogent/wiki/Installing-Cogent

# requires GNU parallel in the same env
# conda install -c conda-forge parallel

# setting up num threads
threads=15

# initializing conda 
eval "$(conda shell.bash hook)"

# getting started ####
conda activate cogent_env

# getting started 
# https://github.com/Magdoll/Cogent/wiki/Running-Cogent#findbig

# creating a directory to perform the tasks
if [[ ! -d fam_find ]] ; then
    mkdir fam_find

    # entering the dir
    cd fam_find

    # copying cluster output file and gunzipping
    cp /home/alorenzetti/quillaja_bucket/quillaja_isoseq/short_read_polish/polished/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster_hq_polished.fasta .

    # setting the appropriate file name
    mv m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster_hq_polished.fasta isoseq_flnc.fasta

    # generating bins
    run_preCluster.py --cpus=$threads > run_preCluster_log.log 2> run_preCluster_err.log

    # generating batch commands to run family finding
    generate_batch_cmd_for_Cogent_family_finding.py --cpus=$threads --cmd_filename=cmd preCluster.cluster_info.csv preCluster_out quillaja_fam

    # running automatically generated bash script
    bash cmd > cmd_log.log 2> cmd_err.log

    # getting results
    printf "Partition\tSize\tMembers\n" > final.partition.txt
    ls preCluster_out/*/*partition.txt | xargs -n1 -i sed '1d; $d' {} | cat >> final.partition.txt

    # remove preclustering stuff, don't need them anymore
    rm -rf preCluster_out/
    rm -rf preCluster.*

    # performing coding genome reconstruction
    # generating batch commands
    # generate_batch_cmd_for_Cogent_reconstruction.py quillaja_fam > cmd_reconstruct.txt

    # # getting a list of commands to run in parallel
    # sed -i 's/reconstruct_contig.py //' cmd_reconstruct.txt
    ls quillaja_fam > cmd_reconstruct.txt

    # running in parallel 
    # if [[ ! -d reconstruct_log ]] ; then mkdir reconstruct_log ; fi
    parallel -a cmd_reconstruct.txt --jobs $threads reconstruct_contig.py quillaja_fam/{} -p {} > reconstruct_log.log 2>&1

    # getting a summary of results
    python /home/alorenzetti/compiledPrograms/Cogent/Cogent/helper_scripts/tally_Cogent_contigs_per_family.py quillaja_fam quillaja quillaja_output

    # creating a fake genome for the next step
    cat quillaja_fam/*/cogent2.renamed.fasta > cogent.fake_genome.fasta

    # going back to parent dir
    cd ..
fi

# collapsing transcripts in absence of a genome
# https://github.com/Magdoll/Cogent/wiki/Tutorial%3A-Using-Cogent-to-collapse-redundant-transcripts-in-absence-of-genome
if [[ ! -d cupcake ]] ; then
    # creating a new directory for the results
    mkdir cupcake
    cd cupcake

    # copying fake genome to current directory
    cp ../fam_find/cogent.fake_genome.fasta .

    # collapsing redundant isoforms
    minimap2 -ax splice -t $threads -uf --secondary=no \
             cogent.fake_genome.fasta ../fam_find/isoseq_flnc.fasta > hq_transcripts.fasta.sam 2> clusters_vs_fakegenome_err.log

    # collapsing from sam
    sort -k 3,3 -k 4,4n hq_transcripts.fasta.sam > hq_transcripts.fasta.sorted.sam

    collapse_isoforms_by_sam.py --input ../fam_find/isoseq_flnc.fasta -s hq_transcripts.fasta.sorted.sam \
                                -c 0.95 -i 0.85 \
                                -o hq_transcripts.fasta > collapse_isoforms_by_sam_log.log 2> collapse_isoforms_by_sam_err.log

    # getting abundance 
    # this step needs the cluster_report.csv
    # file output by isoseq3 cluster 
    get_abundance_post_collapse.py hq_transcripts.fasta.collapsed /home/alorenzetti/quillaja_bucket/quillaja_isoseq/isoseq_refine_cluster/m64168e_210807_154604_ccs_lima_refine_onlyPolyA_cluster.cluster_report.csv > get_abundance_post_collapse_log.log 2> get_abundance_post_collapse_err.log

    # filtering to obtain the final file
    #filter_away_subset.py hq_transcripts.fasta.collapsed > filter_away_subset_log.log 2> filter_away_subset_err.log

    # going back to parent dir
    cd ..
fi
