# alorenzetti 20211108

# description ####
# this script will take
# results from TRAPID
# annotation and incorporate
# them in this collection.
# they are gonna be required
# for functional enrichment analysis
# and so forth 

# getting started ####
# creating a tx2gene object ####
# for now, transcripts isoforms will be ignored
# the criterion for inclusion is:
# PB.1.1, PB.1.2, and PB.1.n will both comprise the gene PB.1
transcriptome = readDNAStringSet(filepath = "cupcake/hq_transcripts.fasta.collapsed.rep.fa",
                                 format = "fasta")

transcriptomeTbl = tibble(transcript_id = names(transcriptome) %>% str_replace("\\|.*$", ""),
                          transcript_nt = as.character(transcriptome))

tx2gene = tibble(TXNAME = names(transcriptome) %>% str_replace("\\|.*$", ""),
                 GENEID = names(transcriptome) %>% str_replace("\\.[0-9]{1,}\\|.*$", ""))

# TRAPID ####
# parsing files downloaded from
# TRAPID Plaza 4.5 dicots
trapid = list()

# reading interpro
trapid$interpro = read_tsv(file = "data/transcripts_interpro_exp4325.txt") %>% 
  dplyr::select(-`#counter`) %>% 
  dplyr::rename(interpro_id = "interpro",
                interpro_description = "description")

# reading gene ontology
trapid$go = read_tsv(file = "data/transcripts_go_exp4325.txt") %>% 
  dplyr::select(-`#counter`,-evidence_code,-is_hidden) %>%
  dplyr::rename(go_id = "go",
                go_description = "description")

# reading gene families
trapid$gf = read_tsv(file = "data/transcripts_gf_exp4325.txt") %>% 
  dplyr::select(-`#counter`)

# reading structure (full length info)
trapid$structural = read_tsv(file = "data/structural_data_exp4325.txt") %>% 
  dplyr::rename(transcript_id = "#transcript_id",
                trapid_status = "meta_annotation") %>% 
  dplyr::select(transcript_id,
                trapid_status) %>% 
  mutate(trapid_status = case_when(trapid_status == "No Information" ~ "Not full/Not Partial",
                                   trapid_status == "Quasi Full Length" ~ "Partial",
                                   TRUE ~ trapid_status))

# parsing go codes
gocodes = trapid$go %>%
  dplyr::select(go_id, go_description) %>% 
  distinct()

# parsing interpro codes
ipcodes = trapid$interpro %>%
  dplyr::select(interpro_id, interpro_description) %>% 
  distinct()

# Codan ####
codan = list()

codan$full = tibble(transcript_id = readDNAStringSet(filepath = "codan/codan_full_pred/ORF_sequences.fasta",
                                                     format = "fasta") %>% 
                      names(),
                    codan_status = "Full Length",
                    codan_nt = readDNAStringSet(filepath = "codan/codan_full_pred/ORF_sequences.fasta",
                                                format = "fasta") %>% 
                      as.character(),
                    codan_aa = readAAStringSet(filepath = "codan/codan_full_pred/PEP_sequences.fa",
                                               format = "fasta") %>% 
                      as.character())

codan$partial = tibble(transcript_id = readDNAStringSet(filepath = "codan/codan_partial_pred/ORF_sequences.fasta",
                                                        format = "fasta") %>% 
                         names(),
                       codan_status = "Partial",
                       codan_nt = readDNAStringSet(filepath = "codan/codan_partial_pred/ORF_sequences.fasta",
                                                   format = "fasta") %>% 
                         as.character(),
                       codan_aa = NA_character_)

codan$nonFullnonPartial = tibble(transcript_id = readDNAStringSet(filepath = "codan/non_full_non_partial_cds_transcript_list.fa",
                                                                  format = "fasta") %>% 
                                   names(),
                                 codan_status = "Not full/Not Partial",
                                 codan_nt = NA_character_,
                                 codan_aa = NA_character_)

codan$all = bind_rows(codan$full,
                      codan$partial,
                      codan$nonFullnonPartial)

# Kraken ####
kraken = list()

kraken$all = read_tsv(file = "decon/decon_tx/hq_transcripts_fasta_collapsed_rep_output.txt",
                      col_names = F) %>% 
  dplyr::rename(class = "X1",
                transcript_id = "X2",
                kraken_taxid = "X3") %>% 
  filter(class == "C") %>% 
  mutate(transcript_id = str_replace(transcript_id, "^(.*?)\\|.*$", "\\1")) %>% 
  dplyr::select(transcript_id,
                kraken_taxid)

# if interested in finding out
# what these taxids correspond to
# (needs further improvement since some taxids are returning NA)
# preparing ncbi taxid database
if(!file.exists("data/accessionTaxa.sql")){
  getNamesAndNodes(outDir = "data")
  read.names.sql(nameFile = "data/names.dmp", sqlFile = "data/accessionTaxa.sql", overwrite = T)
  read.nodes.sql(nodeFile = "data/nodes.dmp", sqlFile = "data/accessionTaxa.sql", overwrite = T)
}

# getting taxonomy for kraken taxids
classifiedTX = getTaxonomy(ids = kraken$all$kraken_taxid,
                           sqlFile = "data/accessionTaxa.sql") %>% 
  as_tibble(rownames = "kraken_taxid") %>% 
  mutate(kraken_taxid = str_replace(string = kraken_taxid,
                                    pattern = " {1,}",
                                    replacement = ""))

minimizeNAtaxa = function(df = df){
  newdf = df %>% 
    filter(!is.na(species)) %>% 
    dplyr::select(kraken_taxid,
                  kraken_taxon = species)
  
  df = df %>% 
    filter(is.na(species))
  
  taxa = c("species",
           "genus",
           "family",
           "order",
           "class",
           "phylum",
           "superkingdom")
  
  for(i in 1:dim(df)[1]){
    for(j in taxa){
      if(is.na(df[i,j]) %>% as.logical()){
        if(j == "superkingdom"){
          kraken_taxid = df[i,"kraken_taxid"] %>% unlist()
          kraken_taxon = "Classified without taxon"
          
          newdf = bind_rows(newdf,
                            tibble(kraken_taxid = kraken_taxid,
                                   kraken_taxon = kraken_taxon))
        }
        next
      }
      else{
        kraken_taxid = df[i,"kraken_taxid"] %>% unlist()
        kraken_taxon = df[i,j] %>% unlist()
        
        newdf = bind_rows(newdf,
                          tibble(kraken_taxid = kraken_taxid,
                                 kraken_taxon = kraken_taxon))
        break
      }
    }
  }
  
  newdf = newdf %>%
    mutate(kraken_taxid = as.numeric(kraken_taxid)) %>% 
    distinct()
  
  return(newdf)
}

# minimizing number of unknown classified transcripts
kraken$taxid2taxa = minimizeNAtaxa(df = classifiedTX) 

# unifying taxon names with taxid obj
kraken$all = left_join(kraken$all,
                       kraken$taxid2taxa,
                       by = "kraken_taxid")

kraken$final = left_join(tx2gene, 
                         kraken$all,
                         by = c("TXNAME" = "transcript_id")) %>% 
  mutate(kraken_taxon = case_when(is.na(kraken_taxon) ~ "Unclassified",
                                  TRUE ~ kraken_taxon)) %>% 
  dplyr::select(transcript_id = "TXNAME",
                kraken_taxid,
                kraken_taxon)

# tps and p450 prediction ####
tps_p450 = list()

tps_p450$tps$mono = read_tsv(file = "find_tps_p450/find_tps/monoTP_dir/high_confidence_final_results.csv") %>% 
  dplyr::rename(transcript_id = "Sequence",
                score = "Highest score") %>% 
  mutate(tps_type = "Monoterpene")

tps_p450$tps$di = read_tsv(file = "find_tps_p450/find_tps/diTP_dir/high_confidence_final_results.csv") %>% 
  dplyr::rename(transcript_id = "Sequence",
                score = "Highest score") %>% 
  mutate(tps_type = "Diterpene")

tps_p450$tps$sesqui = read_tsv(file = "find_tps_p450/find_tps/sesquiTP_dir/high_confidence_final_results.csv") %>% 
  dplyr::rename(transcript_id = "Sequence",
                score = "Highest score") %>% 
  mutate(tps_type = "Sesquiterpene",
         transcript_id = as.character(transcript_id),
         score = as.numeric(score))

tps_p450$tps$all = bind_rows(tps_p450$tps$mono,
                             tps_p450$tps$di,
                             tps_p450$tps$sesqui) %>% 
  mutate(tps_p450 = "tps")

tps_p450$p450 = read_table(file = "find_tps_p450/find_p450/p450_table.txt",
                           comment = "#",
                           col_names = F) %>% 
  dplyr::rename(transcript_id = "X1",
                score = "X6") %>% 
  dplyr::select(transcript_id,
                score) %>% 
  mutate(tps_p450 = "p450")

tps_p450$all = bind_rows(tps_p450$tps$all,
                         tps_p450$p450)

# reference genome alignment ####
qsapRefGen = read_tsv(file = "cogent_vs_refGen_qsap/hq_cogent_vs_refGen_qsap.sam",
                      comment = "@",
                      col_names = F) %>% 
  dplyr::select(X1,X2) %>% 
  filter(X2 != 4) %>% 
  pull(X1) %>% 
  str_replace("\\|.*$", "") %>% 
  unique() %>% 
  tibble(transcript_id = .,
         qsap_ref_aln = "yes")

qbraRefGen = read_tsv(file = "cogent_vs_refGen_qbra/hq_cogent_vs_refGen_qbra.sam",
                      comment = "@",
                      col_names = F) %>% 
  dplyr::select(X1,X2) %>% 
  filter(X2 != 4) %>% 
  pull(X1) %>% 
  str_replace("\\|.*$", "") %>% 
  unique() %>% 
  tibble(transcript_id = .,
         qbra_ref_aln = "yes")

# unifying everything ####
funcat = tx2gene %>% 
  dplyr::rename(transcript_id = "TXNAME",
                gene_id = "GENEID")

# completeness
# adding trapid completeness
funcat = left_join(funcat, trapid$structural,
                   by = "transcript_id")

# adding codan completeness
funcat = left_join(funcat, codan$all,
                   by = "transcript_id")

# adding trapid gene families
funcat = left_join(funcat,
                   trapid$gf,
                   by = "transcript_id")

# adding trapid interpro
funcat = left_join(funcat,
                   trapid$interpro %>% 
                     group_by(transcript_id) %>% 
                     summarise(interpro_id = paste0(interpro_id, collapse = ";"),
                               interpro_description = paste0(interpro_description, collapse = ";")),
                   by = "transcript_id")

# adding trapid go
funcat = left_join(funcat,
                   trapid$go %>% 
                     group_by(transcript_id) %>% 
                     summarise(go_id = paste0(go_id, collapse = ";"),
                               go_description = paste0(go_description, collapse = ";")),
                   by = "transcript_id")

# adding kraken classification
funcat = left_join(funcat,
                   kraken$final,
                   by = "transcript_id")

# adding tps and p450 info
funcat = left_join(funcat,
                   tps_p450$all,
                   by = "transcript_id")

# adding transcript sequence
funcat = left_join(funcat,
                   transcriptomeTbl,
                   by = "transcript_id")

# adding ref genome alignment info
# qsap
funcat = left_join(funcat,
                   qsapRefGen,
                   by = "transcript_id")

# qbra
funcat = left_join(funcat,
                   qbraRefGen,
                   by = "transcript_id")

# adjusting NAs
funcat = funcat %>% 
  dplyr::mutate(across(.cols = ends_with("ref_aln"),
                       .fns = ~ {case_when(is.na(.x) ~ "no",
                                           TRUE ~ .x)}))

# making two versions of funcat obj
# one of them without sequences
funcatSeq = funcat
funcat = funcat %>% 
  dplyr::select(-codan_aa, -codan_nt, -transcript_nt)
