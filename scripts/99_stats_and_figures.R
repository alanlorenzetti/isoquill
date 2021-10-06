# alorenzetti 202108

# description ####
# this companion script
# will get stats
# using the log files
# from pacbio isoseq
# processing

# setting working directory ####
setwd("~/gdrive/quillaja_transcriptome/isoquill/")

# loading libs ####
library(pacman)

packs = c("tidyverse",
          "readtext",
          "Biostrings",
          "ggpubr")

p_load(char = packs)

setwd("~/gdrive/quillaja_transcriptome/isoquill/")

# creating a directory to store plots
if(!dir.exists("plots")){dir.create("plots")}

# reading ccs report files
ccs = list()

ccs[["reports"]] = readtext(file = "ccs_12_parts/*.txt")

ccs[["stats"]] = ccs[["reports"]] %>% 
  mutate(text = str_replace_all(string = text, pattern = "\\n", replacement = " ")) %>% 
  mutate(text = str_replace_all(string = text, pattern = "\\s{2,}", replacement = " ")) %>% 
  mutate(totalzmw = str_replace_all(string = text,
                                    pattern = "ZMWs input : (\\d+).* ?",
                                    replacement = "\\1"),
         passzmw = str_replace_all(string = text,
                                   pattern = ".*ZMWs pass filters : (\\d+).* ?",
                                   replacement = "\\1"),
         failzmw = str_replace_all(string = text,
                                   pattern = ".*ZMWs fail filters : (\\d+).* ?",
                                   replacement = "\\1"),
         failsnr = str_replace_all(string = text,
                                   pattern = ".*Below SNR threshold : (\\d+).* ?",
                                   replacement = "\\1"),
         failpass = str_replace_all(string = text,
                                    pattern = ".*Lacking full passes : (\\d+).* ?",
                                    replacement = "\\1"),
         failrq = str_replace_all(string = text,
                                  pattern = ".*CCS below minimum RQ : (\\d+).* ?",
                                  replacement = "\\1")) %>% 
  as_tibble() %>% 
  select(-c(doc_id, text)) %>% 
  mutate(across(.cols = everything(),
                .fns = ~ as.numeric(.x)))

# getting N50 for collapsed transcripts
collapsed = readDNAStringSet(filepath = "cupcake/hq_transcripts.fasta.collapsed.rep.fa",
                             format = "fasta")
n50 = N50(width(collapsed))

# creating size histograms for different classes of
# transcripts according to codan
collapsedClassified = list()
collapsedClassified$Full = readDNAStringSet(filepath = "codan/codan_full_pred/ORF_sequences.fasta",
                                            format = "fasta")
collapsedClassified$Partial = readDNAStringSet(filepath = "codan/codan_partial_pred/ORF_sequences.fasta",
                                            format = "fasta")
collapsedClassified$NC = readDNAStringSet(filepath = "codan/non_full_non_partial_cds_transcript_list.fa",
                                            format = "fasta")

collapsedtib = tibble(name = names(collapsed) %>% str_replace(string = .,
                                                              pattern = "^(.*?)\\|.*$",
                                                              replacement = "\\1"),
                      width = width(collapsed),
                      seq = collapsed %>% as.character()) %>% 
  mutate(class = case_when(name %in% names(collapsedClassified$Full) ~ "Full",
                           name %in% names(collapsedClassified$Partial) ~ "Partial",
                           name %in% names(collapsedClassified$NC) ~ "Non-coding",
                           TRUE ~ NA_character_))

# plotting transcript lengths
# according to class
txlenplot = collapsedtib %>% 
  ggplot(aes(x = factor(class, levels = c("Full", "Partial", "Non-coding")), y = width)) +
  geom_violin(alpha = 0.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1),
               geom = "pointrange", color = "black") +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "Full") +
  ylab("Length") +
  xlab("Class")

ggsave(filename = "plots/txlenclass.png",
       plot = txlenplot,
       units = "in",
       dpi = 300,
       width = 3.5,
       height = 3.5)
