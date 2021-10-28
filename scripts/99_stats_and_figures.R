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
          "ggpubr",
          "ggtext",
          "ggExtra",
          "ggthemes")

p_load(char = packs)

# setting up ggplot theme
theme_set(theme_bw())

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

# loading tx abundance info
abund = read_tsv(file = "cupcake/hq_transcripts.fasta.collapsed.abundance.txt",
                 comment = "#")

# joining abundance info
collapsedtib = left_join(x = collapsedtib,
                         y = abund,
                         by = c("name" = "pbid"))

# plotting transcript lengths
# according to class
compar = list(c("Full", "Partial"),
              c("Full", "Non-coding"),
              c("Partial", "Non-coding"))

txlenplot = collapsedtib %>% 
  ggplot(aes(x = factor(class, levels = c("Full", "Partial", "Non-coding")), y = width)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1),
               geom = "pointrange", color = "black") +
  stat_compare_means(comparisons = compar,
                     label = "p.signif",
                     method = "t.test") +
  ylab("Length") +
  xlab("Class")

# plotting abundance
# according to class
txabundplot = collapsedtib %>% 
  ggplot(aes(x = factor(class, levels = c("Full", "Partial", "Non-coding")), y = count_fl %>% log10())) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  stat_summary(fun.data = "mean_cl_normal", fun.args = list(mult = 1),
               geom = "pointrange", color = "black") +
  stat_compare_means(comparisons = compar,
                     label = "p.signif",
                     method = "wilcox.test") +
  ylab("Log<sub>10</sub>(Abundance)") +
  xlab("Class") +
  theme(axis.title.y = element_markdown())

lenabscat = collapsedtib %>% 
  ggplot(aes(color = factor(class, levels = c("Full", "Partial", "Non-coding")),
             x = count_fl %>% log10(),
             y = width)) +
  geom_point(alpha = 0.1, show.legend = T, stroke = 0, size = 3) +
  geom_smooth(method = "lm") +
  stat_cor(method = "pearson") +
  ylab("Length") +
  xlab("Log<sub>10</sub>(Abundance)") +
  scale_color_tableau() +
  labs(color = "Class") +
  # facet_grid(~ factor(class, levels = c("Full", "Partial", "Non-coding"))) +
  theme(axis.title.x = element_markdown(),
        legend.position = "bottom")

lenabscat = ggMarginal(lenabscat, type = "density", groupColour = T)

panel = ggarrange(plotlist = list(lenabscat,
                          ggarrange(plotlist = list(txlenplot,
                                                    txabundplot),
                                    labels = c("B", "C"))),
          nrow = 2,
          labels = c("A", NA))

ggsave(filename = "plots/classes_abund_len_panel.png",
       plot = panel,
       units = "in",
       dpi = 300,
       width = 5.5,
       height = 8)
