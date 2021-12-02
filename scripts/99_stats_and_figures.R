# alorenzetti 202108

# description ####
# this companion script
# will get stats
# using the log files
# from pacbio isoseq
# processing

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
  dplyr::select(-c(doc_id, text)) %>% 
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

# getting insights on number of genes ####
# and isoforms
genesIsos = funcat %>%
  dplyr::select(gene_id) %>%
  group_by(gene_id) %>%
  summarise(isoforms = n()) %>% 
  ungroup() %>% 
  group_by(isoforms) %>% 
  summarise(count = n())

# how many 10+?
tenplus = genesIsos %>% 
  filter(isoforms >= 10) %>% 
  pull(count) %>% 
  sum()

# combining
genesIsosSub = bind_rows(genesIsos %>%
                           filter(isoforms < 10) %>% 
                           mutate(isoforms = factor(isoforms)),
                         tibble(isoforms = factor("10+"),
                                count = tenplus))

# plotting 
geneIsosPlot = genesIsosSub %>% 
  ggplot(aes(x = factor(isoforms),
             y = count)) +
  geom_col(fill = "white",
           color = "black") +
  xlab("Number of isoforms per gene") +
  ylab("Absolute frequency")

# saving 
ggsave(filename = "plots/isoform_frequency_per_gene.png",
       plot = geneIsosPlot,
       device = "png",
       units = "in",
       dpi = 300,
       width = 4,
       height = 3)

# exploratory analysis panel
# generating a heat map to check
# sample grouping
# defining heatmap col funct
colfunct = circlize::colorRamp2(breaks = c(5,10,20), colors = viridis(3))

# creating matrix
M = assay(vsd)
colnames(M) = colData(vsd)$sample

# getting the 1000
# with bigger variance
select = order(rowVars(assay(vsd)),
               decreasing=TRUE)[1:1000]

# defining annotation colors
legendCols = list(
  tissue = c("leaf" = "#4E79A7",
             "cell_culture" = "#59A14F"),
  treatment = c("regular_light" = adjustcolor(col = "#B07AA1", alpha.f = 0.5),
                "uv_light" = adjustcolor(col = "#B07AA1", alpha.f = 1),
                "lag_phase" = adjustcolor(col = "#76B7B2", alpha.f = 0.2),
                "log_phase" = adjustcolor(col = "#76B7B2", alpha.f = 0.5),
                "stationary_phase" = adjustcolor(col = "#76B7B2", alpha.f = 1)),
  replicate = c("A" = adjustcolor(col = "#BAB0AC", alpha.f = 0.2),
                "B" = adjustcolor(col = "#BAB0AC", alpha.f = 0.5),
                "C" = adjustcolor(col = "#BAB0AC", alpha.f = 1))
)

# defining annotations
tiannot = colData(vsd)$tissue
treatannot = colData(vsd)$treatment
repanoot = colData(vsd)$replicate

annot = HeatmapAnnotation(which = "col",
                          tissue = anno_simple(tiannot,
                                                   col = legendCols$tissue,
                                                   border = T),
                          treatment = anno_simple(treatannot,
                                                col = legendCols$treatment,
                                                border = T),
                          replicate = anno_simple(repanoot,
                                             col = legendCols$replicate,
                                             border = T),
                          annotation_label = c("Source",
                                               "Condition",
                                               "Replicate")
)

# defining legends
legs = list(
  tissue = Legend(title = "Source",
                  at = legendCols$tissue %>% names(),
                  legend_gp = gpar(fill = legendCols$tissue %>% unname()),
                  border = "black",
                  labels = c("Leaf", "Cell culture")),
  treatment = Legend(title = "Condition",
                     at = legendCols$treatment %>% names(),
                     legend_gp = gpar(fill = legendCols$treatment %>% unname()),
                     border = "black",
                     labels = c("Regular light", "UV light",
                                "Lag phase", "Log phase", "Stationary phase")),
  replicate = Legend(title = "Replicate",
                     at = legendCols$replicate %>% names(),
                     legend_gp = gpar(fill = legendCols$replicate %>% unname()),
                     border = "black",
                     labels = c("A", "B", "C"))
)

# plotting heat map
ht = Heatmap(M[select,], col = colfunct,
             top_annotation = annot,
             border = T,
             show_column_names = F,
             show_row_names = F,
             heatmap_legend_param = list(
               title = expression(VST(Abundance)),
               border = T,
               direction = "horizontal"))

drawnht = grid.grabExpr(draw(ht,
                             heatmap_legend_side = "bottom",
                             annotation_legend_list = legs))

# heat map of genes related to the 
# setting up interesting genes
tps_bas_p450 = tibble(gene_id = c("PB.10364",
                                  "PB.19715",
                                  "PB.20962",
                                  "PB.19022",
                                  "PB.9352",
                                  "PB.9353",
                                  "PB.9356",
                                  "PB.11842",
                                  "PB.17050",
                                  "PB.17853"),
                      product = c("bAS (Beta-amyrin synthase)",
                                  "P450 (C-28 oxidase)",
                                  "P450 (C-16alpha oxidase)",
                                  "P450 (C-23 oxidase)",
                                  rep("Terpene synthase", 6))) %>% 
  mutate(rowName = structure(paste0(gene_id, "; ", product), names = gene_id))

# check who is included in M
# because PB.9356 got kicked out on 
# a filter step for not having at least
# 10 TPM summing all samples
tps_bas_p450 = tps_bas_p450 %>% 
  filter(tps_bas_p450$gene_id %in% rownames(M))

# quillaic acid biosynthesis
tps_bas_p450_ht = Heatmap(M[rownames(M) %in% tps_bas_p450$gene_id,], col = colfunct,
                          top_annotation = annot,
                          border = T,
                          show_column_names = F,
                          show_row_names = T,
                          row_labels = tps_bas_p450$rowName[match(rownames(M[rownames(M) %in% tps_bas_p450$gene_id,]),
                                                                  tps_bas_p450$gene_id)],
                          row_names_max_width = unit(10, "cm"),
                          heatmap_legend_param = list(
                            title = expression(VST(Abundance)),
                            border = T,
                            direction = "horizontal"))

drawnht_bas_p450_ht = grid.grabExpr(draw(tps_bas_p450_ht,
                                         heatmap_legend_side = "top",
                                         annotation_legend_side = "bottom",
                                         annotation_legend_list = legs))

ggsave(filename = "plots/tps_bas_p450_heatmap.png",
       plot = drawnht_bas_p450_ht,
       device = "png",
       units = "in",
       dpi = 300,
       height = 5,
       width = 6)

# comparing interesting genes between conditions
vsdBiosynth = assay(vsd) %>%
  as_tibble(rownames = "gene_id") %>% 
  filter(gene_id %in% tps_bas_p450$gene_id)

colnames(vsdBiosynth) = c("gene_id", paste0(colData$tissue, "_", colData$sample))

vsdBiosynthLong = vsdBiosynth %>%
  pivot_longer(cols = cell_culture_C1:leaf_F7,
               names_to = "source",
               names_pattern = "(.*)_.*$",
               values_to = "values")

compar = list(c("leaf", "cell_culture"))

jitCompare = vsdBiosynthLong %>% 
  ggplot(aes(x = factor(source, levels = c("leaf", "cell_culture")),
             y = values,
             color = source)) +
  geom_jitter(size=.5, position = position_jitter(width=.075),
              show.legend = F) +
  stat_compare_means(mapping = aes(label = ..p.signif..),
                     method = "t.test",
                     label.y.npc = "top",
                     label.x.npc = "center",
                     comparisons = compar) +
  facet_wrap(~ factor(gene_id,
                      levels = c("PB.17050",
                                 "PB.20962",
                                 "PB.19022",
                                 "PB.19715",
                                 "PB.9353",
                                 "PB.10364",
                                 "PB.9352",
                                 "PB.11842",
                                 "PB.17853"))) +
  ylim(c(5,20)) +
  ylab("VST(Abundance)") +
  xlab(label = NULL) +
  scale_x_discrete(labels = c("cell_culture" = "Cell culture",
                              "leaf" = "Leaf")) +
  scale_color_manual(values = c("leaf" = "#4E79A7",
                                "cell_culture" = "#59A14F")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
bas_tps_p450_panel = ggarrange(plotlist = list(drawnht_bas_p450_ht,
                                               jitCompare),
                               labels = "AUTO",
                               widths = c(2,1))

ggsave(filename = "plots/tps_bas_p450_panel.png",
       plot = bas_tps_p450_panel,
       device = "png",
       units = "in",
       dpi = 300,
       height = 5,
       width = 8.75)

# pca plot
# PCA
vsdPCA = assay(vsd) %>% 
  as_tibble(rownames = "locus_tag") %>% 
  transpose_tibble(locus_tag, id_col = "sample")

# adapted from:
# https://tbradley1013.github.io/2018/02/01/pca-in-a-tidy-verse-framework/
# creating a dataframe containing all pca info we need
dfpcaComplete = vsdPCA %>% 
  tidyr::nest(data=everything()) %>% 
  dplyr::mutate(pca = purrr::map(data, ~ prcomp(.x %>% dplyr::select(-sample), 
                                                center = T, scale = T)),
                pca_aug = purrr::map2(pca, data, ~ broom::augment(.x, data = .y)))

# computing the explained variance for
# each of the principal components
expVar = dfpcaComplete %>% 
  unnest(pca_aug) %>% 
  summarize_at(.vars = vars(contains(".fittedPC")), .funs = list(~var(.))) %>% 
  gather(key = pc, value = variance) %>% 
  mutate(var_exp = variance/sum(variance),
         cum_var_exp = cumsum(var_exp),
         pc = str_replace(pc, ".fitted", ""))

# preparing df to plot PCA
pcaAugDf = dfpcaComplete %>%
  unnest(pca_aug) %>%
  mutate(tissue = colData$tissue,
         treatment = colData$treatment,
           replicate = colData$replicate)

# finding top three loadings in PC1
#dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = T),1] %>% head(3)

# finding bottom three loadings in PC1
#dfpcaComplete$pca[[1]]$rotation[dfpcaComplete$pca[[1]]$rotation %>% .[,1] %>% order(decreasing = F),1] %>% head(3)

# plotting 2D scatter (alternative version to the 3D scatter)
twodscatterplot = pcaAugDf %>% 
  ggplot(aes(x=.fittedPC1, y=.fittedPC2, color = treatment, shape = tissue)) +
  geom_point(size = 2.5) +
  xlab(paste0("PC1", " (", round(expVar$var_exp[1] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2", " (", round(expVar$var_exp[2] * 100, digits = 2), "%)")) +
  scale_shape_discrete(name="Source",
                       labels = c("leaf" = "Leaf",
                                  "cell_culture" = "Cell culture")) +
  scale_color_manual(name="Condition",
                     values = c("regular_light" = "#FF9DA7",
                                "uv_light" = "#B07AA1",
                                "lag_phase" = "#4E79A7",
                                "log_phase" = "#76B7B2",
                                "stationary_phase" = "#59A14F"),
                     labels = c("regular_light" = "Regular light",
                                "uv_light" = "UV light",
                                "lag_phase" = "Lag phase",
                                "log_phase" = "Log phase",
                                "stationary_phase" = "Stationary phase")) +
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(order = 1))

# euclidean distance plot
colors = circlize::colorRamp2(breaks = c(250,0), colors = c("white", "#4E79A7"))
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste(vsd$tissue, vsd$treatment, vsd$replicate, sep=":") %>% 
  str_replace_all("_", " ") %>% 
  str_replace("leaf", "Leaf") %>% 
  str_replace("cell", "Cell") %>% 
  str_replace("uv", "UV") %>% 
  str_replace("regular", "Regular") %>% 
  str_replace("lag", "Lag") %>% 
  str_replace("log", "Log") %>% 
  str_replace("stationary", "Stationary")

distHt = Heatmap(sampleDistMatrix, col = colors,
                 top_annotation = annot,
                 border = T,
                 show_column_names = F,
                 show_row_names = T,
                 heatmap_legend_param = list(
                   title = "Euclidean distance",
                   border = T,
                   direction = "horizontal"))

drawnDistHt = grid.grabExpr(draw(distHt,
                                 heatmap_legend_side = "bottom"))

# drawing panel
exploratoryPanel = ggarrange(plotlist = list(ggarrange(plotlist = list(twodscatterplot,
                                                                       drawnDistHt),
                                                       nrow = 2,
                                                       ncol = 1,
                                                       labels = c("A", "B"),
                                                       heights = c(1, 1.5)),
                                             drawnht),
                             nrow = 1,
                             ncol = 2,
                             labels = c("", "C"))

ggsave(filename = "plots/exploratoryPanel.png",
       plot = exploratoryPanel,
       device = "png",
       units = "in",
       dpi = 300,
       height = 7,
       width = 10.5)

# writing volcano plots ####
ggsave(filename = "plots/volcano_leaf_and_uv_light_vs_leaf_and_regular_light.png",
       plot = res$leaf_and_uv_light_vs_leaf_and_regular_light$volcano +
         scale_x_continuous(breaks = seq(-80,80,10)) +
         ggtitle(label = NULL),
       device = "png",
       units = "in",
       dpi = 300,
       height = 5,
       width = 6)
       

