# alorenzetti 20211108

# description ####
# this script will perform the
# DE analysis for isoquill transcriptome

# getting started ####
# setting thresholds
qthr = 0.01
lfcthr = 1

# creating a tx2gene object
# for now, transcripts isoforms will be ignored
# the criterion for inclusion is:
# PB.1.1 and PB.1.2 will both comprise the gene PB.1
transcriptome = readDNAStringSet(filepath = "cupcake/hq_transcripts.fasta.collapsed.rep.fa",
                                 format = "fasta")

tx2gene = tibble(TXNAME = names(transcriptome) %>% str_replace("\\|.*$", ""),
                 GENEID = names(transcriptome) %>% str_replace("\\.[0-9]{1,}\\|.*$", ""))

# finding files
libs = paste0(list.dirs(path = "salmon_quant/libs", recursive = F), "/",
              "quant.sf")

# importing quant data
quant = tximport(files = libs,
                 type = "salmon",
                 tx2gene = tx2gene,
                 ignoreAfterBar = T)

# defining colData
sample = str_replace(libs, ".*/(.*?)/quant.*", "\\1")
tissue = c(rep("cell_culture", "9"), rep("leaf", "6"))
treatment = c(rep("lag_phase", 3), rep("log_phase", 3), rep("stationary_phase", 3),
              rep("regular_light", 3), rep("uv_light", 3))
replicate = rep(LETTERS[1:3], 5)
group = paste0(tissue, "_and_", treatment)

colData = tibble(sample = sample,
                 tissue = tissue,
                 treatment = treatment,
                 replicate = replicate,
                 group = group)

# dds object from tximport
#
# the function below takes care
# of implementing the method called
# "original counts and offset"
# proposed by tximport paper
dds = DESeqDataSetFromTximport(txi = quant,
                               colData = colData,
                               design = ~ group)

# filtering low count entries
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

# running DESeq2
dds = DESeq(dds)

# setting function to plot volcano
volcanoPlot=function(parsedRes = parsedRes, title = title){
  volcano = parsedRes %>% 
    ggplot(aes(x = log2FoldChange,
               y = -log10(padj),
               color = status,
               label = label)) +
    geom_point(show.legend = F,
               alpha = 0.75,
               shape = 16,
               stroke = 0) +
    geom_label_repel(show.legend = F,
                     color = "black",
                     alpha = 0.75,
                     max.overlaps = 30) +
    xlab("Log<sub>2</sub> (Fold Change)") +
    ylab("-Log<sub>10</sub> (Adjusted <i>p</i>-value)") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown()) +
    scale_color_manual(values = c("Upregulated" = tab10$red,
                                  "Downregulated" = tab10$blue,
                                  "Not significant" = "darkgrey")) +
    ggtitle(title)
  
  return(volcano)
}

# setting function to generate results
genRes=function(obj = obj, var = var, numerator = numerator, denominator = denominator){
  title = paste(numerator, "vs.", denominator)
  
  results = results(object = obj,
                    contrast = c(var, numerator, denominator),
                    alpha = qthr,
                    lfcThreshold = lfcthr,
                    altHypothesis="greaterAbs") %>% 
    as_tibble(rownames = "gene")
  
  minpadj = results$padj[!is.na(results$padj) & results$padj != 0] %>% min()
  minpadj = minpadj * 0.01
  
  resultsParsed = results %>% 
    mutate(status = case_when(log2FoldChange >= lfcthr & padj < qthr ~ "Upregulated",
                              log2FoldChange <= -lfcthr & padj < qthr ~ "Downregulated",
                              TRUE ~ "Not significant"),
           padj = case_when(padj == 0 ~ minpadj,
                            TRUE ~ padj))
  
  topup = resultsParsed %>%
    filter(status == "Upregulated") %>% 
    arrange(desc(log2FoldChange)) %>% 
    head(10) %>% 
    pull(gene)
  
  topdown = resultsParsed %>%
    filter(status == "Downregulated") %>% 
    arrange(log2FoldChange) %>% 
    head(10) %>% 
    pull(gene)
  
  resultsParsed = resultsParsed %>% 
    mutate(label = case_when(gene %in% c(topup, topdown) ~ gene,
                             TRUE ~ NA_character_))
  
  sigGenes = resultsParsed %>% 
    filter(status == "Upregulated" | status == "Downregulated") %>% 
    pull(gene)
  
  p = volcanoPlot(parsedRes = resultsParsed, title = title)
  
  resultsList = list(results,
                     resultsParsed,
                     p,
                     sigGenes)
  
  names(resultsList) = c("original",
                         "parsed",
                         "volcano",
                         "sigGenes")
  
  return(resultsList)
}

# getting results
res = list()

comps = list(c("leaf_and_uv_light", "leaf_and_regular_light"),
             c("cell_culture_and_log_phase", "cell_culture_and_lag_phase"),
             c("cell_culture_and_stationary_phase", "cell_culture_and_log_phase"),
             c("cell_culture_and_log_phase", "leaf_and_regular_light"),
             c("cell_culture_and_stationary_phase", "leaf_and_regular_light"),
             c("cell_culture_and_lag_phase", "leaf_and_regular_light"))

for(i in 1:length(comps)){
  contrastName = paste0(comps[[i]][1], "_vs_", comps[[i]][2])
  res[[contrastName]] = genRes(obj = dds,
                               var = "group",
                               numerator = comps[[i]][1],
                               denominator = comps[[i]][2])
}

# writing quant dataset imported
# by the tximport package
# the counts are
# "Aggregated transcript counts +
# average transcript length offsets"
quant2write = quant$abundance %>%
  as_tibble(rownames = "gene_id") %>% 
  left_join(x = .,
            y = funcatSeq %>% 
              filter(codan_status == "Full Length") %>% 
              group_by(gene_id) %>% 
              summarise(across(.cols = everything(),
                               .fns = ~ c(.x)[1])),
            by = "gene_id") %>% 
  dplyr::select(c(gene_id:V15, interpro_id, interpro_description)) %>% 
  .[order(quant2write$gene_id %>%
            str_replace("PB.", "") %>%
            as.numeric()),]

colnames(quant2write) = c("gene_id",
                          paste0(colData(vsd)$group, "_",
                                 colData(vsd)$replicate),
                          "interpro_id",
                          "interpro_description")

write.xlsx(x = quant2write,
           file = "results/quant.xlsx",
           overwrite = T)

# exploring data ####
# getting count matrix and
# applying appropriate normalization
# vsd = vst(dds, blind=FALSE)
# 
# # distance matrix
# sampleDists = dist(t(assay(vsd)))
# sampleDistMatrix = as.matrix(sampleDists)
# rownames(sampleDistMatrix) = paste(vsd$tissue, vsd$treatment, vsd$replicate, sep=":")
# colnames(sampleDistMatrix) = NULL
# colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# 
# # PCA
# pcaData = plotPCA(vsd, intgroup=c("tissue", "treatment", "replicate"), returnData=T)
# percentVar = round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=tissue)) +
#   geom_point() +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed() +
#   theme_bw()
# 
# # count matrix heat map
# select = order(rowVars(assay(vsd)),
#                decreasing=TRUE)[1:50]
# df = as.data.frame(colData(vsd)[,c("tissue","treatment")])
# pheatmap(assay(vsd)[select,],
#          cluster_rows=T,
#          show_rownames=T,
#          cluster_cols=T,
#          annotation_col=df, labels_col = paste0(colData$sample, ":", colData$replicate))
# 
# # expression of tps and p450 genes
# # using vsd
# select = rownames(vsd) %in% (tps_p450$all$transcript_id %>% str_replace("\\.[0-9]{1,2}$", ""))
# df = as.data.frame(colData(vsd)[,c("tissue","treatment")])
# rowannot = left_join(x = tibble(gene_id = rownames(vsd)),
#                      y = funcat %>%
#                        dplyr::select(gene_id, tps_p450) %>% 
#                        group_by(gene_id) %>% 
#                        summarise(tps_p450 = tps_p450 %>% unique()) %>% 
#                        drop_na(),
#                      by = "gene_id") %>% 
#   as.data.frame()
# rownames(rowannot) = rowannot$gene_id
# rowannot = rowannot %>% dplyr::select(-gene_id)
# 
# pheatmap(assay(vsd)[select,],
#          cluster_rows=T,
#          show_rownames=T,
#          cluster_cols=T,
#          annotation_col= df,
#          labels_col = paste0(colData$sample, ":", colData$replicate),
#          annotation_row = rowannot %>% filter(select),
#          annotation_colors = list("tps_p450" = c("p450" = tab10$blue, "tps" = tab10$red)),
#          fontsize_row = 5)
# 
# # using raw
# ddscopy = dds
# colnames(ddscopy) = 1:15
# 
# select = rownames(ddscopy) %in% (tps_p450$all$transcript_id %>% str_replace("\\.[0-9]{1,2}$", ""))
# df = as.data.frame(colData(ddscopy)[,c("tissue","treatment")])
# rowannot = left_join(x = tibble(gene_id = rownames(ddscopy)),
#                      y = funcat %>%
#                        dplyr::select(gene_id, tps_p450) %>% 
#                        group_by(gene_id) %>% 
#                        summarise(tps_p450 = tps_p450 %>% unique()) %>% 
#                        drop_na(),
#                      by = "gene_id") %>% 
#   as.data.frame()
# rownames(rowannot) = rowannot$gene_id
# rowannot = rowannot %>% dplyr::select(-gene_id)
# 
# pheatmap(log10(assay(ddscopy)[select,] + 1),
#          cluster_rows=T,
#          show_rownames=T,
#          cluster_cols=T,
#          annotation_col= df,
#          labels_col = paste0(colData$sample, ":", colData$replicate),
#          annotation_row = rowannot %>% filter(select),
#          annotation_colors = list("tps_p450" = c("p450" = tab10$blue, "tps" = tab10$red)),
#          fontsize_row = 5)
# 
# # p450 of interest
# p450c28c16ac23 = c("PB.19715", "PB.20962", "PB.20963", "PB.19022", "PB.19021")
# select = rownames(vsd) %in% p450c28c16ac23
# df = as.data.frame(colData(vsd)[,c("tissue","treatment")])
# rowannot = left_join(x = tibble(gene_id = rownames(vsd)),
#                      y = funcat %>%
#                        dplyr::select(gene_id, tps_p450) %>% 
#                        group_by(gene_id) %>% 
#                        summarise(tps_p450 = tps_p450 %>% unique()) %>% 
#                        drop_na(),
#                      by = "gene_id") %>% 
#   as.data.frame()
# rownames(rowannot) = rowannot$gene_id
# rowannot = rowannot %>% dplyr::select(-gene_id)
# 
# pheatmap(assay(vsd)[select,],
#          cluster_rows=T,
#          show_rownames=T,
#          cluster_cols=T,
#          annotation_col= df,
#          labels_col = paste0(colData$sample, ":", colData$replicate),
#          annotation_row = rowannot %>% filter(select),
#          annotation_colors = list("tps_p450" = c("p450" = tab10$blue, "tps" = tab10$red)),
#          fontsize_row = 5)
