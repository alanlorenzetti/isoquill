# alorenzetti 20211108

# description ####
# this script will load
# required libs and set up a few variables
# and will create directories

# getting started ####
library(pacman)

# it also requires two packages
# from third-party repo
# orthologr requires ncbi-blast+
# installed in the system

# install metablastr from GitHub
# devtools::install_github("HajkD/metablastr")

# install orthologr from GitHub
# devtools::install_github("HajkD/orthologr")

packs = c("tidyverse",
          "DESeq2",
          "ggthemes",
          "readtext",
          "Biostrings",
          "ggpubr",
          "ggtext",
          "ggrepel",
          "ggExtra",
          "viridis",
          "ComplexHeatmap",
          "tximport",
          "openxlsx",
          "RColorBrewer",
          "pheatmap",
          "taxonomizr",
          "devtools",
          "GenomicRanges",
          "GenomicFeatures",
          "Rsamtools",
          "rtracklayer",
          "orthologr",
          "tblhelpr",
          "openxlsx")

p_load(char = packs)

# loading color scheme
tab10 = list()
tab10$blue = ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1]
tab10$red = ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[3]

# setting ggplot2 blank theme
theme_set(theme_bw())

# creating results and plots
# directories
for(i in c("results", "plots")){
  if(!dir.exists(i)){dir.create(i)}
}
