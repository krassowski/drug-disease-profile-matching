# Title     : R dependencies installation script
# Objective : Instal dependencies
# Created by: krassowski
# Created on: 20/12/18

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# signatures scoring
BiocManager::install("GSVA", version = "3.8")
BiocManager::install("limma", version = "3.8")

# PanCancerAtlas - stratification retrival
BiocManager::install("TCGAbiolinks", version = "3.8")

BiocManager::install("PharmacoGx", version = "3.8")
# BiocManager::install("piano", version = "3.8")
# BiocManager::install("fgsea", version = "3.8")

# gsva.R
install.packages("pbapply")
install.packages("pbmcapply")

# ggplot2 extensions
install.packages('ggbeeswarm')
install.packages('ggalluvial')

# why is that not a part of ggplot?
install.packages('latex2exp')
