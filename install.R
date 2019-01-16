# Title     : R dependencies installation script
# Objective : Instal dependencies
# Created by: krassowski
# Created on: 20/12/18

source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA", version = "3.8")
BiocManager::install("limma", version = "3.8")
install.packages("pbapply")
install.packages("pbmcapply")
