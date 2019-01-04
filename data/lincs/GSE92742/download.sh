#!/usr/bin/env bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5Fsig%5Finfo%2Etxt%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742%5FBroad%5FLINCS%5FLevel5%5FCOMPZ%2EMODZ%5Fn473647x12328%2Egctx%2Egz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_metrics.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_metrics.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz
echo "To decompress the h5 file in limited space setting, use: https://unix.stackexchange.com/a/341473/329794"
