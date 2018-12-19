#!/usr/bin/env bash

cancers="BRCA SKCM SARC GBMLGG GBM BLCA COAD ESCA LAML THCA TGCT STAD PCPG PAAD LUSC HNSC"

for cohort in ${cancers}; do
    wget "http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/$cohort/20160128/gdac.broadinstitute.org_$cohort.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0.tar.gz"
done
