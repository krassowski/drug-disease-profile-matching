# Drug-disease profile matching

Multi-omics disease sub-type specific drug repositioning aided with expression signatures from ConnectivityMap.

An MRes research project at Imperial College London.

### Overview

![](/images/profiles_and_expression.png?raw=true)

![](/images/scoring_functions.png?raw=true)

![](/images/stratifications.png?raw=true)


### Abstract

Drug discovery is a costly and difficult process. As a consequence, attempts to guide the selection of
drug candidates using machine learning and biomedical data has increased in recent years. One of
the approaches measures expression of candidate substances in cell lines to determine perturbation
profiles, which are used to generate compound similarity maps (e.g. the Connectivity Map). It was
proposed that the perturbation profiles may be also used to find drug candidates by matching the
profiles against differential expression patterns of diseases.

Previous studies demonstrated merits of multi-omics disease stratification, evaluating predictive ability
of novel clusters for cancer patients survival or analyzing the functional enrichment in the clusters.

In this work I apply perturbagen-disease profile matching to disease sub-types selected by multiple
multi-omics stratification methods, in order to prioritize new drug repositioning candidates. To
determine optimal method of matching the perturbation profiles to diseases, I evaluate performance of
16 methods (scoring functions), both known (six), and novel. I hypothesize that gene-set enrichment
(GSE) based methods may provide overall benefit by incorporation of additional biological information
and availability of stringent significance estimates.

I confirm previous findings demonstrating the ability of profile matching approaches to recover known
breast cancer drugs and highlight novel, previously unreported drug candidates.

I observe limited benefit of stratification for the drug recovery performance, with promising results
demonstrated by XSum and mROAST scoring functions. I show that the GSE-based approaches
require large numbers of samples, high-performance computing (HPC) facilities and may not increase
the chances of drug recovery in certain circumstances.

I examine whether certain scoring functions may be used to recognize whether breast cancer
stratifications are based on meaningful molecular clustering, using only their drug indications-contraindications classification performance. However, I find only limited support for this hypothesis.
Finally, I observe that within this framework, it may not be possible to answer the question if multiomics based stratifications perform better when used to support classification of a drug indications.

### Results

More details - including results - coming soon.


### Acknowledgements

The cells, RNA, DNA and histone pictograms are derivative works based on graphics from Reactome Icon Library (Creative Commons Attribution 4.0 International License)
