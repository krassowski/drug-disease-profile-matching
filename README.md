# Drug-disease profile matching

Multi-omics disease sub-type specific drug repositioning aided with expression signatures from ConnectivityMap.

An MRes research project at Imperial College London.

### Overview

![](/images/profiles_and_expression.png?raw=true)

![](/images/scoring_functions.png?raw=true)

![](/images/stratifications.png?raw=true)


### Abstract

#### Background
Attempts to guide the selection of drug candidates with machine learning has increased in recent years.
One of the popular approaches is so called "guilt-by-association" (GBA), using drug-drug and disease-disease similarity
to guide the drug candidates selection.
Comparison of genes expression profiles (perturbation profiles) after treatment with the candidate substances (perturbagenes)
is a well established approach to generate compound similarity maps for GBA methods. Large-scale projects, such as The Connectivity Map, offer ways of systematic perturbagene screening.

It was proposed that the perturbation profiles may be also used to find drug candidates by matching the profiles against
differential expression patterns of diseases (being a simpler alternative to advanced machine-learning methods).

Previous studies demonstrated merits of multi-omics disease stratification, evaluating predictive ability
of novel clusters for cancer patients survival or analyzing the functional enrichment in the clusters.

#### Introduction

In this work perturbagen-disease profile matching is applied to diseases and disease sub-types selected by multiple
multi-omics stratification methods, in order to prioritize new drug repositioning candidates.

Multiple perturbation profile-disease expression matching methods (scoring functions) are evaluated,
and then applied to cancer cohorts having enough data. 

Gene-set enrichment (GSE) based methods are hypothesized to provide overall benefit by incorporation of additional biological
information and availability of stringent significance estimates.

Finally, a hypothesis that scoring functions may be used to recognize stratifications based on meaningful molecular clustering,
using only drug indications-contraindications classification performance is proposed.

#### Materials and methods

Cancer data from The Cancer Genome Atlas are used, with an extensive case study on breast carcinoma (BRCA) cohort,
limited validation with prostate (PRAD) and skin (STAD) adenocarcinomas, and a pan-cancer analysis.


Indications-contraindications classification performance is used for scoring function evaluation.
Evaluation was performed on 16 scoring functions of which six were proposed in previously works.
Six scoring functions are chosen and applied in further analyses.

Performance of the scoring functions is compared using four previously published
stratifications of breast cancer (including three based on multi-omics data): PAM50, iCluster, PARADIGM and Pan-Gyn.


#### Results
Ability of profile matching approaches to recover known drugs (as previously reported) is confirmed.

A few previously unreported breast-cancer drug candidates are highlighted.

The advantages and disadvantages of proposed indications-contraidications classification use (with many cancer drugs being known carcinogens) is discussed.

GSE-based methods require large numbers of samples, high-performance computing facilities and may not increase the chances of drug recovery in certain circumstances.

While the results obtained with meaningful stratifications do not always perform better than random permutations,
limited benefit of stratification is observed for the drug recovery performance, with promising results from XSum and mROAST scoring functions.

Despite no definite evidence for superiority of multi-omics stratifications use for classification of drug indications/contraindications,
two multi-omics stratifications are highlighted as performing better than others: PARADIGM and Pan-Gyn.


### Acknowledgements

The cells, RNA, DNA and histone pictograms are derivative works based on graphics from [Reactome Icon Library](https://reactome.org/icon-lib) (licensed under [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by-sa/4.0/)).
