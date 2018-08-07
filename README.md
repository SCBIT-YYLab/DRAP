
# DRAP
DRAP is an R package for drug response analysis tailored for preclinical drug testing on patient-derived xenograft models. 

Patient-derived xenograft (PDX) models represent a valuable platform for preclinical drug testing and personalized cancer therapies. After summarizing a series of literatures carrying out drug response studies on PDXs, we classified the emerging PDX preclinical settings into four patterns: 1-A-N, T-1-N, T-A-1 and T-A-N, with the first letter representing the number of tumors, the second representing the number of arms for each tumor, the third representing the number of animals corresponding to one tumor line in each arm. Note that one means single and T/A/N means multiple. This package provides tools for drug response analysis tailored for drug testing on PDX models, involving data visualization, data analysis, and conclusion representation for the four types of PDX trial settings. The data analysis methods include three type of methods: 1) Assess potential differences in tumor volume between arms with statistical methods, a. ANOVA; b. Kruskal-Wallis test; c. Scheirer-Ray-Hare test; d. mixed-design ANOVA; e. linear mixed model (LMM); f. permutation test. 2) Quantify drug efficacy with tumor growth inhibition (TGI) rate; 3) Label drug response level of each animal or each tumor, a. NPDXE.Response, label drug response level with the method built in Novartis Institutes for BioMedical Research PDX encyclopedia; PPTP.Response: label drug response level with the method build in the Pediatric Preclinical Testing Program; c. RC.Response, label drug response level with the method based on relative change of tumor volume.


# Installation

Use devtools to install DRAP package directly from github.

if(!require(devtoos)){
  install.packages('devtools')
}

devtools::install_github('SCBIT-YYLab/DRAP')


# Reference
Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.

Murphy, B., et al. Evaluation of Alternative In Vivo Drug Screening Methodology: A Single Mouse Analysis. Cancer Res 2016;76(19):5798-5809.

Bertotti, A., et al. The genomic landscape of response to EGFR blockade in colorectal cancer. Nature 2015;526(7572):263-267.
