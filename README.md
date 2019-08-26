
# DRAP
DRAP is an R package for drug response analysis and visualization tailored for preclinical drug testing on patient-derived xenograft models. 

Patient-derived xenograft (PDX) models represent a valuable platform for preclinical drug testing and personalized cancer therapies. After summarizing a series of literatures carrying out drug response studies on PDXs, we classified the emerging PDX preclinical settings into four patterns: 1-A-N, T-1-N, T-A-1 and T-A-N, with the first letter representing the number of tumors, the second representing the number of arms for each tumor, the third representing the number of animals corresponding to one tumor line in each arm. Note that one means single and T/A/N means multiple. This package provides tools for drug response analysis tailored for drug testing on PDX models, involving data visualization, data analysis, and conclusion representation for the four types of PDX trial settings. The data analysis methods include three type of methods: 1) Assess potential differences in tumor volume between arms with statistical methods, a. ANOVA; b. Kruskal-Wallis test; c. Scheirer-Ray-Hare test; d. mixed-design ANOVA; e. linear mixed model (LMM); f. permutation test. 2) Quantify drug efficacy with tumor growth inhibition (TGI) rate; 3) Label drug response level of each animal or each tumor, a. NPDXE.Response, label drug response level with the method built in Novartis Institutes for BioMedical Research PDX encyclopedia; PPTP.Response: label drug response level with the method build in the Pediatric Preclinical Testing Program; c. RC.Response, label drug response level with the method based on relative change of tumor volume.



# Installation

Use devtools to install DRAP package directly from github.

if(!require(devtools)) install.packages("devtools")

devtools::install_github('SCBIT-YYLab/DRAP')



# Reference
Li, Q., et al. DRAP: a toolbox for drug response analysis and visualization tailored for preclinical drug testing on patient-derived xenograft models. Journal of Translational Medicine, 2019, 17(1): 39.


