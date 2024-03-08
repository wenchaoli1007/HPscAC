# scImmunoAging
### scImmunoAging: the R implementation of the cell-type-specific transcriptome aging clocks for human PBMC. 

In this study, we established a robust cell-type-specific aging clock, covering monocytes, CD4+ T, CD8+ T, NK and B cells, based on single-cell transcriptomic 
profiles from 1081 PBMC samples from European healthy adults. Our research sheds light on understanding biological age alterations in response to vaccinations 
and diseases, revealing the most relevant cell type and subset of genes that play dual roles in both aging and immune responses to various stimuli.

We describe the HPscAC in the following paper: Cell-Type-Specific Aging Clocks Unveil Inter-Individual Heterogeneity in Immune Aging and Rejuvenation during Infection and Vaccination

## Introduction to HPscAC
We developed cell-type-specific transcriptome aging clocks for human PBMCs using scRNA-seq datasets from five studies, encompassing 1081 healthy individuals of 
European ancestry aged 18 to 97 years. Focusing on the five most prevalent cell types - CD4+ T cells, CD8+ T cells, monocytes, NK cells, and B cells- we build 
independent aging clocks for each using machine learning (LASSO and random forest) and deep learning (PointNet) methods, assessing performance through various 
metrics.

![Workflow of HPscAC](https://github.com/wenchaoli1007/HPscAC/blob/main/data/workflow.png)

## Please install the following packages
glmnet, dplyr, ggplot2, purrr

## Installation of HPscAC
install.packages("devtools")

devtools::install_github("wenchaoli1007/HPscAC")

## References
Cell-Type-Specific Aging Clocks Unveil Inter-Individual Heterogeneity in Immune Aging and Rejuvenation during Infection and Vaccination


