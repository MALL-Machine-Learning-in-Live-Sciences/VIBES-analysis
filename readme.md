# PAPER TITLE

This repository includes code based on the following publication:

*link paper

DOI: link DOI

## Abstract

## Prerequisites:

Several packages have been used throughout this project. So, before running the code make sure that the following packages have been installed:

```{r}
# CRAN packages
install.packages(c('openxlsx', 'dada2', 'dplyr', 'data.table', 'rentrez', 'taxonomizr', 'taxize', 'readxl', 'tidyverse', 'tibble', 'reshape2', 'mlr3', 'mlr3fselect', 'mlr3tuning', 'mlr3learners', 'mlr3filters', 'mlr3pipelines', 'mlr3misc', 'mlr3learners', 'genalg', 'pbapply', 'glmnet', 'pheatmap', 'cowplot', 'magrittr', 'ggplot2', 'viridis','ggpubr', 'forestplot'))

# GitHub packages:
devtools::install_github("david-barnett/microViz")

# Bioconductor packages
BiocManager::install(c('phyloseq','microbiome', 'ConsensusClusterPlus', 'MOFA2', 'ComplexHeatmap'))

```

In addition, the [enaBrowserTools](https://github.com/enasequence/enaBrowserTools) have been used to download data from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home),  [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to make quality control checks on raw sequence data and adapted script from [Valencia](https://github.com/ravel-lab/VALENCIA) to compute CSTs.

## Project workflow
The starting point for Train and Validation 1 cohorts was the clinical data and the Operational Taxonomic Units (OTUs) count tables. The remaining validation cohorts were obtained by processing the 16S rRNA sequences hosted in the European Nucleotide Archive (ENA). Therefore, there are two workflows in [00_preprocess_cohorts](https://github.com/DiegoFE94/BV_Microbiome/tree/main/00_preprocess_cohorts/code)

### 00.Preprocess_cohorts
w1- From [00_preprocess_metadata.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/00_preprocess_metadata.r)-[04_make_pseq.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/04_make_pseq.r)

1. [00_preprocess_metadata.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/00_preprocess_metadata.r) This script extracts the data from the ENA Browser, filters samples and saves the metadata for each validation cohort.

2. [01_download.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/01_download.r) This script downloads the fastq file(s) for each cohort, using the sample names extracted in the previous script as a guide. It also applies to each sample the fastqc program for quality control.




## Motivation
The diagnosis of BF is based on a few and non-precise (bacteria species + physical) parameters. In addition, only a few patients will be benefited by this diagnosis (extreme patients).

Since bacterias live in an ecosystems, we need to use all of them to stratify patients. And, since for medical practice, is unviable obtain whole microbiome profile for every patient, our aim is to i) stratify patients according all microbiome profile and ii) build a ML able to predict this subgroups from a microbiome signature.

## First analysis

Clustering patients according all (as far as posible) microbiome profiles. How many clusters have we obtained? Which are the differences between clusters? Are the clusters really different in terms of microbiome?

In this first approximation we can stablish a new classification scheme of BV patients and provide more biological knowledge.

Methodology:

Clustering (ConsensusClusterPlus)

Pheatmap of all cohorts with cluster labels (k = 4)

Meta-analysis of differential abundances with ALDEX2 or AMCOMBC for multiclass

Results

First approximation with kmeans and 24 species we have obtained 4 (?) clusters (pending of validation!). Normal (N), Dysbiosis (D), Normal-Incipient Dysbiosis (IDN) and Incipient Disbiosis  (IDD).

### Second analysis

Is there any relation between clusters and drug response?

In addition we can use MEFISTO to see any change of microbes species / genus through treatment. Is there any specific bacterial changed after treatment? Does our model includes this bacteria?

Methodology:

MEFISTO

## Third analysis

Build a ML model able to predict cluster labels. First, we will explore several FS techniques and algorithms (or only one?) to build a robust predictive model. In order to validate the best model performance, we will make a two validation steps: i) CV in train subset and ii) external validation. For external validations, we will first label samples with our cluster method, and then predict those labels with the ML model.

Methodology:

mlr2 / ml3 to train models

Will we need remove batch effects? (if necessary, select a methodology)