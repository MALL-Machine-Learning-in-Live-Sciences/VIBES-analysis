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

### 00.Preprocess cohorts
w1 - From [00_preprocess_metadata.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/00_preprocess_metadata.r) - [04_make_pseq.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/04_make_pseq.r) In this first workflow, Validation cohorts 2, 3 and 4 are processed. 

[00_preprocess_metadata.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/00_preprocess_metadata.r) This script extracts the data from the ENA Browser, filters samples and saves the metadata for each validation cohort.

[01_download.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/01_download.r) This script downloads the FASTQ file(s) for each cohort, using the sample names extracted in the previous script as a guide. It also applies to each sample the FastQC program for quality control.

[02_check_primers.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/02_check_primers.r) This script checks the FASTQ files for primers (provided in the original papers of each cohort).

[03_pipeline_16S.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/03_pipeline_16S.r) - [03_pipeline_paired_16S.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/03_pipeline_paired_16S.r) These scripts process the FASTQ (single or paired-end) files following the pipeline established in the DADA2 package. Both scripts return the ASV table and the taxonomic assignment table for each ASV.

[04_make_pseq.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/04_make_pseq.r) In this script, the phyloseq object is built from the ASV table, the taxonomy table and the metadata.

w2 - From [05_ravel_taxa_acquisition.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/05_ravel_taxa_acquisition.r) - [06_get_phyloseqs.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/06_get_phyloseqs.r) In this second workflow, the Discovery and Validation 1 cohorts are processed.

[05_ravel_taxa_acquisition.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/05_ravel_taxa_acquisition.r)] - [05_sriniv_taxa_acquisition.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/05_sriniv_taxa_acquisition.r) In both scripts, from the OTU tables, the different taxonomic levels to which each species belongs are checked and the correct taxonomic table is returned.
 
[06_get_phyloseqs.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/06_get_phyloseqs.r) In this script, the phyloseq object is built from the OTU table, the taxonomy table and the metadata.

[07_get_intersect.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/00_preprocess_cohorts/code/07_get_intersect.r) In the following script we homogenise all phyloseq objects and check which species are shared by the 5 cohorts. The output is the phyloseqs with only the species shared by all cohorts.

### 01.Get Valencias CSTs

[01_prepare_data.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/01_get_valencias/code/01_prepare_data.r) This script transforms the phyloseq objects of each cohort to compute the VALENCIA CSTs. It returns a CSV per cohort adapted for the VALENCIA computation script. 

[02_calculate_valencias.py](https://github.com/DiegoFE94/BV_Microbiome/blob/main/01_get_valencias/code/02_calculate_valencias.py) Computes the VALENCIA CSTs. Returning the scores and the membership of each CST.

### 02.Clusterization

### 03.Machine Learning

### 04.Treatment

### extdata

### Figures
[figures](https://github.com/DiegoFE94/BV_Microbiome/tree/main/figures/code) folder contains scripts to reproduce each figures of the paper. Note that for aesthetics reasons, the figures were edited using illustrator.