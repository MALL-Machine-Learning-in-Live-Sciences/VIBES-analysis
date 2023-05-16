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

In addition, the [enaBrowserTools](https://github.com/enasequence/enaBrowserTools) have been used to download data from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home),  [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to make quality control checks on raw sequence data and adapted script from [VALENCIA](https://github.com/ravel-lab/VALENCIA) to compute CSTs.

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

[00_prepare_data.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/01_get_valencias/code/00_prepare_data.r) This script transforms the phyloseq objects of each cohort to compute the VALENCIA CSTs. It returns a CSV per cohort adapted for the VALENCIA computation script. 

[01_calculate_valencias.py](https://github.com/DiegoFE94/BV_Microbiome/blob/main/01_get_valencias/code/01_calculate_valencias.py) Computes the VALENCIA CSTs. Returning the scores and the membership of each CST.

### 02.Clusterization

[00_cluster_consensus_plus.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/02_cluster/code/00_cluster_consensus_plus.r) This script runs the Cluster Consensus Plus analysis to find the best parameter settings for stratifying patients.

[01_clustering.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/02_cluster/code/01_clustering.r) In this script, clustering is performed on the Discovery cohort and prediction is run on the remaining cohorts. CLR transformation is also performed on all cohorts. The input is the phyloseqs of all the cohorts with the 22 species. The output is the same phyloseqs but with the membership of each cluster added to the metadata.

### 03.Machine Learning

[00_prepare_data.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/03_machine_learning/code/00_prepare_data.r) This script prepare all cohorts to run machine learning analyses. From the phyloseq of each cohort it generates a dataframe with the 22 species and the labels of the cluster to which each sample belongs.

[01_run_ml](https://github.com/DiegoFE94/BV_Microbiome/tree/main/03_machine_learning/code/01_run_ml) This folder holds all the scripts needed to run each dataframe + algorithm combination in parallel. Both direct counts and CLR transformed data were used. Returns the nested cross validation benchmark training.

[02_train_cv_results.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/03_machine_learning/code/02_train_cv_results.r) This script details the performance measures of the training. The input is the various benchmarks from the previous section and returns an overall summary of the training according to the performance measures you specify.

[03_benchmark_external_validation.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/03_machine_learning/code/03_benchmark_external_validation.r) This script runs the external validation on the validation cohorts. It returns a summary of performance measures across all cohorts.

[04_prune_best_glmnet.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/03_machine_learning/code/04_prune_best_glmnet.r) This script refines the best model found in this analysis. This is done in order not to have dependencies when using that model. As it is a ```glmnet``` model we extract the betas of each species.

[05_best_glmnet_external_validation.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/03_machine_learning/code/05_best_glmnet_external_validation.r) In this script we perform again a external validation in all cohorts but only with the best model.

### 04.Treatment 

[00_preprocess_PRJNA302078.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/04_treatment/code/00_preprocess_PRJNA302078.r) In this script, the Validation 4 cohort is prepared with all species for MEFISTO analysis.

[01_run_mefisto.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/04_treatment/code/01_run_mefisto.r) In this script a MEFISTO analysis is executed. It returns a MOFA object with the model resulting from the analysis.

[02_plot_mefisto.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/04_treatment/code/02_plot_mefisto.r) This script is used to graphically display the results of the MEFISTO analysis.

[03_prepare_D0.r](https://github.com/DiegoFE94/BV_Microbiome/blob/main/04_treatment/code/03_prepare_D0.r) In this script, three datasets are prepared from the D0 samples of Validation Cohort 4 (only information from our clusters, only information from the VALENCIA CSTs, and both). Subsequently, an identical ML methodology like [01_run_ml](https://github.com/DiegoFE94/BV_Microbiome/tree/main/03_machine_learning/code/01_run_ml) is applied.

### Figures
[figures](https://github.com/DiegoFE94/BV_Microbiome/tree/main/figures/code) folder contains scripts to reproduce each figures of the paper. Note that for aesthetics reasons, the figures were edited using illustrator.

## Review

## Citation

## Contact

If you have any questions, comments, or suggestions, please feel free to
contact us at:

- Diego Fernández Edreira
  - Email: <diego.fedreira@udc.es>
  - Twitter: [@diego_edreira](https://twitter.com/diego_edreira)
  - GitHub: [DiegoFE94](https://github.com/DiegoFE94/)
- Jose Liñares Blanco
  - Email: <j.linares@udc.es>
  - Twitter: [@8JoseLinares](https://twitter.com/8JoseLinares)
  - GitHub: [jlinaresb](https://github.com/jlinaresb)
- Carlos Fernández Lozano
  - Email: <carlos.fernandez@udc.es>
  - Twitter: [@cafernandezlo](https://twitter.com/cafernandezlo)
  - GitHub: [cafernandezlo](https://github.com/cafernandezlo)