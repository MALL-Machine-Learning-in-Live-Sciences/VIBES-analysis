# The source code of the analyses conducted for the development of VIBES - VagInal Bacterial subtyping using machine learning for Enhanced classification of bacterial vaginosiS. 

Package documentation an use-case is available [here](https://mall-machine-learning-in-live-sciences.github.io/VIBES-docs/).

R package is available [here](https://mall-machine-learning-in-live-sciences.github.io/VIBES/).


## Citation

Article is open access [here](https://www.csbj.org/article/S2001-0370(23)00465-8/fulltext).

```
@article{vibes2023,
title = {VIBES: a consensus subtyping of the vaginal microbiota reveals novel classification criteria},
author = {D. Fernández-Edreira and J. Liñares-Blanco and P. V.-del-Río and C. Fernandez-Lozano},
editor = {Elsevier},
url = {https://www.csbj.org/article/S2001-0370(23)00465-8/fulltext},
doi = {https://doi.org/10.1016/j.csbj.2023.11.050},
issn = {2001-0370},
year = {2023},
date = {2023-11-30},
urldate = {2023-11-30},
journal = {Computational and Structural Biotechnology Journal},
abstract = {This study aimed to develop a robust classification scheme for stratifying patients based on vaginal microbiome. By employing consensus clustering analysis, we identified four distinct clusters using a cohort that includes individuals diagnosed with Bacterial Vaginosis (BV) as well as control participants, each characterized by unique patterns of microbiome species abundances. Notably, the consistent distribution of these clusters was observed across multiple external cohorts, such as SRA022855, SRA051298, PRJNA208535, PRJNA797778, and PRJNA302078 obtained from public repositories, demonstrating the generalizability of our findings. We further trained an elastic net model to predict these clusters, and its performance was evaluated in various external cohorts. Moreover, we developed VIBES, a user-friendly R package that encapsulates the model for convenient implementation and enables easy predictions on new data. Remarkably, we explored the applicability of this new classification scheme in providing valuable insights into disease progression, treatment response, and potential clinical outcomes in BV patients. Specifically, we demonstrated that the combined output of VIBES and VALENCIA scores could effectively predict the response to metronidazole antibiotic treatment in BV patients. Therefore, this study's outcomes contribute to our understanding of BV heterogeneity and lay the groundwork for personalized approaches to BV management and treatment selection.},
note = {Q1, 60/285 BIO-MB, 6 IF},
keywords = {},
pubstate = {published},
tppubtype = {article}
}
```

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

## Cohorts

[SRA022855](https://www.ebi.ac.uk/ena/browser/view/SRA022855) from [Ravel et al.](https://www.pnas.org/doi/10.1073/pnas.1002611107): Vaginal bacterial communities of 396 North American women (Baltimore and Atlanta), between the ages of 12 and 45 (mean age: 30). Women included in the study self-identified as African American (n=104), White (n=98), Hispanic (n=97), and Asian (n=97).


[SRA051298](https://www.ebi.ac.uk/ena/browser/view/SRA051298) from [Srinivasan et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037818):  Vaginal bacterial communities of 220 North American women (Seattle), between the ages of 18 and 57 (mean age: 29). Women included in the study self-identified as African American (n=75), White (n=97), Asian (n=15), and Others (n=33).

[PRJNA208535](https://www.ebi.ac.uk/ena/browser/view/PRJNA208535) from [Ravel et al.](https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-1-29): Vaginal bacterial communities of 25 North American women (Birmingham) over 10 weeks (1657 16S raw samples on repo.), between the ages of 19 and 45 (mean age: 27). Women included in the study self-identified as African American (n=20), and White (n=5).

[PRJNA797778](https://www.ebi.ac.uk/ena/browser/view/PRJNA797778) from [France et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02635-9): Vaginal bacterial communities of 39 North American women (Baltimore) over 10 weeks (220 16S raw samples on repo.), between the ages of 19 and 45. Women included in the study self-identified as African American (n=24), White (n=10), Hispanic (n=4), and Asian (n=1).

[PRJNA302078](https://www.ebi.ac.uk/ena/browser/view/PRJNA302078) from [Xiao et al.](https://pubmed.ncbi.nlm.nih.gov/27253522/): Vaginal bacterial communities of 65 Chinese women (Beijing) over 3 time points: pretreatment, one week after treatment and one month after treatment with metronidazole (201 16S raw samples on repo.), between the ages of 18 and 53.

## Project workflow

The 16S samples of the cohorts Validation 2, 3 and 4 have been downloaded from the ENA and further processed. The corresponding scripts are [00_preprocess_metadata.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/00_preprocess_metadata.r) to [04_make_pseq.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/04_make_pseq.r). 

As explained in the article, the Discovery and Validation 1 cohorts have been processed differently, due to the impossibility of reliably downloading them from the ENA. The corresponding scripts are from [05_ravel_taxa_acquisition.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/05_ravel_taxa_acquisition.r) to [06_get_phyloseqs.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/06_get_phyloseqs.r). 

So, there are two workflows in [00_preprocess_cohorts](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/tree/main/00_preprocesMALL-Machine-Learning-in-Live-Sciencess_cohorts/code) converging in the script [07_get_intersect.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/07_get_intersect.r).

### 00.Preprocess cohorts

[00_preprocess_metadata.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/00_preprocess_metadata.r) This script extracts the data from the ENA Browser, filters samples and saves the metadata for each validation cohort. This script saves a ```.rds``` with metadata by cohort.

[01_download.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/01_download.r) This script downloads the ```FASTQ``` file(s) for each cohort, using the sample names extracted in the previous script as a guide. It also applies to each sample the FastQC program for quality control, so we will have one ```html``` per sample.

[02_check_primers.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/02_check_primers.r) This script checks the ```FASTQ``` files for primers (provided in the original papers of each cohort).

[03_pipeline_16S.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/03_pipeline_16S.r) - [03_pipeline_paired_16S.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/03_pipeline_paired_16S.r) These scripts process the ```FASTQ``` (single or paired-end) files following the pipeline established in the DADA2 package. Both scripts return the ```ASV table``` and the ```taxonomic table``` for each cohort. Both in ```.rds``` format.

[04_make_pseq.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/04_make_pseq.r) In this script, the ```phyloseq``` object is built from the ```ASV table```, the ```taxonomic table``` and the ```metadata```.

[05_ravel_taxa_acquisition.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/05_ravel_taxa_acquisition.r) - [05_sriniv_taxa_acquisition.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/05_sriniv_taxa_acquisition.r) In both scripts, from the ```OTU tables```, the different taxonomic levels to which each species belongs are checked and the correct ```taxonomic tables``` is returned. Both in ```.rds``` format.
 
[06_get_phyloseqs.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/06_get_phyloseqs.r) In this script, the ```phyloseq``` object is built from the ```OTU table```, the ```taxonomy table``` and the ```metadata```.

[07_get_intersect.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/00_preprocess_cohorts/code/07_get_intersect.r) In the following script we homogenise all phyloseq objects and check which species are shared by the 5 cohorts. The output is the ```phyloseqs``` with only the species shared by all cohorts.

### 01.Get Valencias CSTs

[00_prepare_data.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/01_get_valencias/code/00_prepare_data.r) This script transforms the phyloseq objects of each cohort to compute the VALENCIA CSTs. It returns a ```.csv``` per cohort adapted for the VALENCIA computation script. 

[01_calculate_valencias.py](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/01_get_valencias/code/01_calculate_valencias.py) Computes the VALENCIA CSTs. Returning a ```.csv``` per cohort with the scores and the membership of each CST.

### 02.Clusterization

[00_cluster_consensus_plus.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/02_cluster/code/00_cluster_consensus_plus.r) This script runs the Cluster Consensus Plus analysis to find the best parameter settings for stratifying patients. The output is two ```.pdf``` files containing the results of the CCP.

[01_clustering.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/02_cluster/code/01_clustering.r) In this script, clustering is performed on the Discovery cohort and prediction is run on the remaining cohorts. CLR transformation is also performed on all cohorts. The input are the ```phyloseqs``` of all the cohorts with the 22 species. The output is the same ```phyloseqs``` but with the membership of each cluster added to the metadata.

### 03.Machine Learning

[00_prepare_data.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/03_machine_learning/code/00_prepare_data.r) This script prepare all cohorts to run machine learning analyses. From the phyloseq of each cohort it generates a ```.rds``` dataframe with the 22 species and the labels of the cluster to which each sample belongs.

[01_run_ml](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/tree/main/03_machine_learning/code/01_run_ml) This folder holds all the scripts needed to run each dataframe + algorithm combination in parallel. Both direct counts and CLR transformed data were used. Returns a ```.rds``` dataframe with the nested cross validation benchmark training.

[02_train_cv_results.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/03_machine_learning/code/02_train_cv_results.r) This script details the performance measures of the training. The input is the various benchmarks from the previous section and returns a ```.rds``` dataframe with the summary of the training according to the performance measures you specify.

[03_benchmark_external_validation.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/03_machine_learning/code/03_benchmark_external_validation.r) This script runs the external validation on the validation cohorts. It returns a summary ```.rds``` dataframe of performance measures across all cohorts.

[04_prune_best_glmnet.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/03_machine_learning/code/04_prune_best_glmnet.r) This script refines the best model found in this analysis. This is done in order not to have dependencies when using that model. As it is a ```glmnet``` model we extract the betas of each species. The output is a list in format  ```.rds``` with betas, nclasses and nlambdas.

[05_best_glmnet_external_validation.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/03_machine_learning/code/05_best_glmnet_external_validation.r) In this script we perform again a external validation in all cohorts but only with the best model. It returns a summary ```.rds``` dataframe of performance measures across all cohorts.

### 04.Treatment 

[00_preprocess_PRJNA302078.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/04_treatment/code/00_preprocess_PRJNA302078.r) In this script, the Validation 4 cohort is prepared with all species for MEFISTO analysis. It generates two ```.rds``` dataframes.

[01_run_mefisto.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/04_treatment/code/01_run_mefisto.r) In this script a MEFISTO analysis is executed. It returns a ```MOFA``` object with the model resulting from the analysis.

[02_plot_mefisto.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/04_treatment/code/02_plot_mefisto.r) This script is used to graphically display the results of the MEFISTO analysis.

[03_prepare_D0.r](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/blob/main/04_treatment/code/03_prepare_D0.r) In this script, three datasets are prepared from the D0 samples of Validation Cohort 4 (only information from our clusters, only information from the VALENCIA CSTs, and both). It generates three ```.rds``` dataframes. Subsequently, an identical ML methodology like [01_run_ml](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/tree/main/03_machine_learning/code/01_run_ml) is applied.

### Figures

[Figures](https://github.com/MALL-Machine-Learning-in-Live-Sciences/BV_Microbiome/tree/main/figures/code) folder contains scripts to reproduce each figures of the paper. Note that for aesthetics reasons, the figures were edited using illustrator.