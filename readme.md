## Background

The diagnosis of BF is based on a few and non-precise (bacteria species + physical) parameters. In addition, only a few patients will be benefited by this diagnosis (extreme patients).

## Motivation

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