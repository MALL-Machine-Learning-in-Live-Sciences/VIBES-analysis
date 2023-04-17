# RUN MEFISTO
# =======
require(MOFA2)
require(tidyverse)
require(pheatmap)
require(reshape2)
library(cowplot)
library(magrittr)

# Arguments
setwd("~/git/vaginosis-jlb/")
cohort <- "PRJNA302078"   # PRJNA797778 PRJNA302078
inputfile <- file.path("MEFISTO/data/", paste0(cohort, "_microbiome.rds"))
outfile <- paste0("MEFISTO/res/", cohort, ".rds")
time_variable <- "time"
nfactors <- 3

# Load preprocessed data
microbiome = readRDS(inputfile)

# Create MOFA object (only with microbiome data)
sm <- create_mofa_from_df(df = microbiome)

# Specify time covariate
sm <- set_covariates(
  sm,
  covariates = time_variable)

# MOFA and MEFISTO PARAMETERS
data_opts <- get_default_data_options(sm)
data_opts$center_groups <- FALSE

model_opts <- get_default_model_options(sm)
model_opts$num_factors <- nfactors

mefisto_opts <- get_default_mefisto_options(sm)
mefisto_opts$n_grid <- 10
mefisto_opts$start_opt <- 50
mefisto_opts$opt_freq <- 50
# mefisto_opts$new_values <- matrix(
#   min(microbiome$time):max(microbiome$time), nrow =1)    # add as argument!?

train_opts <- get_default_training_options(sm)
train_opts$seed <- 2020

# Prepare MOFA object
sm <- prepare_mofa(
  sm, 
  model_options = model_opts,
  mefisto_options = mefisto_opts,
  training_options = train_opts,
  data_options = data_opts
  )

# Run MOFA
sm <- run_mofa(
  sm, 
  use_basilisk = TRUE,
  outfile = outfile)
