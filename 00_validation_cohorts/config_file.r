# DADA2 config file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd(dir = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome")
# Declare global paths
experiment_name <- "PRJNA797778"
input_dir_path <- paste0(experiment_name,"/")
out_dir_path  <- "out_data"
out_path <- paste(out_dir_path, "/", experiment_name, sep = "")
if (dir.exists(out_path) == FALSE) {
  dir.create(out_path)
  message(paste("Creating", experiment_name, "directory!"))
}
# Run accession pattern from ENA BROWSER 
pattern = "SRR"
# Declare DADA2 Parameters (Care if use single or pair pipeline)
## Filter parameters
dada2_trunclen <- c(249, 249)
dada2_maxn <- 0
dada2_maxee <- c(2, 3)
dada2_truncq <- 2
dada2_trimleft <- 0 # Default: 0. Ion Torrent sequencing: 15

## Learn Errors parameters
dada2_nbases <-  1e+08# Default: 1e+08

## DADA2 core algorithm parameters
dada2_homopolymer_gap_penalty <- -8 # Default: -8.IT and 454: -1
dada2_band_size <- 16 # Default: 16.IT and 454: 32
dada2_pool <- FALSE # Default: FALSE. If set TRUE -> quimera method in "pooled"

## Merge parameters (Only for paired-ends)
dada2_minoverlap <- 12 # Default: 12
dada2_maxmismatch <- 0 # Default: 0
dada2_concatenate <- FALSE # Default: FALSE

## Removing Quimeras parameters
dada2_method <- "consensus" # Default:"consensus" if dada2_pool -> TRUE, set in "pooled"

## Taxonomical assignment
dada2_path_ref_fasta <- "~/Programas/reference_databases/silva_nr99_v138.1_train_set.fa.gz"
dada2_tryrc <- FALSE
dada2_path_ref_fasta_species <- "~/Programas/reference_databases/silva_species_assignment_v138.1.fa.gz"
dada2_tryrc_species <- FALSE
