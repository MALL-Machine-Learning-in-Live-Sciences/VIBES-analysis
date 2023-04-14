# Construct phyloseq object
library(phyloseq)
library(dplyr)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# 1.Prepare paths and load required data
project_name <- "PRJNA208535"
path <- paste0("../extdata/", project_name, "/")
metadata <- readRDS(file = paste0(path, project_name,
                                  "_metadata.rds"))
tax <- readRDS(file = paste0(path, project_name, "_tax_table.rds"))
otu <- readRDS(file = paste0(path, project_name, "_otu_table.rds"))

# 2.Construction of phyloseq object
otu <- phyloseq::t(otu)
ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(tax))

# 3.Shorten the name of our ASVs and store sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# 4.Correctly rename the species for a correct merge later on.
tax <- data.frame(tax_table(ps))
tax$Species[is.na(tax$Species)] <- 0
tax <- tax %>%
  mutate(Species = ifelse(Species != 0, paste0(Genus, "_", Species), Species))
tax <- tax %>%
  mutate(Species = ifelse(Species == 0, NA, Species))
tax <- as.matrix(tax)
tax_table(ps) <- tax

# 5.Save phyloseq object
saveRDS(object = ps, file = paste0("../data/pseqs/", substr(x = project_name, start = 1,
 stop = 9), "_pseq.rds"))
