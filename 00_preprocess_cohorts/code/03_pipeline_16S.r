##### DADA2 16S Single-end Fastq Processing  #####
# 0.Load packages
library(dada2)

# 1.Set paths and load scripts
setwd(dir = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/")
setwd("~/git/BV_Microbiome/")
source(file = "00_preprocess_cohorts/code/config_file.r")

# 2.Declare path to fastq and retain samples names
l <- list.files(path = input_dir_path, pattern = pattern)
l_fastq_fs <- list()
for (i in seq_along(l)) {
  l_fastq_fs[i] <- paste0(input_dir_path, l[i], "/",
                          list.files(path = paste0(input_dir_path, l[i]),
                                     pattern = ".fastq.gz"))
}
fnfs <- unlist(l_fastq_fs, use.names = FALSE)
sample_names <- sapply(strsplit(basename(fnfs), ".fastq.gz"), `[`, 1)

# 3.Create a folder for place filtered files
filter_path <- paste(input_dir_path, "Filtered_FASTQ", sep = "")
if (dir.exists(filter_path) == FALSE) {
  dir.create(filter_path)
  message(paste("Creating", filter_path, "directory!"))
}
filtfs <- file.path(filter_path, paste0(sample_names, "_F_filt.fastq.gz"))
names(filtfs) <- sample_names

# 4.Apply quality filters on sequences
out <- filterAndTrim(fwd = fnfs, filt = filtfs, truncLen = dada2_trunclen,
                     maxN = dada2_maxn, maxEE = dada2_maxee,
                     truncQ = dada2_truncq, trimLeft = dada2_trimleft,
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# 5.Learning nucleotide errors from sequences
errf <- learnErrors(filtfs, multithread = TRUE, nbases = dada2_nbases)
plotErrors(errf, nominalQ = TRUE)

# 6.Aplying core algorithm to infer real biological sequences
# and construct sequence table
dadafs <- dada(filtfs, err = errf, multithread = TRUE,
               HOMOPOLYMER_GAP_PENALTY = dada2_homopolymer_gap_penalty,
               BAND_SIZE = dada2_band_size)
seqtab <- makeSequenceTable(dadafs)
## Inspect dimensions of sequences before remove quimeras
print("Dimensions of sequence table before quimeras removal:")
dim(seqtab)
d1 <- dim(seqtab)
print("Distribution of sequence lengths before quimeras removal:")
table(nchar(getSequences(seqtab)))

# 7.Remove quimeras from sequences
seqtab_nochim <- removeBimeraDenovo(seqtab, method = dada2_method,
                                    multithread = TRUE, verbose = TRUE)
## Inspect dimensions of sequences after remove quimeras
("Dimensions of sequence table after quimeras removal:")
dim(seqtab_nochim)
d2 <- dim(seqtab_nochim)
print("Distribution of sequence lengths after quimeras removal:")
table(nchar(getSequences(seqtab_nochim)))
k <- sum(seqtab_nochim) / sum(seqtab)
k
d3 <- d2[2] / d1[2]
print(paste((1 - round(d3, digits = 3)) * 100,
            "% of the sequences were quimeras"))
print(paste("The abundance of these quimeras only represent",
            (1 - round(x = k, digits = 3)) * 100,
            "% of the total abundance"))

# 8.Showing evolution of sequences from raw to final step
getn <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadafs, getn), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample_names
print("Summary of the first 5 samples along the Pipeline:")
head(track)
saveRDS(object = track, file = paste0(out_path,"/", experiment_name,"_summary.rds" ))

# 9.Taxonomical assignment
taxa <- assignTaxonomy(seqs = seqtab_nochim, refFasta = dada2_path_ref_fasta,
                       tryRC = dada2_tryrc, multithread = TRUE)
taxa <- addSpecies(taxtab = taxa, refFasta = dada2_path_ref_fasta_species,
                   tryRC =  dada2_tryrc_species)

# 10.Save ASV table and taxa table
tax_table <- taxa
otu_table <- seqtab_nochim
saveRDS(object = tax_table, file = paste0(out_path,"/", experiment_name,
                                          "_tax_table.rds"))
saveRDS(object = otu_table, file = paste0(out_path,"/", experiment_name,
                                          "_otu_table.rds"))
