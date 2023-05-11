##### DADA2 16S Paired-end Fastq Processing #####
# Load packages
library(dada2)
# Set paths and load scripts
setwd(dir = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome")
#setwd("/Users/diego/git/BV_Microbiome")
source(file = "config_file.r")
t1 = Sys.time()
library(parallel)
parallel::detectCores()

# 1.Declare path to fastq and retain samples names
l <- list.files(path = input_dir_path, pattern = pattern)
## Forward
l_fastq_fs <- list()
for (i in seq_along(l)) {
  path <- paste0(input_dir_path, l[i],"/")
  l_fastq_fs[i] <- paste0(path, list.files(path = path, pattern = "_1.fastq.gz"))
}
fnfs <- unlist(l_fastq_fs, use.names = FALSE)
## Reverse
l_fastq_rs <- list()
for (i in seq_along(l)) {
  path <- paste0(input_dir_path, l[i],"/")
  l_fastq_rs[i] <- paste0(path, list.files(path = path, pattern = "_2.fastq.gz"))
}
fnrs <- unlist(l_fastq_rs, use.names = FALSE)
sample_names <- sapply(strsplit(basename(fnfs), "_"), `[`, 1)

print(paste0("Applying pipeline to ", length(fnfs), " pairs of sequences"))
# 2.Create a folder for place filtered files
filter_path <- paste(input_dir_path, "Filtered_FASTQ", sep = "")
if (dir.exists(filter_path) == FALSE) {
  dir.create(filter_path)
  message(paste("Creating", filter_path, "directory!"))
}
filtfs <- file.path(filter_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtrs <- file.path(filter_path, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtfs) <- sample_names
names(filtrs) <- sample_names

# 3.Apply quality filters on sequences
print("Starting the filtering phase with:")
print("truncLen:")
dada2_trunclen
print("maxEE:")
dada2_maxee 
print("truncQ:")
dada2_truncq
print("trimLeft:")
dada2_trimleft
out <- filterAndTrim(fwd = fnfs, rev = fnrs, filt = filtfs, filt.rev = filtrs,
                      truncLen = dada2_trunclen, maxN = dada2_maxn,
                      maxEE = dada2_maxee, truncQ = dada2_truncq,
                      trimLeft = dada2_trimleft, rm.phix = TRUE,
                      compress = TRUE, multithread = TRUE)
print(paste0(round(sum(out[,2] / out[,1])/nrow(out) * 100, digits=2), "% of reads passed filtering"))
print("Proceed only if the majority (>50%) of reads passed filtering.")

# 4.Learning nucleotide errors from sequences
print("Starting the parametric error model phase with nbases:")
dada2_nbases
errf <- learnErrors(filtfs, multithread = TRUE, nbases = dada2_nbases)
errr <- learnErrors(filtrs, multithread = TRUE, nbases = dada2_nbases)

# 5.Aplying core algorithm to infer real biological sequences
print("Starting the sample inference phase with:")
print("HOMOPOLYMER_GAP_PENALTY:")
dada2_homopolymer_gap_penalty
print("BAND_SIZE:")
dada2_band_size
print("Pool:")
dada2_pool
dadafs <- dada(filtfs, err = errf, multithread = TRUE, pool = dada2_pool,
               HOMOPOLYMER_GAP_PENALTY = dada2_homopolymer_gap_penalty,
               BAND_SIZE = dada2_band_size, verbose = FALSE)
dadars <- dada(filtrs, err = errr, multithread = TRUE, pool = dada2_pool,
               HOMOPOLYMER_GAP_PENALTY = dada2_homopolymer_gap_penalty,
               BAND_SIZE = dada2_band_size, verbose = FALSE)

# 6.Merge paired reads and construct sequence table
print("Starting the merging phase with:")
print("minOverlap:")
dada2_minoverlap
print("maxMismatch:")
dada2_maxmismatch
print("justConcatenate:")
dada2_concatenate
contigs <- mergePairs(dadaF = dadafs, derepF = filtfs, dadaR = dadars,
                      derepR = filtrs, minOverlap = dada2_minoverlap,
                      maxMismatch = dada2_maxmismatch,
                      justConcatenate = dada2_concatenate, verbose = FALSE)
seqtab <- makeSequenceTable(contigs)
## Inspect dimensions of sequences before remove quimeras
print("Dimensions of sequence table before quimeras removal:")
dim(seqtab)
d1 <- dim(seqtab)
print("Distribution of sequence lengths before quimeras removal:")
table(nchar(getSequences(seqtab)))

# 7.Remove quimeras from sequences
print("Starting the chimera elimination phase with:")
print("method:")
dada2_method
seqtab_nochim <- removeBimeraDenovo(seqtab, method = dada2_method,
                      multithread = TRUE, verbose = TRUE)
## Inspect dimensions of sequences after remove quimeras
print("Dimensions of sequence table after quimeras removal:")
dim(seqtab_nochim)
d2 <- dim(seqtab_nochim)
print("Distribution of sequence lengths after quimeras removal:")
table(nchar(getSequences(seqtab_nochim)))
k <- sum(seqtab_nochim) / sum(seqtab)
d3 <- d2[2] / d1[2]
print(paste((1 - round(d3, digits = 3)) * 100,
                      "% of the sequences were quimeras"))
print(paste("The abundance of these quimeras only represent",
                       (1 - round(x = k, digits = 3)) * 100,
                        "% of the total abundance"))

# 8.Showing evolution of sequences from raw to final step
print("Summary of the first samples throughout all phases:")
getn <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadafs, getn), sapply(dadars, getn),
                    sapply(contigs, getn), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                      "merged", "nonchim")
rownames(track) <- sample_names
head(track)

# 9.Taxonomical assignment
print("Performing Taxonomic Assignment")
taxa <- assignTaxonomy(seqs = seqtab_nochim, refFasta = dada2_path_ref_fasta,
                      tryRC = dada2_tryrc, multithread = TRUE)
taxa <- addSpecies(taxtab = taxa, refFasta = dada2_path_ref_fasta_species,
                      tryRC =  dada2_tryrc_species)

# 10.Save OTU table and taxa table
print("Finished, saving taxa and ASVs table")
tax_table <- taxa
otu_table <- seqtab_nochim
saveRDS(object = tax_table, file = paste0(out_path,"/", experiment_name,
                      "_tax_table.rds"))
saveRDS(object = otu_table, file = paste0(out_path,"/", experiment_name,
                      "_otu_table.rds"))
t2 = Sys.time()
runtime = t2-t1
print("Session time")
print(runtime)
