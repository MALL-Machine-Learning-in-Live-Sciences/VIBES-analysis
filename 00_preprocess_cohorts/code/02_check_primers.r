##### Check primers in downloaded FASTQ #####
# 0.Load packages
library(dada2)
library(ShortRead)
library(Biostrings)

# 1.Set paths
setwd("~/git/BV_Microbiome/extdata/")
experiment_name <- "PRJNA208535"
pattern <- "SRR"
input_dir_path <- paste0(experiment_name)

# 2.Determining the primers and get all orients
FWD <- "CCTACGGGAGGCAGCAG"
REV <- "CCGTCAATTCMTTTRAGT"
# 2.1.Function for get all orients of primers 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# 3.Declare path to fastq and retain samples names
l <- list.files(path = input_dir_path, pattern = pattern)
# 3.1.Forward
l_fastq_fs <- list()
for (i in seq_along(l)) {
  l_fastq_fs[i] <- paste0(input_dir_path,"/", l[i], "/",
                          list.files(path = paste0(input_dir_path,"/", l[i]),
                                     pattern = ".fastq.gz"))
}
fnfs <- unlist(l_fastq_fs, use.names = FALSE)
# 3.2.Reverse
l_fastq_rs <- list()
for (i in seq_along(l)) {
  l_fastq_rs[i] <- paste0(input_dir_path,"/", l[i], "/",
                          list.files(path = paste0(input_dir_path,"/", l[i]),
                                     pattern = ".fastq.gz"))
}
fnrs <- unlist(l_fastq_rs, use.names = FALSE)

fnFs <- fnfs
fnRs <- fnrs
# 4. The presence of ambiguous bases (Ns) in the sequencing reads makes accurate
# mapping of short primer sequences difficult, just to remove those with Ns
fnFs.filtN <- file.path(input_dir_path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(input_dir_path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)


# 5.Count the number of times the primers appear in the forward and reverse
# read, while considering all possible primer orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[11]]), 
      
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[11]]))
