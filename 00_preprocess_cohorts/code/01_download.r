##### Download RNA 16S from ENA #####

# 1.Download enaBrowserTools from: https://github.com/enasequence/enaBrowserTools

# 2.Download FASTQC from: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# 4.Create new directory to allocate 16S Fastq
project <- "PRJNA208535"
data <- readRDS(file = paste0("extdata/",project,"/",project,"_metadata.rds"))
samples <- rownames(data)
dir.create(path = paste0("extdata/", project,"/"))
setwd(dir <- paste0("extdata/", project,"/"))

# 5.Declare path to function enaDataGet 
path_ena <- "/Users/diego/Programas/enaBrowserTools/python3/enaDataGet" 
# 5.1.Run in terminal in ordere to obtain samples
for (i in seq_along(samples)) {
  system(command = paste(path_ena, "-f fastq -m", samples[i]))
}

# 6.Also check profile quality
# 6.1.Declare paths to Fastq and Fastqc program
path_fastq <- "."
path_fastqc <- "~/Programas/FastQC/fastqc"
# 6.2.List all samples and run recursive in terminal
l <- list.files(path_fastq, pattern = "fastq.gz", include.dirs = TRUE, recursive = TRUE, full.names = TRUE)
cmd <- list()
for (i in seq_along(l)){
  cmd[[i]] <- paste0(path_fastqc, " ", l[[i]])
}
for (i in seq_along(cmd)){
  system(cmd[[i]])
}