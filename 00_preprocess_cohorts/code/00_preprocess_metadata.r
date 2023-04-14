library(openxlsx)
# Ravel con muestreas de 2013 y 2021 (6593 muestras) Single V1-V3 region 454 (510 pb)
PRJNA208535 <- data.frame(read.delim(url("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA208535&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true&limit=0"), header = TRUE))
PRJNA208535$first_created <- substr(PRJNA208535$first_created, 1, 4)
PRJNA208535<- PRJNA208535[PRJNA208535$first_created == "2013",]
maintain <- c("run_accession","sample_title", "sample_alias", "library_layout",
              "instrument_platform", "instrument_model", "first_created")
PRJNA208535 <- PRJNA208535[,maintain]
# Unificamos los identificadores para hacer el merge
PRJNA208535$sample_title <- str_sub(PRJNA208535$sample_title, 22,-1)
PRJNA208535$sample_title <- str_replace_all(string = PRJNA208535$sample_title, "B00", "s" )
PRJNA208535$sample_title <- str_replace_all(string = PRJNA208535$sample_title, "B0", "s" )
PRJNA208535$sample_title <- str_replace_all(string = PRJNA208535$sample_title, "B", "s" )
PRJNA208535_CD <- read.xlsx("~/git/BV_Microbiome/extdata/PRJNA208535/40168_2013_28_MOESM1_ESM.xlsx", startRow = 4, sheet = 1)
PRJNA208535_CD$sampleID <- str_replace_all(string = PRJNA208535_CD$sampleID, ".w", "_" )
PRJNA208535_CD$sampleID <- str_replace_all(string = PRJNA208535_CD$sampleID, "d", "_" )
metadata <- merge(x = PRJNA208535, y = PRJNA208535_CD, by.x = "sample_title", by.y = "sampleID")
metadata <- data.frame(metadata[,-2], row.names=metadata[,2])
saveRDS(object = metadata, file = "~/git/BV_Microbiome/extdata/PRJNA208535/PRJNA208535_metadata.rds")


# Ravel 2023 Paired Illumina 
PRJNA797778 <- data.frame(read.delim(url("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA797778&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true&limit=0"), header = TRUE)) # No existe en el ena
PRJNA797778<- PRJNA797778[PRJNA797778$library_strategy == "AMPLICON",]
PRJNA797778$first_created <- substr(PRJNA797778$first_created, 1, 4)
maintain <- c("run_accession","sample_title", "sample_alias", "library_layout",
              "instrument_platform", "instrument_model", "first_created")
PRJNA797778 <- PRJNA797778[,maintain]
samples <- PRJNA797778$run_accession
PRJNA797778_subset <- PRJNA797778[grep("_W", PRJNA797778$sample_alias), ]
PRJNA797778_subset <- data.frame(PRJNA797778_subset[,-1], row.names=PRJNA797778_subset[,1])
saveRDS(object = PRJNA797778_subset, file = "~/git/BV_Microbiome/extdata/PRJNA797778/PRJNA797778_metadata.rds")

# Cohorte fármacos
PRJNA302078 <- data.frame(read.delim(url("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA302078&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true&limit=0"), header = TRUE))
samples <- PRJNA302078$run_accession
PRJNA302078$first_created <- substr(PRJNA302078$first_created, 1, 4)
maintain <- c("run_accession","sample_title", "sample_alias", "library_layout",
              "instrument_platform", "instrument_model", "first_created")
PRJNA302078 <- PRJNA302078[,maintain]
d0_subset <- PRJNA302078[grep("D0", PRJNA302078$sample_alias), ]
sample_names <- strsplit(PRJNA302078$sample_alias, "D", 1)
PRJNA302078$sample_ID <- sapply(strsplit(basename(PRJNA302078$sample_alias), "D"), `[`, 1)
cd <- read.delim(file = "~/git/BV_Microbiome/extdata/PRJNA302078/PRJNA302078_mdata.csv", header = TRUE,sep = ";")
metadata <- merge(x = cd, y = PRJNA302078, by.x = "ID", by.y = "sample_ID")
metadata <- data.frame(metadata[,-3], row.names=metadata[,3])
saveRDS(object = metadata, file = "~/git/BV_Microbiome/extdata/PRJNA302078/PRJNA302078_metadata.rds")