setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# 1.Descargar el enaBrowserTools desde el siguiente link: https://github.com/enasequence/enaBrowserTools

# 2.Descargar FASTQC desde el siguiente link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# 3.Filtrado de las muestras de 2013 del PRJNA208535
PRJNA208535 <- data.frame(read.delim(url("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA208535&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,submission_accession,tax_id,scientific_name,instrument_platform,instrument_model,library_name,nominal_length,library_layout,library_strategy,library_source,library_selection,read_count,base_count,center_name,first_public,last_updated,experiment_title,study_title,study_alias,experiment_alias,run_alias,fastq_bytes,fastq_md5,fastq_ftp,fastq_aspera,fastq_galaxy,submitted_bytes,submitted_md5,submitted_ftp,submitted_aspera,submitted_galaxy,submitted_format,sra_bytes,sra_md5,sra_ftp,sra_aspera,sra_galaxy,cram_index_ftp,cram_index_aspera,cram_index_galaxy,sample_alias,broker_name,sample_title,nominal_sdev,first_created&format=tsv&download=true&limit=0"), header = TRUE))
PRJNA208535$first_created <- substr(PRJNA208535$first_created, 1, 4)
PRJNA208535<- PRJNA208535[PRJNA208535$first_created == "2013",]
maintain <- c("run_accession","sample_title", "sample_alias", "library_layout",
              "instrument_platform", "instrument_model", "first_created")
PRJNA208535 <- PRJNA208535[, maintain]
# Las muestras de Ravel de 2013  son single, regiones V1-V3 y tecnología 454
# Por las regiones tendríamos un tamaño total de ~510 pb.
# Para descargar los fastq lo que nos interesa es el run_accession
samples <- PRJNA208535$run_accession

# 4.Ahora creamos directorio para alojar las muestras descargadas
dir.create(path = "PRJNA208535/")
setwd(dir <- paste0("PRJNA208535"))

# 5.Declaramos el path a la funcion enaDataGet 
path_ena <- "/Users/diego/Programas/enaBrowserTools/python3/enaDataGet" # Cambiar para donde la tengas

# 5.1.Corremos en la terminal donde vamos a descargar cada muestars de "samples"
for (i in seq_along(samples)) {
  system(command = paste(path_ena, "-f fastq -m", samples[i]))
}

# 6. Una vez descargadas habría que comprobar la calidad con el Fastqc
# para saber donde meter el corte. Declaramos path al programa Fastqc y
# a la carpeta donde descargamos los Fastq
path_fastq <- "~/git/BV_Microbiome/extdata/PRJNA302078"
path_fastqc <- "~/Programas/FastQC/fastqc"
# 6.1.Listamos de manera recursiva todos los fastq (Debería haber 1657 sino algo falla)
l <- list.files(path_fastq, pattern = "fastq.gz", include.dirs = TRUE, recursive = TRUE, full.names = TRUE)
# 6.2.Cremamos una lista con los comandos y las muestras para correr fastqc y la corremos
cmd <- list()
for (i in seq_along(l)){
  cmd[[i]] <- paste0(path_fastqc, " ", l[[i]])
}
for (i in seq_along(cmd)){
  system(cmd[[i]])
}