#args <- commandArgs(trailingOnly = TRUE)

#cesga <- ifelse(args[1] == "cesga", TRUE, FALSE)
#ExperimentName <- args[2]
#inputDir <- args[3]
#outDir <- args[4]

source("03_machine_learning/code/ml/configFile.r")

if (cesga == TRUE) {
  exec <- "sbatch"
} else {
  exec <- "sh"
}

# Execute in parallel from Cesga
files <- list.files(path = inputDir)
input_algs <- list.files(path = path_algs, pattern = pattern)

for (i in seq_along(files)) {
  for (j in seq_along(input_algs)) {
    new_file <- paste0(
      exec_path, "/",
      gsub(".r", "", input_algs[j], fixed = TRUE),
      "_",
      gsub(".rds", ".sh", files[i]))
    # Creating new file
    sink(new_file)
    cat("#!/bin/bash \n")
    cat(paste("#SBATCH", "-p", part, "\n"))
    cat(paste("#SBATCH", "-t", time, "\n"))
    cat(paste("#SBATCH", paste0("--mem=", mem), "\n"))
    cat(paste("#SBATCH", "-N", nodes, "\n"))
    cat(paste("#SBATCH", "-n", ntasks, "\n"))
    cat(paste("name=", gsub(".rds", "", files[i]), "\n", sep = ""))
    cat(paste("data=", file.path(inputDir, files[i]), "\n", sep = ""))
    cat(paste0("outdir=", outDir, "/", "\n"))
    cat(paste("Rscript ", models_path,
              input_algs[j], " $data", " $name", " $outdir", sep = ""))
    sink(file = NULL)
    system(paste(exec, new_file))
  }
}
