require(mlr3)
require(mlr3measures)
require(data.table)

# Measures
measures <- list(
	msr("classif.auc", id = "AUC"),
	msr("classif.mbrier", id = "Brier Error"))

# Path to results
res_dir <- "05_preterm/res/BV_Microbiome_preterm/"

# See CV results
models_files <- list.files(res_dir, pattern = "rsmp_")
cv <- lapply(models_files, function(x) {
  
  model <- readRDS(file.path(res_dir, x))
  pred <- model$result$prediction()
  
  # See performance in outer CV (aggregate)
  # thold <- 0.5
  # threshold <- c("responder" = thold,
  #                "non-responder" = 1 - thold)
  # pred <- pred$set_threshold(threshold = threshold)

  scores <- pred$score(measures = measures)
  
  res_cv <- data.frame(
    AUC = scores[1],
    Brier = scores[2],
    algorithm = sapply(strsplit(x, "_"), "[", 2),
    data = sapply(strsplit(gsub(".rds", "", x), "_"), "[", 3),
    transformation = sapply(strsplit(gsub(".rds", "", x), "_"), "[", 4))

  return(res_cv)
})
cv <- data.table::rbindlist(cv)

toplot <- melt(cv)
saveRDS(
  toplot, 
  file.path(res_dir, "summary_performances.rds")
  )
