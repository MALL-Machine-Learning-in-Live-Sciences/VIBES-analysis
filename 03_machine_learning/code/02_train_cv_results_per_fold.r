require(ggplot2)
require(mlr3)
require(mlr3measures)
require(data.table)
require(ggpubr)
# Measures
measures <- list(
  msr("classif.bacc", id = "Balanced Accuracy"),
  msr("classif.auc", id = "AUC"),
  msr("classif.mbrier", id = "Brier Error"))

res_dir <- "04_treatment/res/BV_Microbiome_treatment/"
models_files <- list.files(res_dir, pattern = "rsmp_")


cv <- lapply(models_files, function(x) {
  
  model <- readRDS(file.path(res_dir, x))
  pred <- model$result$predictions()
  # Make a lsit to extract each fold
  prl <-  list()
  for (i in seq_along(pred)) {
    scores <- pred[[i]]$score(measures = measures)
    prl[[i]] <- scores
  }
aa <- unlist(prl,recursive = TRUE)
variable <- names(aa)
value <- as.vector(aa)
# Asign each measure to each fold
folds <- rep(paste0("Fold_",1:10), each = 3)
# Extracta algorithm names
algorithm <- sapply(strsplit(x, "_"), "[", 2)
algorithm <- rep(x = algorithm, 30)
# Extract data naames
data <- sapply(strsplit(gsub(".rds", "", x), "_"), "[", 4)
data <- rep(x = data, 30)

# Build dataframe with all measures
res_c <- data.frame(cbind(algorithm,
           folds,
           data,
           variable,
           value))
  
  return(res_c)
})
cvt <- data.table::rbindlist(cv)
cvt <-  cvt %>% as.tibble(cvt) %>%
  mutate_at("value", as.numeric)
saveRDS(object = cvt,
        file = "04_treatment/res/BV_Microbiome_treatment/summary_performances.rds")
