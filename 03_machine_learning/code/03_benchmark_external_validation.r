require(dplyr)
require(pbapply)
require(mlr3)
require(mlr3measures)
require(data.table)

setwd("/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/run_mlr3/results/BV_Microbiome")
test_dir <- "/mnt/netapp2/Store_uni/home/ulc/co/jlb/git/run_mlr3/data/BV_Microbiome/test"

# Measures
measures <- list(
    msr("classif.acc", id = "Accuracy"),
    msr("classif.mbrier", id = "Brier Error"),
    msr("classif.bacc", id = "Balanced Accuracy"))

files <- list.files(pattern = "rsmp")

# Load models and make predictions
p <- list()
for (i in seq_along(files)) {

    file <- files[i]
    transformation <- gsub(".rds", "", sapply(strsplit(file, "_"), "[", 4))
    model <- readRDS(file)

    # Train the model with all observations
    set.seed(1993)
    model_trained <- model$result$learner$train(model$task)

    # Load test datasets
    test_files <- list.files(test_dir, pattern = transformation)
    test <- lapply(
        test_files,
        function(x) {
            test_data <- readRDS(file.path(test_dir, x))
            test_data$cluster <- as.factor(test_data$cluster)
            return(test_data)
        })
    names(test) <- test_files

    # Predict on new data
    pred1 <- model_trained$predict_newdata(test[[1]])
    pred2 <- model_trained$predict_newdata(test[[2]])
    pred3 <- model_trained$predict_newdata(test[[3]])
    pred4 <- model_trained$predict_newdata(test[[4]])

    # Save metrics in a data frame
    p1 <- as.data.frame(pred1$score(measures=measures))
    p2 <- as.data.frame(pred2$score(measures=measures))
    p3 <- as.data.frame(pred3$score(measures=measures))
    p4 <- as.data.frame(pred4$score(measures=measures))

    p[[i]] <-
        cbind.data.frame(p1, p2, p3, p4) %>%
        t() %>%
        as.data.frame() %>%
        mutate(
            algorithm = model$result$learner$id,
            data = sapply(strsplit(test_files, "_"), "[", 1),
            transformation = transformation)

}

# Join all predictions
ext_val <- data.table::rbindlist(p)
ext_val <- melt(ext_val)
