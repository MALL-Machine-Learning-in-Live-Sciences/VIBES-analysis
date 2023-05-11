require(mlr3)
require(mlr3measures)
require(data.table)

# Measures
measures <- list(msr("classif.acc", id = "Accuracy"),
                 msr("classif.auc", id = "AUCROC"),
                 msr("classif.prauc", id = "PRAUC"),
                 msr("classif.sensitivity", id = "Sensitivity"),
                 msr("classif.specificity", id = "Specificity"))