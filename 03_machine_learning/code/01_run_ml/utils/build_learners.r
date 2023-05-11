setwd(here::here())
source("code/utils/tuning_utils.r")

# Random Forest
# ====
randomForest <- function(inner,
                         measure,
                         method_at,
                         method_afs,
                         term_evals,
                         n_evals_afs,
                         fselector) {
  # Make learner
  learner <- lrn("classif.ranger",
                  num.trees = 1000,
                  num.threads = 10,
                  importance = "impurity",
                  predict_type = "prob")
  learner$encapsulate <- c(train = "evaluate", predict = "evaluate")
  learner$fallback <- lrn("classif.svm", predict_type = "prob")
  # Hyperparameter space
  ps <- ps(
    mtry = p_int(lower = 3, upper = 15),
    min.node.size = p_int(lower = 1, upper = 5),
    max.depth = p_int(lower = 2, upper = 15),
    alpha = p_dbl(lower = 0, upper = 1)
  )
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    n_evals_afs,
                    method_at,
                    fselector,
                    method_afs)
  return(afs)
}


# Glmnet
# ===
glmnet <- function(inner,
                   measure,
                   method_at,
                   method_afs,
                   n_evals_afs,
                   term_evals,
                   fselector) {
  # Make learner
  learner <- lrn("classif.glmnet",
                predict_type = "prob")
  learner$encapsulate <- c(train = "evaluate", predict = "evaluate")
  learner$fallback <- lrn("classif.svm", predict_type = "prob")
  # Hyperparameter space
  ps <- ps(
    alpha = p_dbl(lower = 0, upper = 1),
    s = p_dbl(lower = 0, upper = 1)
  )
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    n_evals_afs,
                    method_at,
                    fselector,
                    method_afs)
  return(afs)
}


# SVM
# ===
svm <- function(inner,
                measure,
                method_at,
                method_afs,
                term_evals,
                n_evals_afs,
                fselector) {
  # Make learner
  learner <- lrn("classif.svm",
                 predict_type = "prob")
  # Hyperparameter space
  ps <- ps(
    cost = p_dbl(lower = 1e-8, upper = 1e8, logscale = FALSE),
    gamma = p_dbl(lower = 1e-8, upper = 1e8, logscale = FALSE),
    kernel = p_fct(levels = c("polynomial", "radial", "sigmoid")),
    type = p_fct(levels = "C-classification")
  )
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    n_evals_afs,
                    method_at,
                    fselector,
                    method_afs)
  return(afs)
}


# XGBoost
xgboost <- function(inner,
                    measure,
                    method_at,
                    method_afs,
                    term_evals,
                    n_evals_afs,
                    fselector) {
  # Make learner
  learner <- lrn("classif.xgboost",
                 nthread = 10,
                 predict_type = "prob")
  # Hyperparameter space
  ps <- ps(
    eta = p_dbl(lower = 0.01, upper = 0.3),
    max_depth = p_int(lower = 3, upper = 10)
  )
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    n_evals_afs,
                    method_at,
                    fselector,
                    method_afs)
  return(afs)
}


# Ligth GBM
lgbm <- function(inner,
                 measure,
                 method_at,
                 method_afs,
                 term_evals,
                 n_evals_afs,
                 fselector) {
  # Make learner
  learner <- lrn("classif.lightgbm",
                 num_threads = 10,
                 predict_type = "prob")
  # Hyperparameter space
  ps <- ps(
    learning_rate = p_dbl(lower = 0.01, upper = 0.3),
    max_depth = p_int(lower = 5, upper = 10)
  )
  # Hyperparameters and features tuner
  afs <- make_tuner(inner,
                    measure,
                    learner,
                    ps,
                    term_evals,
                    n_evals_afs,
                    method_at,
                    fselector,
                    method_afs)
  return(afs)
}