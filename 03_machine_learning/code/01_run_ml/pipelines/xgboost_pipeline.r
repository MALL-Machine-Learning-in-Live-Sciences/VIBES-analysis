# RandomForest function
# =======
setwd(here::here())
source("code/utils/build_learners.r")
source("code/utils/pipeline_utils.r")

xgboost_pipeline <- function(data,
                        dataname,
                        target,
                        positive,
                        removeConstant,
                        normalize,
                        filterFeatures,
                        inner,
                        outer,
                        measure,
                        method_at,
                        method_afs,
                        term_evals,
                        n_evals_afs,
                        fselector,
                        workers,
                        outDir,
                        parallel,
                        resampling,
                        folds,
                        batch_size,
                        seed) {
  set.seed(seed)
  # Make task
  task <- making_task(data,
                      dataname,
                      target,
                      positive)
  # Preprocess
  task <- preprocess(task,
                     removeConstant,
                     normalize,
                     filterFeatures)
  # Learner
  learner <- xgboost(inner,
                     measure,
                     method_at,
                     method_afs,
                     term_evals,
                     n_evals_afs,
                     fselector)

  if (resampling == TRUE) {
    id <- "rsmp"
    # Parallelization
    if (parallel == TRUE) {
            future::plan(list(
                future::tweak("multisession", workers = folds),  # outer
                future::tweak("multisession", workers = batch_size))) # inner
    }
    # Resampling
    rr <- resample(task,
                   learner,
                   resampling = outer,
                   store_models = TRUE)
  } else {
    id <- "autotuner"
    # Parallelization
    if (parallel == TRUE) {
            future::plan(list(
                future::tweak("multisession", workers = batch_size))) # inner
    }
    # Autotuner
    rr <- learner$train(task)
  }

  # Save resampling object
  res <- list(task = task,
              result = rr)
  saveRDS(res,
          file = paste0(outDir, "/", id, "_xgboost_", dataname, ".rds"))

}