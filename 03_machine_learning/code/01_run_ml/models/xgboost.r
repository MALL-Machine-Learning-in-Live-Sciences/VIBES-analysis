setwd(here::here())
source("requirements.r")
source("code/configFile.r")
source("code/pipelines/xgboost_pipeline.r")

args <- commandArgs(trailingOnly = TRUE)

data <- readRDS(args[1])
names(data) <- make.names(names(data))
set.seed(seed)

xgboost_pipeline(
    data = data,
    dataname = args[2],
    target = target,
    positive = positive,
    removeConstant = removeConstant,
    normalize = normalize,
    filterFeatures = filterFeatures,
    inner = inner,
    outer = outer,
    measure = measure,
    method_at = method_at,
    method_afs = method_afs,
    term_evals = term_evals,
    n_evals_afs = n_evals_afs,
    fselector = fselector,
    workers = ntasks,
    outDir = args[3],
    parallel = parallel,
    resampling = resampling,
    folds = folds,
    batch_size = batch_size,
    seed = seed
)