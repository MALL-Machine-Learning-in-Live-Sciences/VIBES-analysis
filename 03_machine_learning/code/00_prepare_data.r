# Prepare data for mlr
# 1.Load cohort pseqs
cohort <- "PRJNA3020"
cl_pseq <- readRDS(file = paste0("02_cluster/data/C4/", cohort,
                                 "_Cluster_Species_22_pseq.rds"))
raw_pseq <- readRDS(file = paste0("00_preprocess_cohorts/data/SpeciesIntersect/",
                                  cohort,"_Species_pseq_22.rds"))
# 2.Extract datasets from pseq
get_dataset <- function(phyobject){
  require(phyloseq)
  clinics = as.data.frame(sample_data(phyobject))
  otus = as.data.frame(get_taxa(phyobject))
  identical(rownames(clinics), rownames(otus))
  data = cbind(otus, clinics)
  return(data)
}
# 2.1 Extract species and cluster feture
cl_df = get_dataset(phyobject = cl_pseq)
maintain <- colnames(cl_df[1:22])
maintain <- c(maintain,"cluster")
cl_df <- cl_df[, maintain]

# 2.2.Add cluster feature and retain with speceis only
raw_df = get_dataset(phyobject = raw_pseq)
identical(rownames(cl_df), rownames(raw_df))
raw_df$cluster <- cl_df$cluster
raw_df <- raw_df[, maintain]

# 3.Save dataframes
saveRDS(object = cl_df,
        file = paste0("03_machine_learning/data/", cohort, "_clr.rds"))
saveRDS(object = raw_df,
        file = paste0("03_machine_learning/data/", cohort, "_raw.rds"))