require(phyloseq)
ps <- readRDS(file = "~/git/BV_Microbiome/00_preprocess_cohorts/data/pseqs/PRJNA3934_pseq.rds")
# 1. Filter data
df <-  tibble::as_tibble(ps@sam_data, rownames = NA)
# We want to keep samples below 37 weeks and not postpartum.
muestras <- df %>%
  group_by(SubjectID) %>%
  filter(GWcoll < 37, GWcoll <= GWdel) %>%
  slice_max(GWcoll) %>%
  ungroup()
ps_filtered <- phyloseq::prune_samples(x = ps,
                                       samples = muestras$SampleID)
# Save filtered phyloseq
saveRDS(object = ps_filtered,
        file = "~/git/BV_Microbiome/00_preprocess_cohorts/data/pseqs/PRJNA3934_filtered_pseq.rds")
