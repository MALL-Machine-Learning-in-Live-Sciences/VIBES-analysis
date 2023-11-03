require(phyloseq)
require(tibble)
require(dplyr)
require(VIBES)
library(stringr)

ps <-  readRDS("00_preprocess_cohorts/data/pseqs/PRJNA3934_filtered_pseq.rds")
vbs <- VIBES::get_clusters(object = ps,
                           column = 7)
df <- as(sample_data(vbs), "data.frame")

df_vibes <- df %>%
  as_tibble() %>%
  dplyr::select(c(SampleID,VCS.I,VCS.II,VCS.III,VCS.IV, Outcome)) 


# Load valencias and select scores
valencias <- as_tibble(read.delim(file = "01_get_valencias/res/PRJNA3934_filtered.csv",
                                  header = TRUE,sep = ","),rownames = NA)
df_val <- valencias %>%
  dplyr::select(c(sampleID,I.A_sim,I.B_sim,II_sim, III.A_sim, III.B_sim,
                  IV.A_sim, IV.B_sim,IV.C0_sim, IV.C1_sim, IV.C2_sim,IV.C3_sim,
                  IV.C4_sim, score)) %>%
  mutate(sampleID = stringr::str_remove(sampleID, "X"))

# Join data
df_all <- merge(x = df_vibes, y = df_val,
                by.x = "SampleID", by.y = "sampleID")
# Prepare data
df_vibes <- df_vibes %>%
  tibble::column_to_rownames(var = "SampleID")

df_all <- df_all %>%
  select(-Outcome, everything(), Outcome) %>%
  tibble::column_to_rownames(var = "SampleID")
  
df_val <- df_all %>%
  select(-VCS.I, -VCS.II, -VCS.III, -VCS.IV )


df_all <- df_all %>%
  tibble::rownames_to_column(var = "row_name") %>%
  arrange(row_name)

df_vibes <- df_vibes %>%
  tibble::rownames_to_column(var = "row_name") %>%
  arrange(row_name)

df_val <- df_val %>%
  tibble::rownames_to_column(var = "row_name") %>%
  arrange(row_name)

identical(df_vibes$row_name, df_val$row_name)
identical(df_vibes$Outcome, df_val$Outcome)
identical(df_vibes$row_name, df_all$row_name)
identical(df_vibes$Outcome, df_all$Outcome)

# Save dataframes
saveRDS(object = df_vibes, file = "05_preterm/data/df_vibes.rds")
saveRDS(object = df_val, file = "05_preterm/data/df_val.rds")
saveRDS(object = df_all, file = "05_preterm/data/df_all.rds")
