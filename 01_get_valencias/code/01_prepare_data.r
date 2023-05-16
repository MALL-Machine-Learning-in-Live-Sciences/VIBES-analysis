##### Prepare data to compute Valencias CSTs #####
# 0.Load packages
require(dplyr)
require(tibble)
require(phyloseq)
require(reshape2)

valencias <- read.delim2(
  "~/git/VALENCIA/CST_centroids_012920.csv",
  header = TRUE,
  sep = ',')

# 1.Load data
setwd("~/git/BV_Microbiome/00_preprocess_cohorts/data/")
cohort <- "Sriniv" #PRJNA2085 PRJNA3020 PRJNA7977 Ravel Sriniv
file <- list.files(pattern = cohort)
pseq <- readRDS(file)

# 2.Extract OTU
otu <- 
  data.frame(otu_table(pseq)) %>% 
  rownames_to_column(var = "ASV") %>% 
  melt()

# 3.Extract taxa
taxa <- 
  data.frame(tax_table(pseq)) %>% 
  rownames_to_column(var = "ASV")

# 4.Join
data <- 
  left_join(otu, taxa, by = "ASV")

# 5.Rename taxonomic values variables
colnames(data)[4:10] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

kk <- grep("^p_", data$Phylum)
length(kk)

if (length(kk) == nrow(data)){
  data <-
    data %>% 
    mutate(
      Kingdom = sapply(strsplit(Kingdom, "__"), "[", 2),
      Phylum = sapply(strsplit(Phylum, "__"), "[", 2),
      Class = sapply(strsplit(Class, "__"), "[", 2),
      Order = sapply(strsplit(Order, "__"), "[", 2),
      Family = sapply(strsplit(Family, "__"), "[", 2),
      Genus = sapply(strsplit(Genus, "__"), "[", 2),
      Species = sapply(strsplit(Species, "__"), "[", 2)
    )
}

data <- 
  data %>% 
  mutate(
    Kingdom = paste0("d_", Kingdom),
    Phylum = paste0("p_", Phylum),
    Class = paste0("c_", Class),
    Order = paste0("o_", Order),
    Family = paste0("f_", Family),
    Genus = paste0("g_", Genus),
    Species = gsub(" ", "_", Species)
  )

# 6.Calculate total counts (= library size)
read_counts <- 
  data %>% 
  group_by(variable) %>% 
  summarise(value = sum(value)) %>% 
  rename(
    sampleID = variable,
    read_count = value
  )

# 7. Grep common features between valencia centroids and data
dom <- names(valencias)[grep("^d_", names(valencias))]
phy <- names(valencias)[grep("^p_", names(valencias))]
cl <- names(valencias)[grep("^c_", names(valencias))]
ord <- names(valencias)[grep("^o_", names(valencias))]
fam <- names(valencias)[grep("^f_", names(valencias))]
gen <- names(valencias)[grep("^g_", names(valencias))]
spec <- setdiff(
  names(valencias)[2:200],
  c(dom, phy, cl, ord, fam, gen))

# 8.Checking!!
l <- length(dom) + length(phy) + length(cl) + length(ord) + length(fam) + length(gen) + length(spec)
stopifnot(l == (ncol(valencias) - 1))
setdiff(dom, data$Kingdom)
setdiff(phy, data$Phylum)
setdiff(cl, data$Class)
setdiff(ord, data$Order)
setdiff(fam, data$Family)
setdiff(gen, data$Genus)
setdiff(spec, data$Species)

# 9.Selecting features
# Kingdom
data_d <- 
  data %>% 
  select(c(ASV, variable, value, Kingdom)) %>% 
  rename(taxa = Kingdom) %>% 
  filter(taxa %in% dom) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

# Phylum
data_p <- 
  data %>% 
  select(c(ASV, variable, value, Phylum)) %>% 
  rename(taxa = Phylum) %>% 
  filter(taxa %in% phy) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

# Class
data_c <- 
  data %>% 
  select(c(ASV, variable, value, Class)) %>% 
  rename(taxa = Class) %>% 
  filter(taxa %in% cl) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

# Order
data_o <- 
  data %>% 
  select(c(ASV, variable, value, Order)) %>% 
  rename(taxa = Order) %>% 
  filter(taxa %in% ord) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

# Family
data_f <- 
  data %>% 
  select(c(ASV, variable, value, Family)) %>% 
  rename(taxa = Family) %>% 
  filter(taxa %in% fam) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

# Genera
data_g <- 
  data %>% 
  select(c(ASV, variable, value,Genus)) %>% 
  rename(taxa = Genus) %>% 
  filter(taxa %in% gen) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

# Species 
data_s <- 
  data %>% 
  select(c(ASV, variable, value, Species)) %>% 
  rename(taxa = Species) %>% 
  filter(taxa %in% spec) %>% 
  group_by(variable, taxa) %>% 
  summarize(value = sum(value)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = "variable",
    names_from = "taxa",
    values_from = "value"
  ) %>% 
  rename(sampleID = variable)

rm(list = setdiff(
  ls(),
  c("cohort", "read_counts", 
    "data_d", "data_p", "data_c", 
    "data_o", "data_f", "data_g",
    "data_s")))

# 10.Join all datasets
data <- 
  left_join(read_counts, data_d, by = "sampleID") %>% 
  left_join(., data_p, by = "sampleID") %>% 
  left_join(., data_c, by = "sampleID") %>% 
  left_join(., data_o, by = "sampleID") %>% 
  left_join(., data_f, by = "sampleID") %>% 
  left_join(., data_g, by = "sampleID") %>% 
  left_join(., data_s, by = "sampleID")

write.csv(
  data, 
  file = file.path(
    "~/git/vaginosis-jlb/00b_get_valencias/data/",
    paste0(cohort, ".csv")),
  row.names = FALSE)
