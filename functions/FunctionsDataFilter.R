# Functions Data Filter
get.phylo.Ravel = function(otu, clin, taxonomyTable,path, id){
  require(phyloseq)
  require(dplyr)
  ## Retain only the access numbers
  otu$X.OTU.ID = as.vector(otu$X.OTU.ID)
  splitted = strsplit(otu$X.OTU.ID, '_')
  ncbi = list()
  for (i in seq_along(splitted)) {
    ncbi[[i]] = paste(splitted[[i]][1], splitted[[i]][2], sep = '_')
  }
  ncbi = unlist(ncbi)
  otu$X.OTU.ID = ncbi
  
  ## Categorize in: low, intermediate, high 
  high = clin[which(clin$Var > 6),]
  low = clin[which(clin$Var < 4),]
  high$Var = "High"
  low$Var = "Low"
  clinics = rbind(high, low)
  clinics = arrange(clinics, X.SampleID)
  clinics = data_frame(clinics)
  otu = data_frame(otu)
  ## Phyloseq format
  clinics <- clinics %>%
    tibble::column_to_rownames("X.SampleID") 
  
  otu <- otu %>%
    tibble::column_to_rownames("X.OTU.ID") 
  
  tax_mat <- as.matrix(taxonomyTable)
  otu_mat <- as.matrix(otu)
  
  # Make the main phyloseq object
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(clinics)
  BV_phyloseq <- phyloseq(OTU, TAX, samples)
  
  saveRDS(BV_phyloseq, file = paste0(path,id,"_phyloseq.rds"))
  return(BV_phyloseq)
}

get.phylo.Sriniv = function(otu, clin, taxonomyTable,target,path,){
  require(phyloseq)
  require(dplyr)
  ## Categorize in: low, intermediate, high 
  clin$X.SampleID = paste0("X",clin$X.SampleID)
  
  if (target == "Nugent"){
    clin <- subset( clin, select = -Amsel )
    high = clin[which(clin$Nugent > 6),]
    low = clin[which(clin$Nugent < 4),]
    high$Nugent = "High"
    low$Nugent = "Low"
    clinics = rbind(high, low)
    clinics = arrange(clinics, X.SampleID)
  } else if (target == "Amsel"){
    clin <- subset( clin, select = -Nugent )
    clinics = clin[-1,]
  }
  clinics = data_frame(clinics)
  otu = data_frame(otu)
  
  ## Phyloseq format
  clinics <- clinics %>%
    tibble::column_to_rownames("X.SampleID") 
  
  otu <- otu %>%
    tibble::column_to_rownames("X.OTU.ID") 
  
  tax_mat <- as.matrix(taxonomyTable)
  otu_mat <- as.matrix(otu)
  
  # Make the main phyloseq object
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(clinics)
  BV_phyloseq <- phyloseq(OTU, TAX, samples)
  
  saveRDS(BV_phyloseq, file = paste0(path,"Sriniv_phyloseq.rds"))
  return(BV_phyloseq)
}


phy.aglomerate = function(phyobject, rank){
  require(phyloseq)
  phy =  tax_glom(phyobject, rank)
  return(phy)
}

relat.abun = function(phyobject){
  require(phyloseq)
  phy = transform_sample_counts(phyobject, function(x) x / sum(x))
  return(phy)
}

prune.OTUs = function(phyobject, Rank, pctg = 0.05, count = 0, vari = 0){ 
  require(phyloseq)
  if (Rank == "Rank1"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank1 %in% c("k__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
    
  } else if(Rank == "Rank2"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank2 %in% c("p__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
    
  } else if(Rank == "Rank3"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank3 %in% c("c__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
    
  } else if(Rank == "Rank4"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank4 %in% c("o__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
    
  } else if(Rank == "Rank5"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank5 %in% c("f__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
    
  }else if (Rank == "Rank6"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank6 %in% c("g__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)

  }else if(Rank == "Rank7"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Rank7 %in% c("s__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
  }
  return(phy)
}
