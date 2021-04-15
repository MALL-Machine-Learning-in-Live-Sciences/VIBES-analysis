# Functions Data Filter

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

prune.OTUs = function(phyobject, Rank, type, pctg = 0.05, count = 0, vari = 0){ 
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
