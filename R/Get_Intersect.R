Ravel = readRDS("~/git/BVMetaGenomics/data/Phyloseqs/Ravel_phyloseq.rds")
Sriniv = readRDS("~/git/BVMetaGenomics/data/Phyloseqs/Sriniv_phyloseq.rds")
#Funciones
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
get.dataset = function(phyobject){
  require(phyloseq)
  clinics = as.data.frame(sample_data(phyobject))
  otus = as.data.frame(t(get_taxa(phyobject)))
  identical(rownames(clinics), rownames(otus))
  data = cbind(otus, clinics)
  return(data)
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
phylo.intersect = function(Ravel, Sriniv, rank){
  get.dataset = function(phyobject){
    require(phyloseq)
    clinics = as.data.frame(sample_data(phyobject))
    otus = as.data.frame(t(get_taxa(phyobject)))
    identical(rownames(clinics), rownames(otus))
    data = cbind(otus, clinics)
    return(data)
  }
  require(phyloseq)
  if (rank == "Rank6"){
    #Intersect both phyloseqs
    Ravel.df = as.data.frame(tax_table(Ravel))
    Sriniv.df = as.data.frame(tax_table(Sriniv))
    inter = intersect(Sriniv.df$Rank6, Ravel.df$Rank6)
    print("Genera shared")
    print(inter)
    # Select OTUs names that belong to those genera
    Ravel.df = Ravel.df[Ravel.df$Rank6 %in% inter,]
    Ravel_names = Ravel.df[order(Ravel.df$Rank6, decreasing = TRUE), ]
    Ravel_names = Ravel_names[6]
    
    Sriniv.df =Sriniv.df[Sriniv.df$Rank6 %in% inter,]
    Sriniv_names = Sriniv.df[order(Sriniv.df$Rank6, decreasing = TRUE), ]
    Sriniv_names = Sriniv_names[6]
    identical(Sriniv_names$Rank6, Ravel_names$Rank6)
    
    # Merge to obtain correspondence betweeen Studies 
    Ravel_names <- cbind(rownames(Ravel_names),Ravel_names)
    colnames(Ravel_names)[1] <- "Ravel_OTUS"
    Sriniv_names <- cbind(rownames(Sriniv_names),Sriniv_names)
    colnames(Sriniv_names)[1] <- "Sriniv_OTUS"
    Names_correspondecy = merge(x = Ravel_names,y = Sriniv_names, all = TRUE)
    NC = Names_correspondecy[order(Names_correspondecy$Sriniv_OTUS, decreasing = FALSE), ]
    
    # Prune phyloseq to maintain only those genera
    Ravel = prune_taxa(Ravel_names$Ravel_OTUS, Ravel)
    Sriniv = prune_taxa(Sriniv_names$Sriniv_OTUS, Sriniv)
    Ravel_data = get.dataset(Ravel)
    Sriniv_data = get.dataset(Sriniv)
    
    # Replace Ravel and Sriniv OTU names with spps names (Rank7)
    # Ravel
    dat = Ravel_data[1:length(NC$Rank6)]
    dat = dat[,sort(names(dat))]
    s = colnames(dat)
    NCR = NC[order(NC$Ravel_OTUS),]
    identical(s, NCR$Ravel_OTUS)
    colnames(dat) = NCR$Rank6
    colnames(dat)  = substr(colnames(dat) , start = 4, stop = 100)
    Ravel.df = cbind(dat, Ravel_data[(length(NC$Rank6)+1):length(names(Ravel_data))])
    
    #Srinivasan
    dat2 = Sriniv_data[1:length(NC$Rank6)]
    dat2 = dat2[,sort(names(dat2))]
    s2 = colnames(dat2)
    NCS = NC[order(NC$Sriniv_OTUS),]
    identical(s2, NCS$Sriniv_OTUS)
    colnames(dat2) = NCS$Rank6
    colnames(dat2)  = substr(colnames(dat2) , start = 4, stop = 100)
    Sriniv.df = cbind(dat2, Sriniv_data[(length(NC$Rank6)+1):length(names(Sriniv_data))])
    
  }else if(rank == "Rank7"){
    #Intersect both phyloseqs
    Ravel.df = as.data.frame(tax_table(Ravel))
    Sriniv.df = as.data.frame(tax_table(Sriniv))
    inter = intersect(Sriniv.df$Rank7, Ravel.df$Rank7)
    print("Species shared")
    print(inter)
    
    # Select OTUs names that belong to those Species
    Ravel.df = Ravel.df[Ravel.df$Rank7 %in% inter,]
    Ravel_names = Ravel.df[order(Ravel.df$Rank7, decreasing = TRUE), ]
    Ravel_names = Ravel_names[7]
    
    Sriniv.df =Sriniv.df[Sriniv.df$Rank7 %in% inter,]
    Sriniv_names = Sriniv.df[order(Sriniv.df$Rank7, decreasing = TRUE), ]
    Sriniv_names = Sriniv_names[7]
    identical(Sriniv_names$Rank7, Ravel_names$Rank7)
    
    # Merge to obtain correspondence betweeen Studies 
    Ravel_names <- cbind(rownames(Ravel_names),Ravel_names)
    colnames(Ravel_names)[1] <- "Ravel_OTUS"
    Sriniv_names <- cbind(rownames(Sriniv_names),Sriniv_names)
    colnames(Sriniv_names)[1] <- "Sriniv_OTUS"
    Names_correspondecy = merge(x = Ravel_names,y = Sriniv_names, all = TRUE)
    NC = Names_correspondecy[order(Names_correspondecy$Sriniv_OTUS, decreasing = FALSE), ]
    
    # Prune phyloseq to maintain only those species
    Ravel = prune_taxa(Ravel_names$Ravel_OTUS, Ravel)
    Sriniv = prune_taxa(Sriniv_names$Sriniv_OTUS, Sriniv)
    Ravel_data = get.dataset(Ravel)
    Sriniv_data = get.dataset(Sriniv)
    
    # Replace Ravel and Sriniv OTU names with spps names (Rank7)
    # Ravel
    dat = Ravel_data[1:length(NC$Rank7)]
    dat = dat[,sort(names(dat))]
    s = colnames(dat)
    NCR = NC[order(NC$Ravel_OTUS),]
    identical(s, NCR$Ravel_OTUS)
    colnames(dat) = NCR$Rank7
    colnames(dat)  = substr(colnames(dat) , start = 4, stop = 100)
    Ravel.df = cbind(dat, Ravel_data[(length(NC$Rank7)+1):length(names(Ravel_data))])
    
    #Srinivasan
    dat2 = Sriniv_data[1:length(NC$Rank7)]
    dat2 = dat2[,sort(names(dat2))]
    s2 = colnames(dat2)
    NCS = NC[order(NC$Sriniv_OTUS),]
    identical(s2, NCS$Sriniv_OTUS)
    colnames(dat2) = NCS$Rank7
    colnames(dat2)  = substr(colnames(dat2) , start = 4, stop = 100)
    Sriniv.df = cbind(dat2, Sriniv_data[(length(NC$Rank7)+1):length(names(Sriniv_data))])
  }else{
    print("Insert valid rank (6 or 7)")
  }
  data_pool <- list(Ravel.df,Sriniv.df)
  c = c("Ravel", "Sriniv")
  names(data_pool) = c
  return(data_pool)
}

# Aglomerate by Genus/Species
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")

#Relative abundance
Ravel_RA = relat.abun(phyobject = Ravel)
Sriniv_RA = relat.abun(phyobject = Sriniv)

#Prune OTUs
Ravel = prune.OTUs(phyobject = Ravel, Rank = "Rank6")
Ravel_RA = prune.OTUs(phyobject = Ravel_RA, Rank = "Rank6")
Sriniv = prune.OTUs(phyobject = Sriniv, Rank = "Rank6")
Sriniv_RA = prune.OTUs(phyobject = Sriniv_RA, Rank = "Rank6")

#Intersect OTUs
df_S_C = phylo.intersect(Ravel = Ravel, Sriniv = Sriniv, rank = "Rank6")
df_S_RA = phylo.intersect(Ravel = Ravel_RA, Sriniv = Sriniv_RA, rank = "Rank6")

saveRDS(object = df_S_C$Ravel, file = "~/git/BVMetaGenomics/data/GenusIntersect/Ravel_Genus_Counts.rds")
saveRDS(object = df_S_C$Sriniv, file = "~/git/BVMetaGenomics/data/GenusIntersect/Sriniv_Genus_Counts.rds")
saveRDS(object = df_S_RA$Ravel, file = "~/git/BVMetaGenomics/data/GenusIntersect/Ravel_Genus_RA.rds")
saveRDS(object = df_S_RA$Sriniv, file = "~/git/BVMetaGenomics/data/GenusIntersect/Sriniv_Genus_RA.rds")



