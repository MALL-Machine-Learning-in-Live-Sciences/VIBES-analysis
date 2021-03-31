
#load('TFM/Datos/diabimmune_t1d_metaphlan.rdata')

taxo = as.vector(as.character(mp_shotgun$Taxonomy))
t = strsplit(taxo, '[|]')

reino = rep('Bacteria', length(t))
filo = list()
clase = list()
orden = list()
familia = list()
genero = list()
especie = list()

for (i in 1:length(t)) {
  filo[[i]] = t[[i]][2]
  clase[[i]] = t[[i]][3]
  orden[[i]] = t[[i]][4]
  familia[[i]] = t[[i]][5]
  genero[[i]] = t[[i]][6]
  especie[[i]] = t[[i]][7]
}

df = data.frame(reino = reino,
           filo = unlist(filo),
           clase = unlist(clase),
           orden = unlist(orden),
           familia = unlist(familia),
           genero = unlist(genero),
           especie = unlist(especie))

saveRDS(df, file = 'TFM/Datos//Taxonomy_df.rds')

# Query example!
act = df[which(df$familia == 'f__Actinomycetaceae'),]
act_s = act[which(!is.na(act$especie)),]$especie
act_s = as.vector(as.character(act_s))



nesp = function(dict, familia){
  
  fam = dict[which(dict$familia == familia), ]
  fam = fam[which(!is.na(fam$especie)),]$especie
  nesp = as.vector(as.character(fam))
  
  
  return(length(nesp))
}


nesp(dict = df, 
     familia = 'f__Micrococcaceae')

num = list()
for (i in 1:length(xx)) {
  
  num[[i]] = nesp(dict = df, familia = xx[i])
}

names(num) = xx














