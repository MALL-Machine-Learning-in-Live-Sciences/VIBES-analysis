## Project: Entropy ##
## Work git directory: ~/git/Entropy ##
## Work project directory: ~/project/Entropy ##
path = "projects/Entropy/"
MetaData = read.delim2("~/projects/Entropy/data/task-nugent-score.txt", header = TRUE, sep = "\t")
names(MetaData)[1] <- "SampleID"
names(MetaData)[2] <- "N_Score"

#Tablas sin RefSeq, uasndo los id que pusieron ellos
Otu97 = read.delim2("~/projects/Entropy/data/otutablegg97.txt", header = TRUE, sep = "\t")
Taxa97 = read.delim2("~/projects/Entropy/data/taxatablegg97.txt", header = TRUE, sep = "\t")

#Tablas usando RefSeq
OtuRef <- read.delim2("~/projects/Entropy/data/otutableRefSeq.txt", header = TRUE, sep = "\t")
TaxaRef = read.delim2("~/projects/Entropy/data/taxatableRefSeq.txt", header = TRUE, sep = "\t")

library(rentrez)
library(XML)
#Usando el numero de acceso de Nucleotide
entrez_dbs()
entrez_db_searchable(db = "nucleotide")
res <- entrez_search(db = "nucleotide", term = "(NR_036982.1[ACCN])")
res$ids
res$count
esums <- entrez_summary(db = "nucleotide", id = res$ids)
tax_id <- extract_from_esummary(esums, "taxid")
#Usando el numero de acceso de Bioproject


entrez_db_searchable(db = "taxonomy")
res2 <- entrez_search(db = "taxonomy", term = "(147802[UID])")
esums2 <- entrez_summary(db = "taxonomy", id = res2$ids)
 
entrez_db_searchable(db = "genome")
res2 <- entrez_search(db = "genome", term = "(lactobacillus iners[orgn])")
esums2 <- entrez_summary(db = "genome", id = res2$ids)


#Prueba 
·BLALALAL

