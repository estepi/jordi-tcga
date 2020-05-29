library(psichomics)
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/")
dir.create("liver")
folder <- getwd()
cohorts <- 
  getFirebrowseCohorts()
# LIHC Liver hepatocellular carcinoma
# Available sample dates
date <- getFirebrowseDates()
# Available data types
dataTypes <- getFirebrowseDataTypes()
#Junction quantification        "junction_quantification"
# Exon quantification  "exon_quantification"   
#Preprocess  "Preprocess" 
#RSEM isoforms "RSEM_isoforms"
# RSEM genes  "RSEM_genes"  
# Junction expression "junction_expression"
#Gene expression   "gene_expression"
#Exon expression  "exon_expression"
#Genes normalized  "genes_normalized" 
folder <- getwd()
data <- loadFirebrowseData(folder=folder,
                           cohort="LIHC",
                           data=c("clinical"),
                           date="2016-01-28")
data <- loadLocalFiles(folder); head(data)
clinical      <- data[[1]]$`Clinical data`
head(clinical)
dim(clinical)#376 1365
write.table(clinical, "lihc_clinical_psichomics.tab", sep="\t", col.names = NA)
#summary
getwd()