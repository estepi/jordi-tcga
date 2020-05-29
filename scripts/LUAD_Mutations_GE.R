setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/mutations/")
#como matcheamos?
mut<-read.delim("RBM10_LUAD.tsv", header = T, sep="\t")
psiVal<-read.table("../tcga/psi_lung_NA_NewN3.tab", 
                   header = T, 
                   stringsAsFactors = F, 
                   row.names = 1)
psiVal[1:10]; dim(psiVal)#1159
sIds<-colnames(psiVal)[7:ncol(psiVal)]
sIdsB<-gsub('\\.', '-',as.character(sIds)); sIdsB[1:10]
sIdsB<-gsub('^X', '',as.character(sIdsB)); sIdsB[1:10]
colnames(psiVal)<-gsub('\\.', '-',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
colnames(psiVal)<-gsub('^X', '',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
colnames(psiVal)
####################################################################################
pt<-read.table("../tcga/lung_pt.files.txt", header = T, sep="\t")
table(pt$cases.0.disease_type)
convTablePT<-data.frame(case=pt$cases.0.submitter_id, file=pt$file_id)
head(convTablePT)
dim(convTablePT)#1044
convTablePT$origin<-rep("pt", nrow(convTablePT))
####################################################################################
stn<-read.table("../tcga/lung_stn.files.txt", header = T, sep="\t")
head(stn); dim(stn)#110,72
table(stn$cases.0.disease_type)
convTableSTN<-data.frame(case=stn$cases.0.submitter_id, file=stn$file_id)
head(convTableSTN)
convTableSTN$origin<-rep("stn", nrow(convTableSTN))
ConvTotal<-rbind(convTablePT, convTableSTN)
head(ConvTotal)#
table(ConvTotal$origin)
#como matcheamos?
psiVal<-read.table("../../tcga/psi_lung_NA_NewN3.tab", 
                   header = T, 
                   stringsAsFactors = F, 
                   row.names = 1)
ii<-match(colnames(psiVal), ConvTotal$file)
ii
head(ConvTotal)
table(ConvTotal$origin)#pt 1044, stn 110
length(sIds)#1153
ii<-match(sIdsB, ConvTotal$file)
length(which(is.na(ii)))#0 NAs
###################################################
dim(subset)#503 -> solo PT
head(mut)
head(ConvTotal)
colnames(mut)<-gsub("\\.","-", colnames(mut))
colnames(mut)
mut<-read.delim("RBM10_LUAD.tsv", header = T, sep="\t")
head(mut)
dim(mut)#231 samples
dim(ConvTotal)#1154
head(ConvTotal)
mut[1:5,]
sID<-matrix(unlist(strsplit(as.character(mut$Sample.ID), "-")),ncol=4, byrow = T)[,1:3]
sID

head(sID)
paste(sID[1,], sep="-")
##################################################
mut$newsId<-apply(sID[,1:3],  1,  function(x){ paste0(x, collapse = "-")  } )
dd<-match(mut$newsId,convTablePT$case )
##################################################
clinical<-read.table("../psichomics/luad_clinical_psichomics.tab", sep="\t", header = T, row.names = 1)
colnames(clinical)
clinical[1:5,31:40]
#TCGA-05-4382
head(convTablePT)

getwd()
write.table(row.names(clinical),"small_clinical.tab", sep="\t", col.names = NA)
#rename colnames mutaciones
colnames(mut)[1]
mut[1:5,1:5]
colnames(mut)[which(lapply(strsplit(colnames(mut), "-") , length)!=7)]
         