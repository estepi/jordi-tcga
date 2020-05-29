library(survminer)
library(survival)
library(ggfortify)
library(ggplot2)
library(psichomics)
###############cargo las libraries
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/descriptive/")
clinical<-read.table( "/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/psichomics/luad_clinical_psichomics.tab", sep="\t", header = T, row.names = 1)
dim(clinical)#521 2591
cancertype="LUAD"
#select columns
#patient description
gender<-grep("patient.gender", colnames(clinical))
vitalStatus<-grep("patient.vital_status", colnames(clinical))
sample_type<-grep("patient.samples.sample.2.sample_type", colnames(clinical))
age<-grep("patient.age_at_initial_pathologic_diagnosis", colnames(clinical))
#tumor description
stage<-grep("patient.stage_event.pathologic_stage_tumor_stage", colnames(clinical))
#patient survival
follow_up.vital_status<-grep("patient.follow_ups.follow_up.vital_status", colnames(clinical))#alive, dead
days_to_death<-grep("patient.days_to_death", colnames(clinical))
days_last_follow<-grep("patient.days_to_last_followup", colnames(clinical))
days_last_known_alive<-grep("patient.days_to_last_known_alive", colnames(clinical))
fup2<-grep(
  "patient.follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment",
  colnames(clinical))
fup3<-  grep(
  "patient.follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment",
           colnames(clinical))
fup4<-grep(
  "patient.follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment",
    colnames(clinical))
newTumourI<-grep(
  "patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment", 
  colnames(clinical))
newTumourE<-grep(
  "patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment",
  colnames(clinical))
###############################################################
canonical_status<-grep("patient.bcr_canonical_check.bcr_patient_canonical_status", colnames(clinical))
treatment_success<-grep("patient.follow_ups.follow_up.followup_treatment_success", colnames(clinical))
therapy_outcome_success<-grep("patient.follow_ups.follow_up.primary_therapy_outcome_success", colnames(clinical))
primary_therapy_outcome_success<-grep("patient.primary_therapy_outcome_success", colnames(clinical))
###################################################3
#rf survival:
RFS<-data.frame(clinical$patient.days_to_death, 
                clinical$patient.days_to_last_followup, 
                clinical$patient.days_to_last_known_alive,
                clinical$patient.follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)
rownames(RFS)<-rownames(clinical)
head(RFS)
############################################################
#armar subset
subset<-clinical[,unique(c(
                   gender,
                  vitalStatus,
                  sample_type,
                  age,
                  stage,
                  follow_up.vital_status,
                  days_to_death,
                  days_last_follow,
                  days_last_known_alive,
                  fup2,
                  fup3,
                  fup4,
                  newTumourI,
                  newTumourE,
                  canonical_status,
                 treatment_success,
                  therapy_outcome_success,
                  primary_therapy_outcome_success))]
############################################################
head(subset);dim(subset)#521 * 19
colnames(subset)<-sub("patient.", "", colnames(subset))
subset$stage_event.pathologic_stage_tumor_stage<-sub("stage ","",
subset$stage_event.pathologic_stage_tumor_stage);
table(subset$stage_event.pathologic_stage_tumor_stage)
subset$aggStage <- subset$stage
subset$aggStage[subset$aggStage=="i"]<-"early"
subset$aggStage[subset$aggStage=="ia"]<-"early"
subset$aggStage[subset$aggStage=="ib"]<-"early"
subset$aggStage[subset$aggStage=="ii"]<-"early"
subset$aggStage[subset$aggStage=="iia"]<-"early"
subset$aggStage[subset$aggStage=="iib"]<-"early"
subset$aggStage[subset$aggStage=="iii"]<-"late"
subset$aggStage[subset$aggStage=="iiia"]<-"late"
subset$aggStage[subset$aggStage=="iiib"]<-"late"
subset$aggStage[subset$aggStage=="iv"]<-"late"
########################################################
pt<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/tcga/lung_pt.files.txt", header = T, sep="\t")
table(pt$cases.0.disease_type);dim(pt)
convTablePT<-data.frame(case=pt$cases.0.submitter_id, file=pt$file_id)
head(convTablePT);dim(convTablePT)#1044
convTablePT$origin<-rep("pt", nrow(convTablePT))
stn<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/tcga/lung_stn.files.txt", header = T, sep="\t")
head(stn); dim(stn)#110,72
table(stn$cases.0.disease_type)
convTableSTN<-data.frame(case=stn$cases.0.submitter_id, file=stn$file_id)
head(convTableSTN)
convTableSTN$origin<-rep("stn", nrow(convTableSTN))
head(convTableSTN)
##############################################
ConvTotal<-rbind(convTablePT, convTableSTN)
head(ConvTotal)#
table(ConvTotal$origin)
#################################################
#como matcheamos?
psiVal<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/tcga/psi_lung_NA_NewN3.tab", 
                   header = T, 
                   stringsAsFactors = F, 
                   row.names = 1)
psiVal[1:10]; dim(psiVal)#1159
sIds<-colnames(psiVal)[7:ncol(psiVal)]
sIdsB<-gsub('\\.', '-',as.character(sIds)); sIdsB[1:10]
sIdsB<-gsub('^X', '',as.character(sIdsB)); sIdsB[1:10]
colnames(psiVal)<-gsub('\\.', '-',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
colnames(psiVal)<-gsub('^X', '',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
####################################################################################
#check normal PSI
psiNormali<-match(convTableSTN$file, colnames(psiVal))
psiN<-as.numeric(t(psiVal[, psiNormali]))
notna<-which(is.na(psiN))
psiNormal<-as.numeric(psiN[-notna])
summary(psiNormal)
#tal vez no es la maximo sino la media para definir ls grupos
####################################################################################
head(ConvTotal)
table(ConvTotal$origin)#pt 1044, stn 110
length(sIds)#1153
ii<-match(sIdsB, ConvTotal$file)
length(which(is.na(ii)))#0 NAs
#con psichomics:
dim(subset)#503 -> solo PT
head(ConvTotal)
rownames(subset)[1:10]
pp<-match(rownames(subset), ConvTotal$case)
ppi<-pp[!is.na(pp)]
length(ppi)
ppo<-which(!is.na(pp))
length(ppo)
finalDF<-cbind(subset[ppo,], ConvTotal[ppi,])
head(finalDF);dim(finalDF)#516 31
colnames(finalDF)
table(finalDF$origin)#only pt
finalDF<-finalDF[finalDF$origin!="stn",]
dim(finalDF)#516
#agrego el PSI:
dim(psiVal)
finalDF$file[1:10]
colnames(psiVal)[1:10]
vv<-match(finalDF$file, colnames(psiVal))
length(vv)#516
values<-psiVal[,vv]
class(values)
dim(values)
colnames(values)
values[values=="NAold"]<-NA
values[values=="NAnew3"]<-NA
class(values);dim(finalDF);colnames(values)
finalV<-as.numeric(t(values ))
#remove NAs
finalV[1:10]
head(finalV)
names(finalV)<-colnames(values)
#########################################################
df<-data.frame(sample=names(finalV),psi=as.numeric(finalV))
head(df); dim(df)#516, 21 NAs
cc<-match(df$sample, finalDF$file)
finalV<-data.frame(df, finalDF[cc,])
dim(finalV)#516 33
finalV<-finalV[-which(is.na(finalV$psi)),]#494
dim(finalV)#495
finalV$status<-rep(1,  nrow(finalV))
table(finalV$origin)#495 todas pt
######################################################
#aca computamos el optimo sin remover NAs
clinical0<-data.frame(finalV$days_to_last_follow, 
                      finalV$days_to_death,
                      finalV$stage,
                      finalV$gender)
dim(clinical0)#495 4
names(clinical0) <- c("patient.days_to_last_followup", 
                      "patient.days_to_death",
                      "patient.stage_event.pathologic_stage",
                      "patient.gender")
timeStart  <- "days_to_death"
event      <- "days_to_death"
psi<-as.numeric(finalV$psi)
#########################################
opt <- optimalSurvivalCutoff(clinical0, psi, censoring="right", 
                             event="days_to_death", 
                             timeStart="days_to_death")
(optimalCutoff <- opt$par)    # Optimal exon inclusion level 69.40
(optimalPvalue <- opt$value)  # Respective p-value 0.00751
label     <- labelBasedOnCutoff(psi, round(optimalCutoff, 2), 
                                label="PSI values")
survTerms <- processSurvTerms(clinical0, censoring="right",
                              event="days_to_death", timeStart="days_to_death",
                              group=label, scale="years")
surv <- survfit(survTerms)
pvalue <- testSurvival(survTerms)
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
#########################################################
#agrrgo mutaciones:
clinical[1:5,1:5]#rownames are the sample name to match:
#mutaciones:
setwd("../mutations/")
rbm10<-read.delim("rbm10_cbioportal.tsv", sep="\t", header = T)
rbm10[1:5,1:5]
rbm10<-rbm10[-c(1:2),]
dim(rbm10)#21
tcgasamples<-grep("^TCGA", rbm10$Sample.ID)
rbm10S<-rbm10[tcgasamples,]
dim(rbm10S)#19
sampleIDMA<-  
matrix(unlist(strsplit(as.character(rbm10S$Sample.ID), "-")), ncol=4, byrow = T)[,1:3]
head(sampleIDMA)
sID<-apply(sampleIDMA, 1, function(x){paste(x[1],x[2],x[3], sep="-")})
#########################################################
ii<-match(sID, rownames(finalV))
ii
finalV$rbm10Mut<-rep("no", nrow(finalV))
finalV$rbm10Mut[ii]<-"yes"
finalV$color<-rep("skyblue", nrow(finalV))
finalV$color[finalV$rbm10Mut=="yes"]<-"skyblue"
finalV$color[finalV$rbm10Mut=="no"]<-"tomato"
##################other source:
#plot PSI:
sID
library(ggplot2)
cancertype<-"LUAD"
png(paste(cancertype, "by_RBM10_mut_violin.png", sep="_"), 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalV, aes(x=rbm10Mut, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.05),
              size=.5, 
              color= finalV$color   ) + 
  stat_summary(fun.y=median, geom="point", 
               size=3, color=c("red", "blue")  )+
  stat_compare_means(label.y = 100, label.x=2)
dev.off()
getwd()