install.packages("kableExtra")

library(ggplot2)
library(survival)
library(ggfortify)
library(knitr)

install.packages("ggpubr")#no
install.packages("maxstat")#no
install.packages("mvtnorm")#no
install.packages("survminer")#no
devtools::install_github("wilkelab/cowplot")#no

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("psichomics")
library(psichomics) #no
###############cargo las libraries
clinical<-read.table( "../../data/luad_clinical_psichomics.tab", sep="\t", header = T, row.names = 1)
dim(clinical)#521 2591
cancertype="LUAD"
#select columns
#patient description
gender<-grep("patient.gender", colnames(clinical))
vitalStatus<-grep("patient.vital_status", colnames(clinical))
sample_type<-grep("patient.samples.sample.2.sample_type", colnames(clinical))
age<-grep("patient.age_at_initial_pathologic_diagnosis", colnames(clinical))
#phenotype
muts<-grep("mutation", colnames(clinical)); length(muts)#4
colnames(clinical)[muts]
######################################################################
#specific
alk<-grep("alk", colnames(clinical),
          ignore.case = TRUE); length(alk)#2
#new mutations
EGFR<-grep("EGFR",colnames(clinical),ignore.case = TRUE); 
length(EGFR)#2
KRAS<-grep("KRAS",colnames(clinical),ignore.case = TRUE); length(KRAS)#3
HER2<-grep("HER2",colnames(clinical),ignore.case = TRUE); 
length(HER2)#no
BRAF<-grep("BRAF",colnames(clinical),ignore.case = TRUE);
length(BRAF)#no
DDR2<-grep("DDR2",colnames(clinical),ignore.case = TRUE); 
length(DDR2)#no
AKT<-grep("AKT",colnames(clinical),ignore.case = TRUE); 
length(DDR2)#no
FGFR<-grep("FGFR",colnames(clinical),ignore.case = TRUE);
length(FGFR)#no
NOTCH<-grep("NOTCH",colnames(clinical),ignore.case = TRUE);
length(NOTCH)#no
TP53<-grep("TP53",colnames(clinical),ignore.case = TRUE); 
length(TP53)#no
CDH10<-grep("CDH10",colnames(clinical),ignore.case = TRUE); 
length(CDH10)#no
NFE2L2<-grep("NFE2L2",colnames(clinical),ignore.case = TRUE); length(NFE2L2)#no
KMT2C<-grep("KMT2C",colnames(clinical),ignore.case = TRUE); length(KMT2C)#no
CDKN2A<-grep("CDKN2A",colnames(clinical),ignore.case = TRUE); length(CDKN2A)#no
FAT1<-grep("FAT1",colnames(clinical),ignore.case = TRUE); length(FAT1)#no
PI3KCA<-grep("PI3KCA",colnames(clinical),ignore.case = TRUE); length(PI3KCA)#no
######################################################################
trans<-grep("translocation", colnames(clinical)); length(trans)#2
#tumor description
stage<-grep("patient.stage_event.pathologic_stage_tumor_stage", colnames(clinical))
neoplasm_cancer_status<-grep("patient.person_neoplasm_cancer_status", colnames(clinical))
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
############################################################
#armar subset
subset<-clinical[,unique(c(gender,
                  vitalStatus,
                  sample_type,
                  age,
                  muts,
                  EGFR,#2
                  KRAS,
                  alk,
                  trans,
                  stage,
                  neoplasm_cancer_status,
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
head(subset);dim(subset)#521 * 27
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
########################################################
pt<-read.table("../../data//lung_pt.files.txt", header = T, sep="\t")
table(pt$cases.0.disease_type);dim(pt)
convTablePT<-data.frame(case=pt$cases.0.submitter_id, file=pt$file_id)
head(convTablePT);dim(convTablePT)#1044
convTablePT$origin<-rep("pt", nrow(convTablePT))
stn<-read.table("../../data/lung_stn.files.txt", header = T, sep="\t")
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
psiVal<-read.table("../../data/psi_lung_NA_NewN3.tab", 
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
#ahora reploteo por grupos con colores adhoc, etc
#cambia el pvalue porque ahora el OS se redefine como:
#survival plots:
#
survival<-finalV
#OS ##########################################################
NAs<-
data.frame(survival$days_to_death,
                          survival$days_to_last_followup,
                         survival$days_to_last_known_alive);head(NAs)
SUMNAS<-rowSums(apply(NAs, 2,is.na) ); length(SUMNAS); table(SUMNAS)#19 muestras que son triple NA
#
#saco las triple NAs ->
survivalNotNas<-survival[SUMNAS!=3,]
dim(survival)#495
dim(survivalNotNas)#476
495-476#19
#
survivalNotNas
survivalNotNas$days_to_death[is.na(survivalNotNas$days_to_death)]<-0
survivalNotNas$days_to_last_followup[is.na(survivalNotNas$days_to_last_followup)]<-0
survivalNotNas$days_to_last_known_alive[is.na(survivalNotNas$days_to_last_known_alive)]<-0
#
survivalNotNas$daysOS<-
  unlist(apply(data.frame(survivalNotNas$days_to_death,
                          survivalNotNas$days_to_last_followup,
                          survivalNotNas$days_to_last_known_alive),
               1,max));
survivalNotNas$daysOS
#
#RFS ##########################################################
daysRFSDF<-
  data.frame(
    survivalNotNas$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)
###################################################
SUMNASRFS<-rowSums(apply(daysRFSDF, 2,is.na) );
table(SUMNASRFS)
#ahora corrijo con 0
survivalNotNas$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment[
  is.na(survivalNotNas$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment)]<-0

survivalNotNas$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment[
  is.na(survivalNotNas$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment)]<-0

survivalNotNas$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment[
  is.na(survivalNotNas$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment)]<-0

survivalNotNas$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment[
  is.na(survivalNotNas$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment)]<-0

survivalNotNas$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment[
  is.na(survivalNotNas$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)]<-0

survivalNotNas$daysRFS<-
  unlist(apply(data.frame(
    survivalNotNas$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment),
    1,max))
#############################################
dim(survivalNotNas)#476
head(df)
df<-data.frame(survivalNotNas$daysRFS,SUMNASRFS) 
df[df$SUMNASRFS==1,]
dim(daysRFSDF)#476 muestras 5 columnas, busco aquellas donde NAs=5 -> 307
#si hay 0 es porque habia quintuple NAS
#como ya corregi por OS, remplaz como dice el texto
###################################################
table(survivalNotNas$status)#1 475
survivalNotNas$status[survivalNotNas$days_to_death==0]<-0
#339 151 antes
table(survivalNotNas$status)
#359 117 ahora
survivalNotNas$days<-survivalNotNas$daysRFS
survivalNotNas$days[survivalNotNas$daysRFS==0]<-survivalNotNas$daysOS[survivalNotNas$daysRFS==0]
#########################################################
#optimo
survivalNotNas$group<-rep("g0", nrow(survivalNotNas))
survivalNotNas$group[survivalNotNas$psi>69.41]<-"g1"
dim(survivalNotNas)#495
table(survivalNotNas$group)
dim(survivalNotNas)
#########################################################
cancertype="LUAD"
getwd()
dir.create("days_analysis2")
setwd("days_analysis2/")
#OS 
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalNotNas)
sum(is.infinite(survivalNotNas$daysOS))
survivalNotNas$daysOS
title<-paste("Optimo 69.41 PSI OS, Remuevo triple NAs en OS")
png(paste(cancertype, "optimo_psi_os.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalNotNas,
           conf.int = F,
           pval = FALSE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
######################################################
title<-paste("Optimo 69.41 PSI OS, Remuevo triple NAs en OS")
png(paste(cancertype, "optimo_psi_os_pval.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalNotNas,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#########################################################
#RFS:
fit <- survfit(Surv(days/30, status) ~ group,  data = survivalNotNas)
survivalNotNas$days
#########################################################
title<-paste("Optimo 69.41 PSI RFS, Remuevo triple NAs en OS")
png(paste(cancertype, "optimo_psi_rfs.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalNotNas,
           conf.int = F,
           pval = FALSE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
######################################################
title<-paste("Optimo 69.41 PSI RFS, Remuevo triple NAs en OS")
png(paste(cancertype, "optimo_psi_rfs_pval.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalNotNas,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#########################################################
#optimize en g1:
clinicalTumor<-data.frame(days_to_last_follow=finalV$days_to_last_follow, 
                          days_to_death=finalV$days_to_death,
                          stage=finalV$stage,
                          gender=finalV$gender,
                          psi=finalV$psi)
head(clinicalTumor)
clinicalTumor<-clinicalTumor[clinicalTumor$psi>=69.41,]
dim(clinicalTumor)#47 (1 tiene triple NA)
names(clinicalTumor) <- c("patient.days_to_last_followup", 
                          "patient.days_to_death",
                          "patient.stage_event.pathologic_stage",
                          "patient.gender",
                          "psi")
timeStart  <- "days_to_death"
event      <- "days_to_death"
psi<-as.numeric(clinicalTumor$psi)
#########################################
opt <- optimalSurvivalCutoff(clinicalTumor,
                             psi, censoring="right", 
                             event="days_to_death", 
                             timeStart="days_to_death")
(optimalCutoff <- opt$par)    # Optimal exon inclusion level 69.76
(optimalPvalue <- opt$value)  # Respective p-value 0.00751
label     <- labelBasedOnCutoff(psi, round(optimalCutoff, 2), 
                                label="PSI values")
survTerms <- processSurvTerms(clinicalTumor, censoring="right",
                              event="days_to_death", timeStart="days_to_death",
                              group=label, scale="months")
surv <- survfit(survTerms)
pvalue <- testSurvival(survTerms)
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
#########################################################
table(survival$group)#29 466
table(survivalNotNas$group)#46 son G1
survivalHigh<-survivalNotNas[survivalNotNas$group=="g1",]
dim(survivalHigh)#46
survivalHigh$group<-rep("g0", nrow(survivalHigh))
survivalHigh$group[survivalHigh$psi>83.41]<-"g1"
dim(survivalHigh)#46
table(survivalHigh$group)#41 5
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalHigh)
#########################################################
title<-paste("PSI 83.41  (optimal psi inside g1 high) OS")
png(paste(cancertype, "psi_os_high.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalHigh,
           conf.int = F,
           pval = FALSE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#################################################################
title<-paste("PSI 83.41  (optimal psi inside g1 high) OS")
png(paste(cancertype, "psi_os_high_pval.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalHigh,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#################################################################
fit <- survfit(Surv(days/30, status) ~ group,  data = survivalHigh)
#########################################################
title<-paste("PSI 83.41  (optimal psi inside g1 high) RFS")
png(paste(cancertype, "psi_rfs_high.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalHigh,
           conf.int = F,
           pval = FALSE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#################################################################
title<-paste("PSI 83.41  (optimal psi inside g1 high) RFS")
png(paste(cancertype, "psi_rfs_high_pval.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalHigh,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("skyblue","tomato"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#################################################################


#########################################################
#extra plots optimize psi only in PSI tumor:
#optimo
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi>16.450]<-"g1"
dim(survival)#495
table(survival$group)#29 466
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
survival$daysOS[is.infinite(survival$daysOS)]<-0
survival$daysOS
setwd("days_analysis/")
title<-paste("PSI 16.41 (tumor/normal) OS")
png(paste(cancertype, "psi_os_tumor_normal.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = FALSE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("dodgerblue","chocolate1"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
######################################################
png(paste(cancertype, "psi_os_tumor_normal_pvalT.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("dodgerblue","chocolate1"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
######################################################
#Optimize inside tumor:
#aca computamos el optimo sin remover NAs
clinicalTumor<-data.frame(finalV$days_to_last_follow, 
                      finalV$days_to_death,
                      finalV$stage,
                      finalV$gender,
                      finalV$psi)
clinicalTumor<-clinicalTumor[clinicalTumor$psi>16.45,]
dim(clinicalTumor)#466 4
names(clinicalTumor) <- c("patient.days_to_last_followup", 
                      "patient.days_to_death",
                      "patient.stage_event.pathologic_stage",
                      "patient.gender",
                      "psi")
timeStart  <- "days_to_death"
event      <- "days_to_death"
psi<-as.numeric(clinicalTumor$psi)
#########################################
opt <- optimalSurvivalCutoff(clinicalTumor,
                             psi, censoring="right", 
                             event="days_to_death", 
                             timeStart="days_to_death")
(optimalCutoff <- opt$par)    # Optimal exon inclusion level 69.76
(optimalPvalue <- opt$value)  # Respective p-value 0.00751
label     <- labelBasedOnCutoff(psi, round(optimalCutoff, 2), 
                                label="PSI values")
survTerms <- processSurvTerms(clinicalTumor, censoring="right",
                              event="days_to_death", timeStart="days_to_death",
                              group=label, scale="months")
surv <- survfit(survTerms)
pvalue <- testSurvival(survTerms)
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
#########################################################
table(survival$group)#29 466
survivalTumor<-survival[survival$group=="g1",]
dim(survivalTumor)#466
survivalTumor$group<-rep("g0", nrow(survivalTumor))
survivalTumor$group[survivalTumor$psi>69.77]<-"g1"
dim(survivalTumor)#466
table(survivalTumor$group)#420 46
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalTumor)
survivalTumor$daysOS[is.infinite(survivalTumor$daysOS)]<-0
survivalTumor$daysOS
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalTumor)
setwd("days_analysis/")
title<-paste("PSI 69.77  (optimal psi inside tumor) OS")
png(paste(cancertype, "psi_os_only_tumor.png", sep="_"),
    width = 6, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survivalTumor,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("dodgerblue","chocolate1"),
           legend = "right",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()

