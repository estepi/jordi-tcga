library(survminer)
library(survival)
library(ggfortify)
library(ggplot2)
library(psichomics)
###############cargo las libraries
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/descriptive/")
clinical<-read.table( "/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc_clinical_psichomics.tab", sep="\t", header = T, row.names = 1)
dim(clinical)#376 1365
write.table(colnames(clinical), file="colnames_lihc.tab", col.names = NA)
cancertype="LIHC"
#select columns
#patient description
gender<-grep("patient.gender", colnames(clinical))
vitalStatus<-grep("patient.vital_status", colnames(clinical))
sample_type<-grep("patient.samples.sample.2.sample_type", colnames(clinical))
age<-grep("patient.age_at_initial_pathologic_diagnosis", colnames(clinical))
#phenotype
muts<-grep("mutation", colnames(clinical)); length(muts)#0
######################################################################
#specific
alk<-grep("alk", colnames(clinical),ignore.case = TRUE); length(alk)#0
#new mutations
EGFR<-grep("EGFR",colnames(clinical),ignore.case = TRUE); length(EGFR)#0
KRAS<-grep("KRAS",colnames(clinical),ignore.case = TRUE); length(KRAS)#0
HER2<-grep("HER2",colnames(clinical),ignore.case = TRUE); length(HER2)#no
BRAF<-grep("BRAF",colnames(clinical),ignore.case = TRUE); length(BRAF)#no
DDR2<-grep("DDR2",colnames(clinical),ignore.case = TRUE); length(DDR2)#no
AKT<-grep("AKT",colnames(clinical),ignore.case = TRUE); length(DDR2)#no
FGFR<-grep("FGFR",colnames(clinical),ignore.case = TRUE); length(FGFR)#no
NOTCH<-grep("NOTCH",colnames(clinical),ignore.case = TRUE); length(NOTCH)#no
TP53<-grep("TP53",colnames(clinical),ignore.case = TRUE); length(TP53)#no
CDH10<-grep("CDH10",colnames(clinical),ignore.case = TRUE); length(CDH10)#no
NFE2L2<-grep("NFE2L2",colnames(clinical),ignore.case = TRUE); length(NFE2L2)#no
KMT2C<-grep("KMT2C",colnames(clinical),ignore.case = TRUE); length(KMT2C)#no
CDKN2A<-grep("CDKN2A",colnames(clinical),ignore.case = TRUE); length(CDKN2A)#no
FAT1<-grep("FAT1",colnames(clinical),ignore.case = TRUE); length(FAT1)#no
PI3KCA<-grep("PI3KCA",colnames(clinical),ignore.case = TRUE); length(PI3KCA)#no
######################################################################
#otras_
mut<-read.table("../lihc/IntOGen-DriverGenes_HC.tsv", header = T, stringsAsFactors = F)
head(mut)
dim(mut)#67

for (i in 1:nrow(mut))
{
  #print(mut$Symbol[i])
  gene<-mut$Symbol[i]
  if( length(grep(gene, colnames(clinical),ignore.case = TRUE ) )!=0 )
     {         print(gene)    }
}
#alb albumin
#atm treatment
######################################################################
trans<-grep("translocation", colnames(clinical)); length(trans)#0
#tumor description
stage<-grep("patient.stage_event.pathologic_stage_tumor_stage", colnames(clinical))
neoplasm_cancer_status<-grep("patient.person_neoplasm_cancer_status", colnames(clinical))
#patient survival
follow_up.vital_status<-grep("patient.follow_ups.follow_up.vital_status", colnames(clinical))#alive, dead
colnames(clinical)[319]
#
days_to_death<-grep("patient.days_to_death", colnames(clinical))#ok
days_last_follow<-grep("patient.days_to_last_followup", colnames(clinical))#ok
colnames(clinical)[114]
######################################################################
######################################################################
days_last_known_alive<-grep("patient.days_to_last_known_alive", colnames(clinical))
#no existe
#podria cambiarse por:
clinical$patient.follow_ups.follow_up.2.days_to_last_followup

fup2<-grep(
  "patient.follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment",
  colnames(clinical))
#no existe
fup3<-  grep(
  "patient.follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment",
           colnames(clinical))
#no existe
fup4<-grep(
  "patient.follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment",
    colnames(clinical))
#no existe
newTumourI<-grep(
  "patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment", 
  colnames(clinical))
#no existe
newTumourE<-grep(
  "patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment",
  colnames(clinical))
#SI
###############################################################
canonical_status<-grep("patient.bcr_canonical_check.bcr_patient_canonical_status", colnames(clinical))
#SI
treatment_success<-grep("patient.follow_ups.follow_up.followup_treatment_success", colnames(clinical))
#NO
therapy_outcome_success<-grep("patient.follow_ups.follow_up.primary_therapy_outcome_success", colnames(clinical))
#NO
primary_therapy_outcome_success<-grep("patient.primary_therapy_outcome_success", colnames(clinical))
#NO
###################################################

clinical$patient.follow_ups.follow_up.days_to_death
RFS<-data.frame(clinical$patient.days_to_death, 
                clinical$patient.days_to_last_followup, 
                clinical$patient.follow_ups.follow_up.2.days_to_death,
                clinical$patient.follow_ups.follow_up.2.days_to_last_followup,
                clinical$patient.follow_ups.follow_up.2.new_tumor_events.new_tumor_event.2.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.2.new_tumor_events.new_tumor_event.days_to_new_tumor_event_additional_surgery_procedure,
                clinical$patient.follow_ups.follow_up.2.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)
rownames(RFS)<-rownames(clinical)
################################################
clinical$patient.follow_ups.follow_up.3.vital_status#alive
clinical$patient.follow_ups.follow_up.4.vital_status#alive
############################################################
#armar subset
muts,EGFR,KRAS,alktrans,days_last_known_alive,fup2,fup3,fup4,newTumourI,treatment_success
therapy_outcome_success,                  primary_therapy_outcome_success

subset<-clinical[,unique(c(gender,
                  vitalStatus,
                  sample_type,
                  age,
                  stage,
                  neoplasm_cancer_status,
                  follow_up.vital_status,
                  days_to_death,
                  days_last_follow,
                  newTumourE,
                  canonical_status))]
############################################################
head(subset);dim(subset)#576 * 12
colnames(subset)<-sub("patient.", "", colnames(subset))
subset$stage_event.pathologic_stage_tumor_stage<-sub("stage ","",
subset$stage_event.pathologic_stage_tumor_stage);
table(subset$stage_event.pathologic_stage_tumor_stage)
subset$aggStage <- subset$stage
subset$aggStage[subset$aggStage=="i"]<-"early"
subset$aggStage[subset$aggStage=="ii"]<-"early"
subset$aggStage[subset$aggStage=="iii"]<-"late"
subset$aggStage[subset$aggStage=="iiia"]<-"late"
subset$aggStage[subset$aggStage=="iiib"]<-"late"
subset$aggStage[subset$aggStage=="iiic"]<-"late"
subset$aggStage[subset$aggStage=="iv"]<-"late"
subset$aggStage[subset$aggStage=="iva"]<-"late"
subset$aggStage[subset$aggStage=="ivb"]<-"late"
########################################################
pt<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc/liver_pt.files.txt", header = T, sep="\t")
table(pt$cases.0.disease_type);dim(pt)
convTablePT<-data.frame(case=pt$cases.0.submitter_id, file=pt$file_id)
colnames(pt)
dim(pt)#86 columnas
head(convTablePT);dim(convTablePT)#371,3
convTablePT$origin<-rep("pt", nrow(convTablePT))
stn<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc/liver_stn.files.txt", header = T, sep="\t")
head(stn); dim(stn)#50,74
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
psiVal<-read.table("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/lihc/lihc_tcga_numb9_notNAs.tab", 
                   header = T, 
                   stringsAsFactors = F, 
                   row.names = 1)
psiVal[1:10]; dim(psiVal)#426
sIds<-colnames(psiVal)[7:ncol(psiVal)]
sIdsB<-gsub('\\.', '-',as.character(sIds)); sIdsB[1:10]
sIdsB<-gsub('^X', '',as.character(sIdsB)); sIdsB[1:10]
colnames(psiVal)<-gsub('\\.', '-',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
colnames(psiVal)<-gsub('^X', '',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
####################################################################################
#check normal PSI
psiNormali<-match(convTableSTN$file, colnames(psiVal))
psiN<-as.numeric(t(psiVal[, psiNormali]))
psiN
notna<-which(is.na(psiN))
psiNormal<-as.numeric(psiN[-notna])
summary(psiNormal)
g0M<-max(psiNormal)
g0M #46.38 #usamos este summary para definir el max de G0
#################################################
####################################################################################
head(ConvTotal)
table(ConvTotal$origin)#pt 371, stn 50
length(sIds)#420
ii<-match(sIdsB, ConvTotal$file)
length(which(is.na(ii)))#0 NAs
#con psichomics:
dim(subset)#376 -> solo PT
head(ConvTotal)
rownames(subset)[1:10]
pp<-match(rownames(subset), ConvTotal$case)
ppi<-pp[!is.na(pp)]
length(ppi)#370
ppo<-which(!is.na(pp))
length(ppo)#370
finalDF<-cbind(subset[ppo,], ConvTotal[ppi,])
head(finalDF);dim(finalDF)#370 16
colnames(finalDF)
table(finalDF$origin)#370 only pt
#agrego el PSI:
dim(psiVal)
finalDF$file[1:10]
colnames(psiVal)[1:10]
vv<-match(finalDF$file, colnames(psiVal))
length(vv)#370
values<-psiVal[,vv]
class(values)
dim(values)
colnames(values)
values[values=="NAold"]<-NA
values[values=="NAnew3"]<-NA
class(values);dim(finalDF);colnames(values)
finalV<-as.numeric(t(values ))
finalV[1:10]
head(finalV)
names(finalV)<-  colnames(values)
#########################################################
df<-data.frame(sample=names(finalV),psi=as.numeric(finalV))
head(df); dim(df)
cc<-match(df$sample, finalDF$file)
length(which(is.na(df$psi)))#84
finalV<-data.frame(df, finalDF[cc,]);dim(finalV)
finalV<-finalV[-which(is.na(finalV$psi)),];dim(finalV)#286
summary(finalV$psi)
dim(finalV)
370-286#84 se pierden por NA
#########################################################
finalV$status<-rep(1,  nrow(finalV))
finalV$status[is.na(finalV$clinical.patient.days_to_death)]<-0
finalV$status[finalV$days_to_death==0]<-0
colnames(finalV)<-sub("_event.pathologic_stage_tumor_stage", "", colnames(finalV))
#################################################################
#plots
finalV$stage
isna<-which(is.na(finalV$stage))
finalVNotNA<-finalV[-isna,]
dim(finalVNotNA)#265
##############################
png(paste(cancertype, "by_stager.png", sep="_"), 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalVNotNA, aes(x=stage, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") 
dev.off()
##################################################
png(paste(cancertype, "by_stage_gender_bp.png", sep="_"), 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalVNotNA, aes(x=stage, y=psi, fill=gender))+ 
  geom_boxplot()   + 
  theme_bw()
dev.off()
##################################################
#by Agg Stage
png(paste(cancertype, "by_aggStage.png", sep="_"), 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalVNotNA, aes(x=aggStage, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") 
dev.off()
##################################################
png(paste(cancertype, "by_aggStage_gender_bp.png", sep="_"), 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalVNotNA, aes(x=aggStage, y=psi, fill=gender))+ 
  geom_boxplot()   + 
  theme_bw()
dev.off()
##################################################
colnames(finalV)<-gsub('-', '',as.character(colnames(finalV) ) );
colnames(finalV)
dim(finalV)#286 22
finalV<-finalV[!is.na(finalV$psi),]
dim(finalV)#286
###################################################
#MUTATIONS
#########################################################################
#ALK
#treatment success
###########################################################
#primary_therapy_outcome_success
################################################################
finalV
finalV$person_neoplasm_cancer_status<-as.character(finalV$person_neoplasm_cancer_status)#SI
finalV$person_neoplasm_cancer_status[is.na(finalV$person_neoplasm_cancer_status)]<-"unknown"
table(finalV$person_neoplasm_cancer_status)
finalVNotUnknown<-finalV[finalV$person_neoplasm_cancer_status!="unknown",]
dim(finalVNotUnknown)#264
################################################################
getwd()
png(paste(cancertype,
          "by_neoplasm_cancer_status.png", sep="_"),
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalVNotUnknown, 
       aes(x=person_neoplasm_cancer_status, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
################################################################
###########################################
###########################################
finalV$samples.sample.2.sample_type<-as.character(finalV$samples.sample.2.sample_type)
finalV$samples.sample.2.sample_type[is.na(finalV$samples.sample.2.sample_type)]<-"unknown"
table(finalV$samples.sample.2.sample_type)
png(paste(cancertype, "_by_sample_type.png", sep="_"),
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalV, 
       aes(x=samples.sample.2.sample_type, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
#########################################################
finalV$age<-as.factor(finalV$age)
agena<-which(is.na(finalV$age))
finalVNotNA<-finalV[-agena,]
#########################################################
png(paste(cancertype, "_by_age.png", sep="_"),
    width = 4, height = 4, units = 'in', res = 300)
ggplot(finalVNotNA, 
       aes(x=age, 
           y=psi)) + 
           geom_point()+
  theme_bw()
dev.off()
#########################################################
#survival plots:
survival<-finalV
survival$days_to_death
#1st step ##########################################################
survival$days_to_death[is.na(survival$days_to_death)]<-0
survival$days_to_last_followup[is.na(survival$days_to_last_followup)]<-0
#2nd step ##########################################################
#3rd step ##########################################################
survival$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment[is.na(survival$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)]<-0
#OS ##########################################################
survival$days==survival$daysOS
survival$days==survival$daysRFS


survival$daysOS<-
  unlist(apply(data.frame(survival$days_to_death,
                          survival$days_to_last_followup),
               1,max));
survival$daysOS
#RFS ##########################################################
survival$daysRFS<-
  unlist(apply(data.frame(survival$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment),
    1,max))
###################################################
table(survival$status)
survival$status[survival$days_to_death==0]<-0
table(survival$status)#221 65
survival$days<-survival$daysRFS
survival$days[survival$daysRFS==0]<-survival$daysOS[survival$daysRFS==0]
#########################################################
#optimo
clinical0<-data.frame(finalV$days_to_last_follow, 
                      finalV$days_to_death,
                      finalV$stage,
                      finalV$gender)
dim(clinical0)
head(clinical0)
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
(optimalCutoff <- opt$par)    # Optimal exon inclusion level 34.39
(optimalPvalue <- opt$value)  # Respective p-value
label     <- labelBasedOnCutoff(psi, round(optimalCutoff, 2), 
                                label="PSI values")
survTerms <- processSurvTerms(clinical0, censoring="right",
                              event="days_to_death", timeStart="days_to_death",
                              group=label, scale="years")
surv <- survfit(survTerms)
pvalue <- testSurvival(survTerms)
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
#########################################################
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi>34.39]<-"g1"
#########################################################

#########################################################
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("Optimo PSI OS")
png(paste(cancertype, "optimo_psi_os.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","black"),
           legend = "bottom",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
#############################################################RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("Optimo PSI RFS")
png(paste(cancertype, "optimo_psi_rfs.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","black"),
           legend = "bottom",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
############################################################
#removemos outliers:

survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi>60]<-"g1"
survivalOutlier<-survival[survival$group=="g1",]
dim(survivalOutlier)#106
boxplot(finalV$days_to_death)
png(paste(cancertype, "group1_days_boxplot.png", sep="_"),
      width = 4, height = 4, units = 'in', res = 300)
  boxplot(survivalOutlier$days_to_death)
  abline(h=c(1000,2000), col="red")
  dev.off()
dim(survivalOutlier[survivalOutlier$days_to_death>1000,])#2
rn<-rownames(survivalOutlier)[survivalOutlier$days_to_death>1000]
survival[rn, ]
ii<-match(rn, rownames(survival))
survival[ii, ]
survivalWOOutl<-survival[-ii,]
dim(survival)#496
dim(survivalWOOutl)#494
#########################################################
#repito survival optimo psi
ii<-match(rn, rownames(finalV))
finalVWO<-finalV[-ii,]
#optimo
clinical0<-data.frame(finalVWO$days_to_last_follow, 
                      finalVWO$days_to_death,
                      finalVWO$stage,
                      finalVWO$gender)
dim(clinical0)
names(clinical0) <- c("patient.days_to_last_followup", 
                      "patient.days_to_death",
                      "patient.stage_event.pathologic_stage",
                      "patient.gender")
timeStart  <- "days_to_death"
event      <- "days_to_death"
psi<-as.numeric(survivalWOOutl$psi)
opt <- optimalSurvivalCutoff(clinical0, psi, censoring="right", 
                             event="days_to_death", 
                             timeStart="days_to_death")
(optimalCutoff <- opt$par)    # Optimal exon inclusion level 50.82
(optimalPvalue <- opt$value)  # Respective p-value
label     <- labelBasedOnCutoff(psi, round(optimalCutoff, 2), 
                                label="PSI values")
survTerms <- processSurvTerms(clinical0, censoring="right",
                              event="days_to_death", timeStart="days_to_death",
                              group=label, scale="years")
surv <- survfit(survTerms)
pvalue <- testSurvival(survTerms)
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
#########################################################
############################################################
dim(survivalOutlier[survivalOutlier$days_to_death>1000,])#6
rn<-rownames(survivalOutlier)[survivalOutlier$days_to_death>1000]
survival[rn, ]
ii<-match(rn, rownames(survival))
survival[ii, ]
survivalWOOutl<-survival[-ii,]
dim(survival)#494
dim(survivalWOOutl)#488
#########################################################
#repito survival optimo psi
ii<-match(rn, rownames(finalV))
finalVWO<-finalV[-ii,]
#optimo
clinical0<-data.frame(finalVWO$days_to_last_follow, 
                      finalVWO$days_to_death,
                      finalVWO$stage,
                      finalVWO$gender)
dim(clinical0)
names(clinical0) <- c("patient.days_to_last_followup", 
                      "patient.days_to_death",
                      "patient.stage_event.pathologic_stage",
                      "patient.gender")
timeStart  <- "days_to_death"
event      <- "days_to_death"
psi<-as.numeric(survivalWOOutl$psi)
opt <- optimalSurvivalCutoff(clinical0, psi, censoring="right", 
                             event="days_to_death", 
                             timeStart="days_to_death")
(optimalCutoff <- opt$par)    # Optimal exon inclusion level 50.82
(optimalPvalue <- opt$value)  # Respective p-value
label     <- labelBasedOnCutoff(psi, round(optimalCutoff, 2), 
                                label="PSI values")
survTerms <- processSurvTerms(clinical0, censoring="right",
                              event="days_to_death", timeStart="days_to_death",
                              group=label, scale="years")
surv <- survfit(survTerms)
pvalue <- testSurvival(survTerms)
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
############################################################
#########################################################
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi<20]<-"g0"
survival$group[survival$psi>20&survival$psi<50]<-"g1*"
table(survival$group)
survival$group[survival$psi>50&survival$psi<68 ]<-"g2*"
survival$group[survival$psi>68]<-"g3*"
sum(table(survival$group))#490
table(survival$group)
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("All groups OS")
png("All_groups_OS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","red","blue","black"),
           legend = "bottom",
           legend.labs = c("g0","g1","g2","g3"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("All groups RFS")
png("All_groups_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","red","blue","black"),
           legend = "bottom",
           legend.labs = c("g0","g1","g2","g3"),
           title=title)
dev.off()
############################################################
#define groups
survival$group[survival$psi<20]<-"g0"
survival$group[survival$psi>20]<-"g1"
table(survival$group)
85+405
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("G1 vs !G1 OS")
png("G1_OS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","red"),
           legend = "bottom",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("G1 vs !G1  RFS")
png("G1_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","red"),
           legend = "bottom",
           legend.labs = c("g0","g1"),
           title=title)
dev.off()
############################################################
survival$group[survival$psi<50]<-"g0"
survival$group[survival$psi>50]<-"g2"
table(survival$group)
395+95
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("G2 vs !G2 OS")
png("G2_OS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","blue"),
           legend = "bottom",
           legend.labs = c("g0","g2"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("G2 vs !G2  RFS")
png("G2_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","blue"),
           legend = "bottom",
           legend.labs = c("g0","g2"),
           title=title)
dev.off()
############################################################
survival$group[survival$psi<68]<-"g0"
survival$group[survival$psi>68]<-"g3"
table(survival$group)
476+14
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("G3 vs !G3 OS")
png("G3_OS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","black"),
           legend = "bottom",
           legend.labs = c("g0","g3"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("G3 vs !G3  RFS")
png("G3_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("pink","black"),
           legend = "bottom",
           legend.labs = c("g0","g3"),
           title=title)
dev.off()
############################################################
#other G1
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi>20&survival$psi<50 ]<-"g1"
table(survival$group)
85+81+14+310
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("G1* vs !G1* OS")
png("G1*_OS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("grey","red"),
           legend = "bottom",
           legend.labs = c("g0*","g1*"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("G1* vs !G1*  RFS")
png("G1*_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("grey","red"),
           legend = "bottom",
           legend.labs = c("g0*","g1*"),
           title=title)
dev.off()
############################################################
#other G2
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi>50&survival$psi<68 ]<-"g2"
table(survival$group)
409+81
85+14+310
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("G2* vs !G2* OS")
png("G2*_OS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("grey","blue"),
           legend = "bottom",
           legend.labs = c("g0*","g2*"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("G2* vs !G2*  RFS")
png("G2*_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survival,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("grey","blue"),
           legend = "bottom",
           legend.labs = c("g0*","g2*"),
           title=title)
dev.off()
############################################################
#pairwise comparison
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi<20]<-"g0"
survival$group[survival$psi>20&survival$psi<50]<-"g1*"
survival$group[survival$psi>50&survival$psi<68 ]<-"g2*"
survival$group[survival$psi>68]<-"g0"
table(survival$group)
sum(table(survival$group))
dim(survival)#490
survivalS<-survival[survival$group!="g0",]
dim(survivalS)#392= 258+134
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalS)
title<-paste("G1* vs G2* OS")
png("G1*_vs_G2*_OS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","blue"),
           legend = "bottom",
           legend.labs = c("g1*","g2*"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survivalS)
title<-paste("G1* vs G2*  RFS")
png("G1*_vs_G2*_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","blue"),
           legend = "bottom",
           legend.labs = c("g1*","g2*"),
           title=title)
dev.off()
############################################################
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi<20]<-"g0"
survival$group[survival$psi>20&survival$psi<50]<-"g0"
survival$group[survival$psi>50&survival$psi<68 ]<-"g2*"
survival$group[survival$psi>68]<-"g3*"
table(survival$group)
sum(table(survival$group))
dim(survival)#490
survivalS<-survival[survival$group!="g0",]
dim(survivalS)#95= 81+14
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalS)
title<-paste("G2* vs G3* OS")
png("G2*_vs_G3*_OS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("blue","grey"),
           legend = "bottom",
           legend.labs = c("g2*","g3*"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survivalS)
title<-paste("G2* vs G3*  RFS")
png("G2*_vs_G3*_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("blue","grey"),
           legend = "bottom",
           legend.labs = c("g2*","g3*"),
           title=title)
dev.off()
########################################################################################################################
#pairwise comparison
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi<20]<-"g0"
survival$group[survival$psi>20&survival$psi<50]<-"g1*"
survival$group[survival$psi>50&survival$psi<68 ]<-"g2*"
survival$group[survival$psi>68]<-"g0"
table(survival$group)
sum(table(survival$group))
dim(survival)#496
survivalS<-survival[survival$group!="g0",]
dim(survivalS)#391= 310+81
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalS)
title<-paste("G1* vs G2* OS")
png("G1*_vs_G2*_OS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","blue"),
           legend = "bottom",
           legend.labs = c("g1*","g2*"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survivalS)
title<-paste("G1* vs G2*  RFS")
png("G1*_vs_G2*_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","blue"),
           legend = "bottom",
           legend.labs = c("g1*","g2*"),
           title=title)
dev.off()
############################################################
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi<20]<-"g0"
survival$group[survival$psi>20&survival$psi<50]<-"g1*"
survival$group[survival$psi>50&survival$psi<68 ]<-"g0"
survival$group[survival$psi>68]<-"g3*"
table(survival$group)
sum(table(survival$group))
dim(survival)#490
survivalS<-survival[survival$group!="g0",]
dim(survivalS)#324= 310 +14
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survivalS)
title<-paste("G1* vs G3* OS")
png("G1*_vs_G3*_OS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","grey"),
           legend = "bottom",
           legend.labs = c("g1*","g3*"),
           title=title)
dev.off()
############################################################
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survivalS)
title<-paste("G1* vs G3*  RFS")
png("G1*_vs_G3*_RFS_LUSC.png")
ggsurvplot(fit, 
           data = survivalS,
           conf.int = F,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           palette = c("red","grey"),
           legend = "bottom",
           legend.labs = c("g1*","g3*"),
           title=title)
dev.off()
############################################################
