library(survminer)
library(survival)
library(ggfortify)
library(ggplot2)

setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/psichomics/lusc")
clinical<-read.table("luas_clinical_psichomics.tab", sep="\t", header = T, row.names = 1)
head(clinical); rownames(clinical); dim(clinical)
table(clinical$patient.gender)
#select columns
gender<-grep("patient.gender", colnames(clinical))
muts<-grep("mutation", colnames(clinical))
alk<-grep("alk", colnames(clinical))
trans<-grep("translocation", colnames(clinical))
days_to_death<-grep("patient.days_to_death", colnames(clinical))
days_last_follow<-grep("patient.days_to_last_followup", colnames(clinical))
days_last_known_alive<-grep("patient.days_to_last_known_alive", colnames(clinical))
days2 <-grep("patient.follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment", colnames(clinical))
days3 <-grep("patient.follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment", colnames(clinical))
days4 <-grep("patient.follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment", colnames(clinical))
newTumourI<-grep("patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment", colnames(clinical))
newTumourE<-grep("patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment", colnames(clinical))
stage<-grep("patient.stage_event.pathologic_stage_tumor_stage", colnames(clinical))
vitalStatus<-grep("patient.vital_status", colnames(clinical))
neoplasm_cancer_status<-grep("patient.person_neoplasm_cancer_status", colnames(clinical))
follow_up.vital_status<-grep("patient.follow_ups.follow_up.vital_status", colnames(clinical))
canonical_status<-grep("patient.bcr_canonical_check.bcr_patient_canonical_status", colnames(clinical))
treatment_success<-grep("patient.follow_ups.follow_up.followup_treatment_success", colnames(clinical))
therapy_outcome_success<-grep("patient.follow_ups.follow_up.primary_therapy_outcome_success", colnames(clinical))
primary_therapy_outcome_success<-grep("patient.primary_therapy_outcome_success", colnames(clinical))
sample_type<-grep("patient.samples.sample.2.sample_type", colnames(clinical))
age<-grep("patient.age_at_initial_pathologic_diagnosis", colnames(clinical))

############################################################
fup2<-grep(
  "patient.follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment",
  colnames(clinical))
fup3<-grep("patient.follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment",
           colnames(clinical))
fup4<-grep("patient.follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment",
           colnames(clinical))
fup5<-grep("patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment",
           colnames(clinical))
fup6<-grep("patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment", 
           colnames(clinical))
###################################################3
#overall survival:
RFS<-data.frame(clinical$patient.days_to_death, 
                clinical$patient.days_to_last_followup, 
                clinical$patient.days_to_last_known_alive,
                clinical$patient.follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
                clinical$patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)
RFS[1:5,1:5]
rownames(RFS)<-rownames(clinical)
############################################################
subset<-clinical[,c(stage,gender,muts,alk,trans,
                    days_to_death,
                    days_last_follow,
                    days_to_last_known_alive,
                    vitalStatus,
                    neoplasm_cancer_status,
                    follow_up.vital_status,
                                        canonical_status,
                                        treatment_success,
                                        therapy_outcome_success,
                                        primary_therapy_outcome_success,
                                        sample_type,
                    fup2,fup3,fup4,fup5,fup6)]
rownames(subset)
#enlarge using OVDays and RFS
hist(subset$patient.days_to_death)
clinical[1:5,1:5]
rownames(clinical)[1:10]
dim(subset)#503 * 27
head(subset)
subset$patient.stage_event.pathologic_stage_tumor_stage<-sub("stage ","",
                                                             subset$patient.stage_event.pathologic_stage_tumor_stage);
table(subset$patient.stage_event.pathologic_stage_tumor_stage)
#########################################################
colnames(subset)<-sub("patient.", "", colnames(subset))
subset$aggStage<-subset$stage
subset$aggStage[subset$aggStage=="i"]<-"early"
subset$aggStage[subset$aggStage=="ia"]<-"early"
subset$aggStage[subset$aggStage=="ib"]<-"early"
subset$aggStage[subset$aggStage=="ii"]<-"early"
subset$aggStage[subset$aggStage=="iia"]<-"early"
subset$aggStage[subset$aggStage=="iib"]<-"early"
subset$aggStage[subset$aggStage=="iiia"]<-"late"
subset$aggStage[subset$aggStage=="iiib"]<-"late"
subset$aggStage[subset$aggStage=="iv"]<-"late"
########################################################
pt<-read.table("../../tcga/lung_pt.files.txt", header = T, sep="\t")
table(pt$cases.0.disease_type)
convTablePT<-data.frame(case=pt$cases.0.submitter_id, file=pt$file_id)
head(convTablePT)
dim(convTablePT)#1044
convTablePT$origin<-rep("pt", nrow(convTablePT))
stn<-read.table("../../tcga/lung_stn.files.txt", header = T, sep="\t")
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
psiVal[1:10]; dim(psiVal)#1159
sIds<-colnames(psiVal)[7:ncol(psiVal)]
sIdsB<-gsub('\\.', '-',as.character(sIds)); sIdsB[1:10]
sIdsB<-gsub('^X', '',as.character(sIdsB)); sIdsB[1:10]
colnames(psiVal)<-gsub('\\.', '-',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
colnames(psiVal)<-gsub('^X', '',as.character(colnames(psiVal))); colnames(psiVal)[1:10]
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
dim(finalDF)#500 19
head(finalDF)
colnames(finalDF)
table(finalDF$origin)#only pt
#agrego el PSI:
dim(psiVal)
finalDF$file[1:10]
colnames(psiVal)[1:10]
vv<-match(finalDF$file, colnames(psiVal))
length(vv)#500
values<-psiVal[,vv]
values[1:10]
class(values)
dim(values)
colnames(values)
values[values=="NAold"]<-NA
values[values=="NAnew3"]<-NA
class(values);dim(finalDF);colnames(values)
finalV<-as.numeric(t(values ))
finalV[1:10]
head(finalV)
names(finalV)<-colnames(values)
#########################################################
df<-data.frame(sample=names(finalV),psi=as.numeric(finalV))
head(df); dim(df)
cc<-match(df$sample, finalDF$file)
finalV<-data.frame(df, finalDF[cc,])
dim(finalV)
finalV<-finalV[-which(is.na(finalV$psi)),]#494
dim(finalV)
finalV$status<-rep(1,  nrow(finalV))
finalV$status[is.na(finalV$clinical.patient.days_to_death)]<-0
#################################################################
colnames(finalV)[3]<-"stage"
colnames(finalV)[4]<-"gender"
table(finalV$aggStage)
finalV$status<-rep(1,  nrow(finalV))
finalV$status[is.na(finalV$days_to_death)]<-0
finalV$status[finalV$days_to_death==0]<-0
#################################################################
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/descriptive/")
pdf("LUSC_by_stage.pdf", width = 4, height = 4)
ggplot(finalV, aes(x=stage, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") 
dev.off()
png("boxplotPSI_LUSC.png")
par(mfrow=c(1,2))
boxplot(psiEarly)
boxplot(psiLate)
dev.off()
t.test(psiEarly, psiLate)# NS
wilcox.test(psiEarly, psiLate)#NS
################################################
colnames(finalV)<-gsub('patient.', '',as.character(colnames(finalV) ) );
colnames(finalV)<-gsub('-', '',as.character(colnames(finalV) ) );
colnames(finalV)
#remove NAs
dim(finalV)#494 22
finalV<-finalV[!is.na(finalV$psi),]
dim(finalV)#494
finalV$egfr_mutation_performed[is.na(finalV$egfr_mutation_performed)]<-"no"
table(finalV$egfr_mutation_performed)##473 no, 21 yes
finalV$egfr_mutation_result<-as.character(finalV$egfr_mutation_result)
finalV$egfr_mutation_result[is.na(finalV$egfr_mutation_result)]<-"ns"
finalVE<-finalV[finalV$egfr_mutation_performed=="yes",]
table(finalVE$egfr_mutation_result)
table(finalV$egfr_mutation_performed)
table(finalVE$egfr_mutation_performed)
dim(finalVE)#21 -> yes
finalVE$gender_stage<-paste(finalVE$gender, finalVE$aggStage, sep=".")
head(finalVE$gender_stage)
table(finalVE$gender)#5 female, 16 male
table(finalVE$gender_stage)
finalVE$egfr_mutation_result_agg<-as.character(finalVE$egfr_mutation_result)
finalVE$egfr_mutation_result_agg[which(!is.na(finalVE$egfr_mutation_result_agg))]<-"Y"
finalVE$egfr_mutation_result_agg[which(is.na(finalVE$egfr_mutation_result_agg))]<-"N"
table(finalVE$egfr_mutation_performed)
table(finalVE$egfr_mutation_result)
table(finalVE$egfr_mutation_result_agg)
table(finalV$kras_mutation_found); 
finalV[1:10,]
finalKras<-finalV[finalV$kras_mutation_found=="yes",]
dim(finalKras)
levels(finalKras$kras_mutation_found)
table(finalKras$kras_mutation_found)# 0 21 # solos los yes
levels(finalKras$kras_mutation_result)
table(finalKras$kras_mutation_result)
class(finalV$egfr_mutation_performed)
table(finalV$egfr_mutation_performed)
table(finalV$egfr_mutation_result)
row.names(finalV)
head(finalV[finalV$egfr_mutation_performed=="yes",])
egfr<-finalV[finalV$egfr_mutation_performed=="yes",]
dim(egfr)#21,22
table(egfr$egfr_mutation_result)
levels(egfr$egfr_mutation_result)
finalV[1:5,]
finalV$kras_mutation_found[is.na(finalV$kras_mutation_found)]<-"no"
finalV$kras_mutation_result<-as.character(finalV$kras_mutation_result)
finalV$kras_mutation_result[is.na(finalV$kras_mutation_result)]<-"ns"
class(finalV$kras_mutation_found)
table(finalV$kras_mutation_found)
head(finalV[finalV$egfr_mutation_performed=="yes",])
kras<-finalV[finalV$kras_mutation_found=="yes",]
###########################################################
#survival plots:
survival<-finalV
survival$days_to_death
#1st step ##########################################################
survival$days_to_death[is.na(survival$days_to_death)]<-0
survival$days_to_last_followup[is.na(survival$days_to_last_followup)]<-0
survival$days_to_last_known_alive[is.na(survival$days_to_last_known_alive)]<-0
#2nd step ##########################################################
survival$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment[is.na(survival$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment)]<-0
survival$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment[is.na(survival$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment)]<-0
survival$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment[is.na(survival$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment)]<-0
survival$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment[is.na(survival$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment)]<-0
#3rd step ##########################################################
survival$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment[is.na(survival$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)]<-0
#OS ##########################################################
survival$daysOS<-
  unlist(apply(data.frame(survival$days_to_death,
                          survival$days_to_last_followup,
                          survival$days_to_last_known_alive),
               1,max));survival$daysOS
#RFS ##########################################################
survival$daysRFS<-
  unlist(apply(data.frame(
    survival$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
    survival$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
    survival$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
    survival$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
    survival$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment),
    1,max))
###################################################
dim(survival)#494 25
length(which(is.na(survival$aggStage)))#4
notStage<-which(is.na(survival$aggStage))
survival<-survival[-notStage,];dim(survival)#490
survival$combined<-paste(survival$psiV, survival$aggStage)
#######################################################
table(survival$status)
survival$status[survival$days_to_death==0]<-0
table(survival$status)#339 151
survival$psiV<-rep("low", nrow(survival))
survival$psiV[survival$psi>65]<-"high"
survival$days<-survival$daysRFS
survival$days[survival$daysRFS==0]<-survival$daysOS[survival$daysRFS==0]
survivalEarly<-survival[survival$aggStage=="early",]
survivalLate<-survival[survival$aggStage=="late",]
#########################################################
#optimo
survival$group<-rep("g0", nrow(survival))
survival$group[survival$psi>50]<-"g1"
#########################################################
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survival)
title<-paste("Optimo PSI OS")
png("optimo_psi_LUSC.png")
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
#RFS
fit <- survfit(Surv(days/30, status) ~ group,  data = survival)
title<-paste("Optimo PSI RFS")
png("Optimo_psi_RFS_LUSC.png")
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
getwd()
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
