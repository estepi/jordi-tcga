library(survminer)
library(survival)
library(ggfortify)
library(ggplot2)
library(psichomics)
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a1/jordi/psichomics/lusc")
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
days_to_last_known_alive
days_last_known_alive
subset<-clinical[,c(stage,gender,muts,alk,trans,
                    days_to_death,
                    days_last_follow,
                    days_last_known_alive,
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
dir.create("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/descriptive2/")
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/descriptive2/")
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
#cuidado: la estimacion del PSI optimo es sin remover NAs:
table(finalV$origin)#494 todas pt
######################################################
#aca computamos el optimo sin remover NAs
clinical0<-data.frame(finalV$days_to_last_follow, 
                      finalV$days_to_death,
                      finalV$stage,
                      finalV$gender)
dim(clinical0)#494 4
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
                              group=label, scale="months")
surv <- survfit(survTerms)#50.83
pvalue <- testSurvival(survTerms)#0.37
plotSurvivalCurves(surv, pvalue=pvalue, mark=FALSE)
#########################################################
#ahora corrijo los NAs
survival<-finalV
dim(finalV)#494
dim(clinical0)#494
dim(survival)#494
#OS ##########################################################
NAs<-
  data.frame(survival$days_to_death,
             survival$days_to_last_followup,
             survival$days_to_last_known_alive);head(NAs)
SUMNAS<-rowSums(apply(NAs, 2,is.na) ); length(SUMNAS); table(SUMNAS)
#11 muestras que son triple NA
#
#saco las triple NAs ->
survivalNotNas<-survival[SUMNAS!=3,]
dim(survival)#483
dim(survivalNotNas)#494-483#11
#1st step OS ##########################################################
survivalNotNas$days_to_death[is.na(survivalNotNas$days_to_death)]<-0
survivalNotNas$days_to_last_followup[is.na(survivalNotNas$days_to_last_followup)]<-0
survivalNotNas$days_to_last_known_alive[is.na(survivalNotNas$days_to_last_known_alive)]<-0
##########################################################
survivalNotNas$daysOS<-
  unlist(apply(data.frame(survivalNotNas$days_to_death,
                          survivalNotNas$days_to_last_followup,
                          survivalNotNas$days_to_last_known_alive),
               1,max));survivalNotNas$daysOS
#2nd step RFS ##########################################################
daysRFSDF<-
  data.frame(
    survivalNotNas$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment)
SUMNASRFS<-rowSums(apply(daysRFSDF, 2,is.na) );
table(SUMNASRFS)#355 que tienen quintuple NA -> esos deberian ser remplazados por los dias de OS
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
#######################################################
survivalNotNas$daysRFS<-
  unlist(apply(data.frame(
    survivalNotNas$follow_ups.follow_up.2.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.3.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.4.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment,
    survivalNotNas$new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment),
    1,max))
#######################################################
table(survivalNotNas$daysRFS)#355 que tienen 0 -> van a tener el valor de OS -> todas vienen de los NAs
dim(survivalNotNas)#483
######################################################
#optimo
survPSI<-survivalNotNas
survPSI$group<-rep("g0", nrow(survPSI))
survPSI$group[survPSI$psi>=50.83]<-"g1"
dim(survPSI)#483
table(survPSI$group)
###############################################################################
getwd()
setwd("../../ML/")
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survPSI)
cancertype="LUSC"
title<-paste("Optimo 50.83 PSI OS.")
png(paste(cancertype, "optimo_psi_os.png", sep="_"),
    width = 4, height = 6, units = 'in', res = 300)
ggsurvplot(fit, 
           data = survPSI,
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
#######################################################
#try to do the plot
survPSI<-survivalNotNas
FindOptimalPSI<-function(survPSI, cancertype)
{
  df<-data.frame(matrix(NA, 
                        ncol=2, 
                        nrow= length(seq(min(survPSI$psi),max(survPSI$psi), by=0.05))))
  colnames(df)<-c("psi","pval")
  ###################################################################
  for ( i in 1:length(  seq(min(survPSI$psi),max(survPSI$psi), by=0.05))  )  {
    print(i)
    psi=seq(min(survPSI$psi),max(survPSI$psi), by=0.05)[i]
    df[i,1]<-psi
    print(psi)
    survPSI$group<-rep("g0", nrow(survPSI))
    survPSI$group[survPSI$psi>psi]<-"g1"
    fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survPSI)
    pval<-surv_pvalue(fit, data = NULL, method = "survdiff",
                      test.for.trend = FALSE, combine = FALSE)[2]
    df[i,2]<-pval
  }
  title<-paste(cancertype, "PSI OPTIMO",df[which.min(df$pval),1], "pval", round(df[which.min(df$pval),2], digits = 4))
  file<-paste(paste(cancertype, "PSIOPTIMO",df[which.min(df$pval),1], "pval", round(df[which.min(df$pval),2], digits = 4),sep="_"),
              "png", sep=".")
  png(file, width = 4, height = 4, units = 'in', res = 300)
  ggplot(data=df, 
         aes(x=psi, y=pval, group=1)) + 
    geom_line(colour="grey")+
    theme_bw()+
    geom_hline(yintercept=range(0.05, 0.1), color='coral', size=0.5)+
    ggtitle(title)
  dev.off()
  getwd()
}
###################################################################
cancertype
survPSI<-survivalNotNas
survPSI<-survivalNotNas[survivalNotNas$psi<=18.28,] #g0
survPSI<-survivalNotNas[survivalNotNas$psi>=18.28,] #g1
#############################################################################
#repito los survival plots:
survPSI<-survivalNotNas
#############################################################################
survPSI$group<-rep("g0", nrow(survPSI))#control
plot(sort(survPSI$psi))
survPSI$group[survPSI$psi>=18.28]<-"g1"
dim(survPSI)#483
table(survPSI$group)#68 415
###############################################################################
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survPSI)
title<-paste("Optimo 18.28 PSI OS")
png(paste(cancertype,
          "optimo_psi_os.png", sep="_"),
    width = 4, 
    height = 6, 
    units = 'in', 
    res = 300)
ggsurvplot(fit, 
           data = survPSI,
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
##############################################################################
#split DF according PSI optimo:
survPSI<-survivalNotNas[survivalNotNas$psi<=18.28,]
#############################################################################
survPSI$group<-rep("g0", nrow(survPSI))#control
plot(sort(survPSI$psi))
survPSI$group[survPSI$psi>=11.88]<-"g1"
dim(survPSI)#
table(survPSI$group)#28 48
###############################################################################
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survPSI)
title<-paste("Optimo G0 11.88 PSI OS")
png(paste(cancertype,
          "optimo_psi_G0_os.png", sep="_"),
    width = 4, 
    height = 6, 
    units = 'in', 
    res = 300)
ggsurvplot(fit, 
           data = survPSI,
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

###########################################################
#split DF according PSI optimo:
survPSI<-survivalNotNas[survivalNotNas$psi>=18.28,]
#############################################################################
survPSI$group<-rep("g0", nrow(survPSI))#control
plot(sort(survPSI$psi))
survPSI$group[survPSI$psi>=21.15]<-"g1"
dim(survPSI)#
table(survPSI$group)#30 385
###############################################################################
fit <- survfit(Surv(daysOS/30, status) ~ group,  data = survPSI)
title<-paste("Optimo G1 21.15 PSI OS")
png(paste(cancertype,
          "optimo_psi_G1_os.png", sep="_"),
    width = 4, 
    height = 6, 
    units = 'in', 
    res = 300)
ggsurvplot(fit, 
           data = survPSI,
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



