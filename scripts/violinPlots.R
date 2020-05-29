library(ggplot2)
#violin plots for outliers:
getwd()
setwd("../")
dir.create("descriptive")
setwd("/media/emancini/e323fb45-bfcd-44b1-9dc8-cde403de4b5a/jordi/descriptive/")
clinicalTable<-read.table("subsetClinical_psi.tab", header = T, row.names = 1)
dim(clinicalTable)#495,21
colnames(clinicalTable)
clinicalTable[1:5,1:5]
table(clinicalTable$stage_event.pathologic_stage_tumor_stage)
table(clinicalTable$aggStage)
################################################
clinicalTable$aggStage[is.na(clinicalTable$aggStage)]<-"late"
table(clinicalTable$aggStage)
png("LUAD_by_stage.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, aes(x=aggStage, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") 
dev.off()
################################################
table(clinicalTable$gender)
png("LUAD_by_gender.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, aes(x=gender, y=psi))+ 
  geom_violin()   + 
  theme_bw()+
  geom_jitter(shape=16, 
              position=position_jitter(0.2),
              size=.5, color="darkgrey")     + 
  stat_summary(fun.y=median, geom="point", size=1, color="red")  + 
  stat_compare_means(label.y = 100, label.x=2)+
  theme(legend.position="none") 
dev.off()
################################################
clinicalTable$follow_ups.follow_up.followup_treatment_success<-as.character(clinicalTable$follow_ups.follow_up.followup_treatment_success)
clinicalTable$follow_ups.follow_up.followup_treatment_success[is.na(clinicalTable$follow_ups.follow_up.followup_treatment_success)]<-"unknown"
table(clinicalTable$follow_ups.follow_up.followup_treatment_success)
png("LUAD_by_followup_treatment_success.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, aes(x=follow_ups.follow_up.followup_treatment_success, y=psi))+ 
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
###############################################
clinicalTable$follow_ups.follow_up.primary_therapy_outcome_success<-as.character(clinicalTable$follow_ups.follow_up.primary_therapy_outcome_success)
clinicalTable$follow_ups.follow_up.primary_therapy_outcome_success[is.na(clinicalTable$follow_ups.follow_up.primary_therapy_outcome_success)]<-"unknown"
table(clinicalTable$follow_ups.follow_up.primary_therapy_outcome_success)
png("LUAD_by_primary_therapy_outcome_success.png" , 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
       aes(x=follow_ups.follow_up.primary_therapy_outcome_success, y=psi))+ 
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
clinicalTable$person_neoplasm_cancer_status<-as.character(clinicalTable$person_neoplasm_cancer_status)
clinicalTable$person_neoplasm_cancer_status[is.na(clinicalTable$person_neoplasm_cancer_status)]<-"unknown"
table(clinicalTable$person_neoplasm_cancer_status)
png("LUAD_by_neoplasm_cancer_status.png", 
    width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
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
clinicalTable$primary_therapy_outcome_success<-as.character(clinicalTable$primary_therapy_outcome_success)
clinicalTable$primary_therapy_outcome_success[is.na(clinicalTable$primary_therapy_outcome_success)]<-"unknown"
table(clinicalTable$primary_therapy_outcome_success)
png("LUAD_by_primary_therapy_outcome_success.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
       aes(x=primary_therapy_outcome_success, y=psi))+ 
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
###########################################
clinicalTable$samples.sample.2.sample_type<-as.character(clinicalTable$samples.sample.2.sample_type)
clinicalTable$samples.sample.2.sample_type[is.na(clinicalTable$samples.sample.2.sample_type)]<-"unknown"
table(clinicalTable$samples.sample.2.sample_type)
png("LUAD_by_sample_type.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
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
colnames(clinicalTable)
clinicalTable$samples.sample.2.sample_type<-as.character(clinicalTable$samples.sample.2.sample_type)
clinicalTable$samples.sample.2.sample_type[is.na(clinicalTable$samples.sample.2.sample_type)]<-"unknown"
table(clinicalTable$samples.sample.2.sample_type)
png("LUAD_by_sample_type.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
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
################################################################
clinicalTable$age<-as.numeric(clinicalTable$age)
clinicalTable$age
png("LUAD_by_age.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
       aes(x=age, 
           y=psi,
           color=gender)) + geom_point()
dev.off()
#############################################################3
clinical[1:2,alk]
clinicalTable$eml4_alk_translocation_method<-as.character(clinicalTable$eml4_alk_translocation_method)
clinicalTable$eml4_alk_translocation_method[is.na(clinicalTable$eml4_alk_translocation_method)]<-"unknown"
table(clinicalTable$eml4_alk_translocation_method)
#############################################################
clinical[1:2,alk]
clinicalTable$eml4_alk_translocation_performed<-as.character(clinicalTable$eml4_alk_translocation_performed)
clinicalTable$eml4_alk_translocation_performed[is.na(clinicalTable$eml4_alk_translocation_performed)]<-"unknown"
table(clinicalTable$eml4_alk_translocation_performed)
#############################################################
png("LUAD_by_alk.png", width = 4, height = 4, units = 'in', res = 300)
ggplot(clinicalTable, 
       aes(x=eml4_alk_translocation_performed, y=psi))+ 
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
