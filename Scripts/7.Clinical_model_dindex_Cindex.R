
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Clinical_models1.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/clinical_features_included_censored_As_well_new_PROB.RData")

######3 Including  censored data for C-index and D-INDEX caluculation and comprasion across the models.
model1=clinical_models$model1
model3=clinical_models$model3

pcsi=data.frame(clinical_features$pcsi_clinical)
icgc=data.frame(clinical_features$icgc_clinical)
tcga=data.frame(clinical_features$tcga_clinical)
icgc_arr=data.frame(clinical_features$icgc_arr_clinical)
ouh=data.frame(clinical_features$ouh_clinical)

pcsi$Age=as.numeric(as.character(pcsi$Age))
pcsi$Sex = as.character(pcsi$Sex)
pcsi$T_status=gsub(" ", "", pcsi$T_status, fixed = TRUE)
pcsi$pred_prob=as.numeric(as.character(pcsi$pred_prob))

pcsi_cl_pred1=1-predict(model1,pcsi, na.action = na.exclude)
pcsi_cl_pred2= pcsi$pred_prob[which(pcsi$ID %in% names(pcsi_cl_pred1))]
pcsi_cl_pred3=1-predict(model3,pcsi, na.action = na.exclude)

######### Validation cohorts including censored samples as well
########## ICGC

icgc$Age=as.numeric(as.character(icgc$Age))
icgc$Sex = as.character(icgc$Sex)
icgc$T_status = as.character(icgc$T_status)

icgc_cl_pred1=1-predict(model1, icgc, na.action = na.exclude)
icgc_cl_pred2= icgc$pred_prob[which(icgc$ID %in% names(icgc_cl_pred1))]
icgc_cl_pred3=1-predict(model3, icgc, na.action = na.exclude)

########## TCGA

tcga$Age=as.numeric(as.character(tcga$Age))
tcga$Sex = as.character(tcga$Sex)
tcga$T_status = as.character(tcga$T_status)
tcga$pred_prob=as.numeric(as.character(tcga$pred_prob))

tcga_cl_pred1=1-predict(model1, tcga, na.action = na.exclude)
tcga_cl_pred2= tcga$pred_prob[which(tcga$ID %in% names(tcga_cl_pred1))]
tcga_cl_pred3=1-predict(model3, tcga, na.action = na.exclude)

########## ICGC array

icgc_arr$Age=as.numeric(as.character(icgc_arr$Age))
icgc_arr$Sex = as.character(icgc_arr$Sex)
icgc_arr$T_status = as.character(icgc_arr$T_status)
icgc_arr$pred_prob=as.numeric(as.character(icgc_arr$pred_prob))

icgc_arr_cl_pred1=1-predict(model1, icgc_arr, na.action = na.exclude)
icgc_arr_cl_pred2= icgc_arr$pred_prob[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))]
icgc_arr_cl_pred3=1-predict(model3, icgc_arr, na.action = na.exclude)

########## OUH

ouh$Age=as.numeric(as.character(ouh$Age))
ouh$Sex = as.character(ouh$Sex)

ouh_cl_pred1 = 1-predict(model1, ouh, na.action = na.exclude)
ouh_cl_pred2= ouh$pred_prob
ouh_cl_pred3=1-predict(model3, ouh, na.action = na.exclude)

###########################################################################
###########################################################################
###########################################################################
### Concordance indices and Dindex calcualtion 

## Clinical model
dindex_ouh <- D.index(x=ouh_cl_pred1, surv.time=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS)), surv.event=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_ouh <- concordance.index(x=ouh_cl_pred1, surv.time=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS)), surv.event=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_pcsi <- D.index(x=pcsi_cl_pred1, surv.time=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS)), surv.event=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_pcsi <- concordance.index(x=pcsi_cl_pred1, surv.time=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS)), surv.event=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_tcga <- D.index(x=tcga_cl_pred1, surv.time=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS)), surv.event=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_tcga <- concordance.index(x=tcga_cl_pred1, surv.time=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS)), surv.event=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_icgc_array <- D.index(x=icgc_arr_cl_pred1 , surv.time=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS)), surv.event=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_icgc_array <- concordance.index(x=icgc_arr_cl_pred1 , surv.time=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS)), surv.event=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");

##PCOSP
dindex_ouh1 <- D.index(x=ouh_cl_pred2, surv.time=as.numeric(as.character(ouh$OS)), surv.event=as.numeric(as.character(ouh$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_ouh1 <- concordance.index(x=ouh_cl_pred2, surv.time=as.numeric(as.character(ouh$OS)), surv.event=as.numeric(as.character(ouh$OS_Status)), na.rm=TRUE, method="noether");
dindex_pcsi1 <- D.index(x=pcsi_cl_pred2, surv.time=as.numeric(as.character(pcsi$OS[which(pcsi$ID %in% names(pcsi_cl_pred1))])), surv.event=as.numeric(as.character(pcsi$OS_Status[which(pcsi$ID %in% names(pcsi_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_pcsi1 <- concordance.index(x=pcsi_cl_pred2, surv.time=as.numeric(as.character(pcsi$OS[which(pcsi$ID %in% names(pcsi_cl_pred1))])), surv.event=as.numeric(as.character(pcsi$OS_Status[which(pcsi$ID %in% names(pcsi_cl_pred1))])), na.rm=TRUE, method="noether");
dindex_tcga1 <- D.index(x=tcga_cl_pred2, surv.time=as.numeric(as.character(tcga$OS[which(tcga$ID %in% names(tcga_cl_pred1))])), surv.event=as.numeric(as.character(tcga$OS_Status[which(tcga$ID %in% names(tcga_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_tcga1 <- concordance.index(x=tcga_cl_pred2, surv.time=as.numeric(as.character(tcga$OS[which(tcga$ID %in% names(tcga_cl_pred1))])), surv.event=as.numeric(as.character(tcga$OS_Status[which(tcga$ID %in% names(tcga_cl_pred1))])), na.rm=TRUE, method="noether");
dindex_icgc_array1 <- D.index(x=icgc_arr_cl_pred2 , surv.time=as.numeric(as.character(icgc_arr$OS[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), surv.event=as.numeric(as.character(icgc_arr$OS_Status[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_icgc_array1 <- concordance.index(x=icgc_arr_cl_pred2 , surv.time=as.numeric(as.character(icgc_arr$OS[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), surv.event=as.numeric(as.character(icgc_arr$OS_Status[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), na.rm=TRUE, method="noether");

## PCOSP+Clinical model
dindex_ouh2 <- D.index(x=ouh_cl_pred3, surv.time=as.numeric(as.character(ouh[names(ouh_cl_pred3),]$OS)), surv.event=as.numeric(as.character(ouh[names(ouh_cl_pred3),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_ouh2 <- concordance.index(x=ouh_cl_pred3, surv.time=as.numeric(as.character(ouh[names(ouh_cl_pred3),]$OS)), surv.event=as.numeric(as.character(ouh[names(ouh_cl_pred3),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_pcsi2 <- D.index(x=pcsi_cl_pred3, surv.time=as.numeric(as.character(pcsi$OS[which(pcsi$ID %in% names(pcsi_cl_pred1))])), surv.event=as.numeric(as.character(pcsi$OS_Status[which(pcsi$ID %in% names(pcsi_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_pcsi2 <- concordance.index(x=pcsi_cl_pred3, surv.time=as.numeric(as.character(pcsi$OS[which(pcsi$ID %in% names(pcsi_cl_pred1))])), surv.event=as.numeric(as.character(pcsi$OS_Status[which(pcsi$ID %in% names(pcsi_cl_pred1))])), na.rm=TRUE, method="noether");
dindex_tcga2 <- D.index(x=tcga_cl_pred3, surv.time=as.numeric(as.character(tcga$OS[which(tcga$ID %in% names(tcga_cl_pred1))])), surv.event=as.numeric(as.character(tcga$OS_Status[which(tcga$ID %in% names(tcga_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_tcga2 <- concordance.index(x=tcga_cl_pred3, surv.time=as.numeric(as.character(tcga$OS[which(tcga$ID %in% names(tcga_cl_pred1))])), surv.event=as.numeric(as.character(tcga$OS_Status[which(tcga$ID %in% names(tcga_cl_pred1))])), na.rm=TRUE, method="noether");
dindex_icgc_array2 <- D.index(x=icgc_arr_cl_pred3 , surv.time=as.numeric(as.character(icgc_arr$OS[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), surv.event=as.numeric(as.character(icgc_arr$OS_Status[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_icgc_array2 <- concordance.index(x=icgc_arr_cl_pred3 , surv.time=as.numeric(as.character(icgc_arr$OS[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), surv.event=as.numeric(as.character(icgc_arr$OS_Status[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), na.rm=TRUE, method="noether");

####################################### Meta estimates calculations
## Clinical model
dindex_meta <- combine.est(c( dindex_pcsi$d.index, dindex_tcga$d.index, dindex_ouh$d.index,dindex_icgc_array$d.index),c(dindex_pcsi$se, dindex_tcga$se, dindex_ouh$se,dindex_icgc_array$se),na.rm = TRUE)
con_meta <- combine.est(c( con_pcsi$c.index, con_tcga$c.index, con_ouh$c.index,con_icgc_array$c.index),c(con_pcsi$se, con_tcga$se, con_ouh$se,con_icgc_array$se),na.rm = TRUE)
dindex_seq <- combine.est(c( dindex_pcsi$d.index, dindex_tcga$d.index),c(dindex_pcsi$se, dindex_tcga$se),na.rm = TRUE)
con_seq <- combine.est(c( con_pcsi$c.index, con_tcga$c.index),c(con_pcsi$se, con_tcga$se),na.rm = TRUE)
dindex_micro <- combine.est(c( dindex_ouh$d.index,dindex_icgc_array$d.index),c( dindex_ouh$se,dindex_icgc_array$se),na.rm = TRUE)
con_micro <- combine.est(c(  con_ouh$c.index,con_icgc_array$c.index),c( con_ouh$se,con_icgc_array$se),na.rm = TRUE)
dindex_meta_lower <- dindex_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta$se
con_meta_lower <- con_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta$se
dindex_seq_lower <- dindex_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq$se
con_seq_lower <- con_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq$se
dindex_micro_lower <- dindex_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro$se
con_micro_lower <- con_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro$se
dindex_meta_upper <- dindex_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta$se
con_meta_upper <- con_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta$se
dindex_seq_upper <- dindex_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq$se
con_seq_upper <- con_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq$se
dindex_micro_upper <- dindex_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro$se
con_micro_upper <- con_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro$se
dindex_meta_pval <- combine.test(p=c(dindex_pcsi$p.value, dindex_tcga$p.value, dindex_ouh$p.value,dindex_icgc_array$p.value),w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1),length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
dindex_micro_pval <- combine.test(p=c(dindex_ouh$p.value,dindex_icgc_array$p.value),w=c(length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
dindex_seq_pval <- combine.test(p=c(dindex_pcsi$p.value, dindex_tcga$p.value),w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1)),hetero = FALSE,method="z.transform")
con_meta_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value, con_ouh$p.value,con_icgc_array$p.value),w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1),length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
con_micro_pval <- combine.test(p=c(con_ouh$p.value,con_icgc_array$p.value),w=c(length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
con_seq_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value),w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1)),hetero = FALSE,method="z.transform")

## TSP ensembl model
dindex_meta1 <- combine.est(c( dindex_pcsi1$d.index, dindex_tcga1$d.index, dindex_ouh1$d.index,dindex_icgc_array1$d.index),c(dindex_pcsi1$se, dindex_tcga1$se, dindex_ouh1$se,dindex_icgc_array1$se),na.rm = TRUE)
con_meta1 <- combine.est(c(con_pcsi1$c.index, con_tcga1$c.index, con_ouh1$c.index, con_icgc_array1$c.index),c(con_pcsi1$se, con_tcga1$se, con_ouh1$se,con_icgc_array1$se),na.rm = TRUE)
dindex_seq1 <- combine.est(c( dindex_pcsi1$d.index, dindex_tcga1$d.index),c(dindex_pcsi1$se, dindex_tcga1$se),na.rm = TRUE)
con_seq1 <- combine.est(c( con_pcsi1$c.index, con_tcga1$c.index),c(con_pcsi1$se, con_tcga1$se),na.rm = TRUE)
dindex_micro1 <- combine.est(c( dindex_ouh1$d.index,dindex_icgc_array1$d.index),c( dindex_ouh1$se,dindex_icgc_array1$se),na.rm = TRUE)
con_micro1 <- combine.est(c(  con_ouh1$c.index,con_icgc_array1$c.index),c( con_ouh1$se,con_icgc_array1$se),na.rm = TRUE)
dindex_meta_lower1 <- dindex_meta1$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta1$se
con_meta_lower1 <- con_meta1$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta1$se
dindex_seq_lower1 <- dindex_seq1$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq1$se
con_seq_lower1 <- con_seq1$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq1$se
dindex_micro_lower1 <- dindex_micro1$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro1$se
con_micro_lower1 <- con_micro1$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro1$se
dindex_meta_upper1 <- dindex_meta1$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta1$se
con_meta_upper1 <- con_meta1$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta1$se
dindex_seq_upper1 <- dindex_seq1$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq1$se
con_seq_upper1 <- con_seq1$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq1$se
dindex_micro_upper1 <- dindex_micro1$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro1$se
con_micro_upper1 <- con_micro1$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro1$se
dindex_meta_pval1 <- combine.test(p=c(dindex_pcsi1$p.value, dindex_tcga1$p.value, dindex_ouh1$p.value,dindex_icgc_array1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2),length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
dindex_micro_pval1 <- combine.test(p=c(dindex_ouh1$p.value,dindex_icgc_array1$p.value),w=c(length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
dindex_seq_pval1 <- combine.test(p=c(dindex_pcsi1$p.value, dindex_tcga1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2)),hetero = FALSE,method="z.transform")
con_meta_pval1 <- combine.test(p=c(con_pcsi1$p.value, con_tcga1$p.value, con_ouh1$p.value,con_icgc_array1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2),length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
con_micro_pval1 <- combine.test(p=c(con_ouh1$p.value,con_icgc_array1$p.value),w=c(length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
con_seq_pval1 <- combine.test(p=c(con_pcsi1$p.value, con_tcga1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2)),hetero = FALSE,method="z.transform")

## Clinical model + TSP ensembl model
dindex_meta2 <- combine.est(c( dindex_pcsi2$d.index, dindex_tcga2$d.index, dindex_ouh2$d.index,dindex_icgc_array2$d.index),c(dindex_pcsi2$se, dindex_tcga2$se, dindex_ouh2$se,dindex_icgc_array2$se),na.rm = TRUE)
con_meta2 <- combine.est(c(con_pcsi2$c.index, con_tcga2$c.index, con_ouh2$c.index, con_icgc_array2$c.index),c(con_pcsi2$se, con_tcga2$se, con_ouh2$se,con_icgc_array2$se),na.rm = TRUE)
dindex_seq2 <- combine.est(c( dindex_pcsi2$d.index, dindex_tcga2$d.index),c(dindex_pcsi2$se, dindex_tcga2$se),na.rm = TRUE)
con_seq2 <- combine.est(c( con_pcsi2$c.index, con_tcga2$c.index),c(con_pcsi2$se, con_tcga2$se),na.rm = TRUE)
dindex_micro2 <- combine.est(c( dindex_ouh2$d.index,dindex_icgc_array2$d.index),c( dindex_ouh2$se,dindex_icgc_array2$se),na.rm = TRUE)
con_micro2 <- combine.est(c(  con_ouh2$c.index,con_icgc_array2$c.index),c( con_ouh2$se,con_icgc_array2$se),na.rm = TRUE)
dindex_meta_lower2 <- dindex_meta2$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta2$se
con_meta_lower2 <- con_meta2$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta2$se
dindex_seq_lower2 <- dindex_seq2$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq2$se
con_seq_lower2 <- con_seq2$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq2$se
dindex_micro_lower2 <- dindex_micro2$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro2$se
con_micro_lower2 <- con_micro2$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro2$se
dindex_meta_upper2 <- dindex_meta2$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta2$se
con_meta_upper2 <- con_meta2$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta2$se
dindex_seq_upper2 <- dindex_seq2$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq2$se
con_seq_upper2 <- con_seq2$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq2$se
dindex_micro_upper2 <- dindex_micro2$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro2$se

con_micro_upper2 <- con_micro2$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro2$se
dindex_meta_pval2 <- combine.test(p=c(dindex_pcsi2$p.value, dindex_tcga2$p.value, dindex_ouh2$p.value,dindex_icgc_array2$p.value),w=c(length(pcsi_cl_pred3), length(tcga_cl_pred3),length(ouh_cl_pred3),length(icgc_arr_cl_pred3)),hetero = FALSE,method="z.transform")
dindex_micro_pval2 <- combine.test(p=c(dindex_ouh2$p.value,dindex_icgc_array2$p.value),w=c(length(ouh_cl_pred3),length(icgc_arr_cl_pred3)),hetero = FALSE,method="z.transform")
dindex_seq_pval2 <- combine.test(p=c(dindex_pcsi2$p.value, dindex_tcga2$p.value),w=c(length(pcsi_cl_pred3), length(tcga_cl_pred3)),hetero = FALSE,method="z.transform")
con_meta_pval2 <- combine.test(p=c(con_pcsi2$p.value, con_tcga2$p.value, con_ouh2$p.value,con_icgc_array2$p.value),w=c(length(pcsi_cl_pred3), length(tcga_cl_pred3),length(ouh_cl_pred3),length(icgc_arr_cl_pred3)),hetero = FALSE,method="z.transform")
con_micro_pval2 <- combine.test(p=c(con_ouh2$p.value,con_icgc_array2$p.value),w=c(length(ouh_cl_pred3),length(icgc_arr_cl_pred3)),hetero = FALSE,method="z.transform")
con_seq_pval2 <- combine.test(p=c(con_pcsi2$p.value, con_tcga2$p.value),w=c(length(pcsi_cl_pred3), length(tcga_cl_pred3)),hetero = FALSE,method="z.transform")


############# Plotting Concordance index and comparison across models
#pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Manuscript-Figures/Clinical-c-index.pdf")
r.mean <- c( con_tcga$c.index,con_pcsi$c.index, con_pcsi$c.index,con_ouh$c.index,con_seq$estimate,con_micro$estimate,con_meta$estimate,NA, NA,con_tcga1$c.index,con_pcsi1$c.index, con_icgc_array1$c.index, con_ouh1$c.index, con_seq1$estimate, con_micro1$estimate, con_meta1$estimate, NA, NA,con_tcga2$c.index,con_pcsi2$c.index, con_icgc_array2$c.index, con_ouh2$c.index,con_seq2$estimate, con_micro2$estimate, con_meta2$estimate)
r.lower <- c( con_tcga$lower,con_pcsi$lower, con_pcsi$lower,con_ouh$lower,con_seq_lower, con_micro_lower, con_meta_lower,NA, NA,con_tcga1$lower,con_pcsi1$lower, con_icgc_array1$lower, con_ouh1$lower, con_seq_lower1, con_micro_lower1, con_meta_lower1, NA, NA,con_tcga2$lower,con_pcsi2$lower, con_icgc_array2$lower, con_ouh2$lower, con_seq_lower2, con_micro_lower2, con_meta_lower2)
r.upper <- c( con_tcga$upper,con_pcsi$upper, con_pcsi$upper,con_ouh$upper,con_seq_upper, con_micro_upper, con_meta_upper, NA,NA, con_tcga1$upper,con_pcsi1$upper, con_icgc_array1$upper, con_ouh1$upper, con_seq_upper1, con_micro_upper1, con_meta_upper1, NA, NA,con_tcga2$upper,con_pcsi2$upper, con_icgc_array2$upper, con_ouh2$upper, con_seq_upper2, con_micro_upper2, con_meta_upper2)
r.pval <- round(c(con_tcga$p.value, con_pcsi$p.value,con_icgc_array$p.value,con_ouh$p.value, con_seq_pval, con_micro_pval, con_meta_pval, NA, NA,con_tcga1$p.value, con_pcsi1$p.value,con_icgc_array1$p.value,con_ouh1$p.value,  con_seq_pval1, con_micro_pval1, con_meta_pval1,NA, NA,con_tcga2$p.value, con_pcsi2$p.value,con_icgc_array2$p.value,con_ouh2$p.value,  con_seq_pval2, con_micro_pval2, con_meta_pval2),2)

r.pval1 <- c(sprintf("%.1E", con_tcga$p.value),sprintf("%.1E", con_pcsi$p.value), sprintf("%.1E", con_icgc_array$p.value),sprintf("%.1E",con_ouh$p.value), sprintf("%.1E",con_seq_pval ),sprintf("%.1E", con_micro_pval), sprintf("%.1E", con_meta_pval),NA,NA,sprintf("%.1E", con_tcga1$p.value),sprintf("%.1E", con_pcsi1$p.value), sprintf("%.1E", con_icgc_array1$p.value),sprintf("%.1E",con_ouh1$p.value), sprintf("%.1E",con_seq_pval1 ),sprintf("%.1E", con_micro_pval1), sprintf("%.1E", con_meta_pval1),NA,NA,sprintf("%.1E", con_tcga2$p.value),sprintf("%.1E", con_pcsi2$p.value), sprintf("%.1E", con_icgc_array2$p.value),sprintf("%.1E",con_ouh2$p.value), sprintf("%.1E",con_seq_pval2),sprintf("%.1E", con_micro_pval2), sprintf("%.1E", con_meta_pval2))

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","PCSI", "ICGC-array", "OUH","Sequencing","Micro-array","Overall","", "PCOSP","TCGA","PCSI", "ICGC-array", "OUH","Sequencing","Micro-array","Overall"," ","Combined PCOSP and clinical", "TCGA","PCSI", "ICGC-array", "OUH","Sequencing","Micro-array","Overall")

data2 <- 
  structure(list(
    mean  = c(NA,NA,t[,1]),
    lower = c(NA,NA,t[,2]),
    upper = c(NA,NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, NA, -26L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts","Clinical model",rownames(t)),
  c("P values","",r.pval1))

seq=length(tcga_cl_pred1) +  length(pcsi_cl_pred1)
arr=length(icgc_arr_cl_pred1)+ length(ouh_cl_pred1)
length_c=c(length(tcga_cl_pred1), length(pcsi_cl_pred1), length(icgc_arr_cl_pred1), length(ouh_cl_pred1), seq,arr,seq+arr,  
           length(tcga_cl_pred1), length(pcsi_cl_pred1), length(icgc_arr_cl_pred1), length(ouh_cl_pred1), seq,arr,seq+arr,  
           length(tcga_cl_pred1), length(pcsi_cl_pred1), length(icgc_arr_cl_pred1), length(ouh_cl_pred1), seq,arr,seq+arr
)/1000

pdf("/Users/vandanasandhu/Desktop/c.pdf")

fn <- local({
  i = 0
  
  b_clrs =  c(c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"))
  l_clrs =   c(c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"))
  #s_clrs =c(rep("red",10),"green","pink","yellow","orange")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    #fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fn1 <- local({
  i = 0
  
  s_clrs =c("#FF7F00","#1F78B4","grey57","#FF7F00","#1F78B4","grey57","#FF7F00","#1F78B4","grey57")
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})


forestplot(tabletext2,data2,xlab="Concordance index",is.summary=c(TRUE,TRUE,rep(FALSE,4),rep(TRUE,3),FALSE,TRUE, rep(FALSE,4),rep(TRUE,3),FALSE,TRUE,rep(FALSE,4),rep(TRUE,3)),clip=c(0,3.0),cex=10,
           fn.ci_norm = fn,  fn.ci_sum = fn1,zero=0.5,graphwidth=unit(2, "inches"), align=c("l"), new_page = FALSE,txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  
                                                                                                                                    xlab  = gpar(fontfamily = "Helvetica", cex = 1)), col = fpColors(text="black"))

############# Plotting D-index and comparison across models
#pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Manuscript-Figures/Clinical-d-index.pdf")
r.mean <- c( log2(dindex_tcga$d.index), log2(dindex_pcsi$d.index), log2(dindex_icgc_array$d.index),log2(dindex_ouh$d.index),log2(dindex_seq$estimate),log2(dindex_micro$estimate),log2(dindex_meta$estimate),NA, NA,log2(dindex_tcga1$d.index),log2(dindex_pcsi1$d.index), log2(dindex_icgc_array1$d.index), log2(dindex_ouh1$d.index), log2(dindex_seq1$estimate), log2(dindex_micro1$estimate), log2(dindex_meta1$estimate), NA, NA,log2(dindex_tcga2$d.index),log2(dindex_pcsi2$d.index), log2(dindex_icgc_array2$d.index), log2(dindex_ouh2$d.index),log2(dindex_seq2$estimate), log2(dindex_micro2$estimate), log2(dindex_meta2$estimate))
r.lower <- c( log2(dindex_tcga$lower),log2(dindex_pcsi$lower), log2(dindex_icgc_array$lower),log2(dindex_ouh$lower),log2(dindex_seq_lower), log2(dindex_micro_lower), log2(dindex_meta_lower),NA, NA,log2(dindex_tcga1$lower),log2(dindex_pcsi1$lower), log2(dindex_icgc_array1$lower), log2(dindex_ouh1$lower), log2(dindex_seq_lower1), log2(dindex_micro_lower1), log2(dindex_meta_lower1), NA, NA,log2(dindex_tcga2$lower),log2(dindex_pcsi2$lower), log2(dindex_icgc_array2$lower), log2(dindex_ouh2$lower), log2(dindex_seq_lower2), log2(dindex_micro_lower2), log2(dindex_meta_lower2))
r.upper <- c( log2(dindex_tcga$upper),log2(dindex_pcsi$upper), log2(dindex_icgc_array$upper),log2(dindex_ouh$upper),log2(dindex_seq_upper), log2(dindex_micro_upper), log2(dindex_meta_upper), NA,NA, log2(dindex_tcga1$upper),log2(dindex_pcsi1$upper), log2(dindex_icgc_array1$upper), log2(dindex_ouh1$upper), log2(dindex_seq_upper1), log2(dindex_micro_upper1), log2(dindex_meta_upper1), NA, NA,log2(dindex_tcga2$upper),log2(dindex_pcsi2$upper), log2(dindex_icgc_array2$upper), log2(dindex_ouh2$upper), log2(dindex_seq_upper2), log2(dindex_micro_upper2), log2(dindex_meta_upper2))
r.pval <- round(c(dindex_tcga$p.value, dindex_pcsi$p.value,dindex_icgc_array$p.value,dindex_ouh$p.value, dindex_seq_pval, dindex_micro_pval, dindex_meta_pval, NA, NA,dindex_tcga1$p.value, dindex_pcsi1$p.value,dindex_icgc_array1$p.value,dindex_ouh1$p.value,  dindex_seq_pval1, dindex_micro_pval1, dindex_meta_pval1,NA, NA,dindex_tcga2$p.value, dindex_pcsi2$p.value,dindex_icgc_array2$p.value,dindex_ouh2$p.value,  dindex_seq_pval2, dindex_micro_pval2, dindex_meta_pval2),2)

r.pval1 <- c(sprintf("%.1E", dindex_tcga$p.value),sprintf("%.1E", dindex_pcsi$p.value), sprintf("%.1E", dindex_icgc_array$p.value),sprintf("%.1E",dindex_ouh$p.value), sprintf("%.1E",dindex_seq_pval ),sprintf("%.1E", dindex_micro_pval), sprintf("%.1E", dindex_meta_pval),NA,NA,sprintf("%.1E", dindex_tcga1$p.value),sprintf("%.1E", dindex_pcsi1$p.value), sprintf("%.1E", dindex_icgc_array1$p.value),sprintf("%.1E",dindex_ouh1$p.value), sprintf("%.1E",dindex_seq_pval1 ),sprintf("%.1E", dindex_micro_pval1), sprintf("%.1E", dindex_meta_pval1),NA,NA,sprintf("%.1E", dindex_tcga2$p.value),sprintf("%.1E", dindex_pcsi2$p.value), sprintf("%.1E", dindex_icgc_array2$p.value),sprintf("%.1E",dindex_ouh2$p.value), sprintf("%.1E",dindex_seq_pval2),sprintf("%.1E", dindex_micro_pval2), sprintf("%.1E", dindex_meta_pval2))



t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","PCSI", "ICGC-array", "OUH","Sequencing","Micro-array","Overall","", "PCOSP","TCGA","PCSI", "ICGC-array", "OUH","Sequencing","Micro-array","Overall"," ","Combined PCOSP and clinical", "TCGA","PCSI", "ICGC-array", "OUH","Sequencing","Micro-array","Overall")

data2 <- 
  structure(list(
    mean  = c(NA,NA,t[,1]),
    lower = c(NA,NA,t[,2]),
    upper = c(NA,NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, NA, -26L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts","Clinical model",rownames(t)),
  c("P values","",r.pval1))
pdf("/Users/vandanasandhu/Desktop/d.pdf")
fn <- local({
  i = 0
  
  b_clrs =  c(c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"))
  l_clrs =   c(c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"),c("#FF7F00","#FF7F00","#1F78B4","#1F78B4"))
  #s_clrs =c(rep("red",10),"green","pink","yellow","orange")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    #fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fn1 <- local({
  i = 0
  
  s_clrs =c("#FF7F00","#1F78B4","grey57","#FF7F00","#1F78B4","grey57","#FF7F00","#1F78B4","grey57")
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

forestplot(tabletext2,data2,xlab="Log2 D-index",is.summary=c(TRUE,TRUE,rep(FALSE,4),rep(TRUE,3),FALSE,TRUE, rep(FALSE,4),rep(TRUE,3),FALSE,TRUE,rep(FALSE,4),rep(TRUE,3)),clip=c(-2,3.0),cex=9,
           fn.ci_norm = fn,  fn.ci_sum = fn1,zero=0,graphwidth=unit(2, "inches"), align=c("l"), new_page = FALSE,txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  
                                                                                                                                  xlab  = gpar(fontfamily = "Helvetica", cex = 1)), col = fpColors(text="black"))

















