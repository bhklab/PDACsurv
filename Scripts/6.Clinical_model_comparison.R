load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/clinical_features_ICGC_prob1.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/cohorts_subtypes_Average_data.RData")
library(verification)
library(reportROC)
library(glmnet)
library(ggthemes)
############################## ###########################################

pcsi=data.frame(clinical_features$pcsi_clinical)
icgc=data.frame(clinical_features$icgc_clinical)
tcga=data.frame(clinical_features$tcga_clinical)
icgc_arr=data.frame(clinical_features$icgc_arr_clinical)
ouh=data.frame(clinical_features$ouh_clinical)

icgc$Age=as.numeric(as.character(icgc$Age))
icgc$Sex = as.character(icgc$Sex)
icgc$binary_grp = as.numeric(as.character(icgc$binary_grp))
icgc$T_status=gsub(" ", "", icgc$T_status, fixed = TRUE)

##### Model building comparisons
model1 <- glm(binary_grp ~ Age + Sex + T_status+ N + M + Grade, data=icgc , na.action = na.exclude,family = binomial(link = "logit"))
model3 <- glm(binary_grp ~ Age + Sex + T_status+ N + M + Grade + pred_prob, data=icgc , na.action = na.exclude,family = binomial(link = "logit"))

summary(model1); 
summary(model3);
anova(model1);
anova(model3);
####### Predicting probabilities in validation cohorts

########## 1. PCSI

pcsi$Age=as.numeric(as.character(pcsi$Age))
pcsi$Sex = as.character(pcsi$Sex)
pcsi$binary_grp = as.numeric(as.character(pcsi$binary_grp))
pcsi$T_status = as.character(pcsi$T_status)


pcsi_cl_pred1=predict(model1, pcsi, na.action = na.exclude,type="response")
pcsi_cl_pred2=1-pcsi$pred_prob[which(pcsi$ID %in% names(pcsi_cl_pred1))]

########## 2. TCGA

tcga$Age=as.numeric(as.character(tcga$Age))
tcga$Sex = as.character(tcga$Sex)
tcga$binary_grp = as.numeric(as.character(tcga$binary_grp))
tcga$T_status = as.character(tcga$T_status)
tcga$pred_prob=as.numeric(as.character(tcga$pred_prob))
#tcga$Grade=as.character(tcga$Grade)

tcga_cl_pred1=predict(model1, tcga, na.action = na.exclude,type="response")
tcga_cl_pred2=1-tcga$pred_prob[which(tcga$ID %in% names(tcga_cl_pred1))]

########## 3. ICGC array

icgc_arr$Age=as.numeric(as.character(icgc_arr$Age))
icgc_arr$Sex = as.character(icgc_arr$Sex)
icgc_arr$binary_grp = as.numeric(as.character(icgc_arr$binary_grp))
icgc_arr$T_status = as.character(icgc_arr$T_status)


icgc_arr1=icgc_arr[order(icgc_arr$ID),]

icgc_arr_cl_pred1 = predict(model1, icgc_arr, na.action = na.exclude,type="response")
icgc_arr_cl_pred2 = 1-icgc_arr$pred_prob[which(icgc_arr1$ID %in% names(icgc_arr_cl_pred1))]

########## 4.OUH

ouh$Age=as.numeric(as.character(ouh$Age))
ouh$Sex = as.character(ouh$Sex)
ouh$binary_grp = as.numeric(as.character(ouh$binary_grp))

ouh_cl_pred1 = predict(model1, ouh, na.action = na.exclude,type="response")
ouh_cl_pred2 = 1-ouh$pred_prob


######## Model AUCS ESTIMATION ####################
###########################################################################
pcsi_roc1=reportROC(pcsi[names(pcsi_cl_pred1),]$binary_grp,pcsi_cl_pred1); pcsi_roc_se1=pcsi_roc1$AUC.SE
tcga_roc1=reportROC(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred1); tcga_roc_se1=tcga_roc1$AUC.SE
ouh_roc1=reportROC(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred1); ouh_roc_se1=ouh_roc1$AUC.SE
icgc_array_roc1 = reportROC(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred1); icgc_array_roc_se1 = icgc_array_roc1$AUC.SE

pcsi_roc2=reportROC(pcsi[names(pcsi_cl_pred1),]$binary_grp,pcsi_cl_pred2); pcsi_roc_se2=pcsi_roc2$AUC.SE
tcga_roc2=reportROC(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred2); tcga_roc_se2=tcga_roc2$AUC.SE
ouh_roc2=reportROC(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred2); ouh_roc_se2=ouh_roc2$AUC.SE
icgc_array_roc2 = reportROC(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred2);icgc_array_roc_se2 = icgc_array_roc2$AUC.SE

pcsi_roc1_pval= roc.area(pcsi[names(pcsi_cl_pred1),]$binary_grp,pcsi_cl_pred1)$p.value  
tcga_roc1_pval= roc.area(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred1)$p.value  
ouh_roc1_pval= roc.area(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred1)$p.value  
icgc_array_roc1_pval= roc.area(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred1)$p.value  

pcsi_roc2_pval= roc.area(pcsi[names(pcsi_cl_pred1),]$binary_grp,pcsi_cl_pred2)$p.value  
tcga_roc2_pval= roc.area(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred2)$p.value  
ouh_roc2_pval= roc.area(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred2)$p.value  
icgc_array_roc2_pval= roc.area(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred2)$p.value  


######## Plotting barplot
df.plot=data.frame(cohort=c("1.TCGA","1.TCGA","2.PCSI","2.PCSI","3.ICGC-array","3.ICGC-array","4.OUH","4.OUH"), 
                   model=c(rep(c("Clinicopathological model","PCOSP"),4)), 
                   AUC=c(tcga_roc1$AUC, tcga_roc2$AUC ,pcsi_roc1$AUC, pcsi_roc2$AUC,icgc_array_roc1$AUC,
                         icgc_array_roc2$AUC,ouh_roc1$AUC,ouh_roc2$AUC),
                   AUC_up      =        c(tcga_roc1$AUC + tcga_roc1$AUC.SE,
                                          tcga_roc2$AUC + tcga_roc2$AUC.SE,
                                          pcsi_roc1$AUC + pcsi_roc1$AUC.SE,
                                          pcsi_roc2$AUC + pcsi_roc2$AUC.SE,
                                          icgc_array_roc1$AUC +   icgc_array_roc1$AUC.SE,
                                          icgc_array_roc2$AUC +   icgc_array_roc2$AUC.SE,
                                          ouh_roc1$AUC +  ouh_roc1$AUC.SE,
                                          ouh_roc2$AUC +  ouh_roc2$AUC.SE),
                   
                   AUC_low      =        c(tcga_roc1$AUC - tcga_roc1$AUC.SE,
                                           tcga_roc2$AUC - tcga_roc2$AUC.SE,
                                           pcsi_roc1$AUC - pcsi_roc1$AUC.SE,
                                           pcsi_roc2$AUC - pcsi_roc2$AUC.SE,
                                           icgc_array_roc1$AUC -   icgc_array_roc1$AUC.SE,
                                           icgc_array_roc2$AUC -   icgc_array_roc2$AUC.SE,
                                           ouh_roc1$AUC -  ouh_roc1$AUC.SE,
                                           ouh_roc2$AUC -  ouh_roc2$AUC.SE))

ggplot(df.plot, aes(x=cohort, y=AUC, fill=model)) +
  geom_bar( stat='identity',width = 0.8, position = position_dodge(width = 0.85)) +
  xlab(NULL) + scale_fill_manual(values=c("palevioletred1", "darkgrey"))+
  geom_errorbar(aes(ymin = AUC_low, ymax = AUC_up),
                width = 0.2,
                linetype = "dotted",
                position = position_dodge(width = 0.9),
                color="black", size=1) +
theme_hc() + theme(legend.position="none",axis.text = element_text(size = 12),axis.title  = element_text(size = 12))


########## Meta-estimate Calculations
model1_meta = combine.est(c(pcsi_roc1$AUC, tcga_roc1$AUC, ouh_roc1$AUC,icgc_array_roc1$AUC), c( pcsi_roc_se1, tcga_roc_se1, ouh_roc_se1,icgc_array_roc_se1),na.rm=TRUE,hetero=TRUE)$estimate
model1_seq = combine.est(c(pcsi_roc1$AUC, tcga_roc1$AUC), c( pcsi_roc_se1, tcga_roc_se1),na.rm=TRUE,hetero=TRUE)$estimate
model1_micro =combine.est(c( ouh_roc1$AUC,icgc_array_roc1$AUC), c( ouh_roc_se1,icgc_array_roc_se1),na.rm=FALSE,hetero=TRUE)$estimate

model2_meta = combine.est(c(pcsi_roc2$AUC, tcga_roc2$AUC, ouh_roc2$AUC,icgc_array_roc2$AUC), c( pcsi_roc_se2, tcga_roc_se2, ouh_roc_se2,icgc_array_roc_se2),na.rm=TRUE,hetero=TRUE)$estimate
model2_seq = combine.est(c(pcsi_roc2$AUC, tcga_roc2$AUC), c( pcsi_roc_se2, tcga_roc_se2),na.rm=TRUE,hetero=TRUE)$estimate
model2_micro =combine.est(c( ouh_roc2$AUC,icgc_array_roc2$AUC), c( ouh_roc_se2,icgc_array_roc_se2),hetero=TRUE,na.rm=TRUE)$estimate

m11=combine.est(c(pcsi_roc1$AUC, tcga_roc1$AUC, ouh_roc1$AUC,icgc_array_roc1$AUC), c( pcsi_roc_se1, tcga_roc_se1, ouh_roc_se1,icgc_array_roc_se1),na.rm=TRUE,hetero=TRUE)
model1_meta_pval = pnorm((m11$estimate -0.5)/m11$se, lower.tail = m11$estimate < 0.5) * 2

m11=combine.est(c(pcsi_roc1$AUC, tcga_roc1$AUC), c( pcsi_roc_se1, tcga_roc_se1),na.rm=TRUE,hetero=TRUE)
model1_seq_pval <-  pnorm((m11$estimate -0.5)/m11$se, lower.tail = m11$estimate < 0.5) * 2

m11=combine.est(c( ouh_roc1$AUC,icgc_array_roc1$AUC), c( ouh_roc_se1,icgc_array_roc_se1),na.rm=FALSE,hetero=TRUE)
model1_micro_pval <-pnorm((m11$estimate -0.5)/m11$se, lower.tail = m11$estimate < 0.5) * 2

m12=combine.est(c(pcsi_roc2$AUC, tcga_roc2$AUC, ouh_roc2$AUC,icgc_array_roc2$AUC), c( pcsi_roc_se2, tcga_roc_se2, ouh_roc_se2,icgc_array_roc_se2),na.rm=TRUE,hetero=TRUE)
model2_meta_pval = pnorm((m12$estimate -0.5)/m12$se, lower.tail = m12$estimate < 0.5) * 2

m12=combine.est(c(pcsi_roc2$AUC, tcga_roc2$AUC), c( pcsi_roc_se2, tcga_roc_se2),na.rm=TRUE,hetero=TRUE)
model2_seq_pval <-  pnorm((m12$estimate -0.5)/m12$se, lower.tail = m12$estimate < 0.5) * 2

m12=combine.est(c( ouh_roc2$AUC,icgc_array_roc2$AUC), c( ouh_roc_se2,icgc_array_roc_se2),hetero=TRUE,na.rm=TRUE)
model2_micro_pval <-pnorm((m12$estimate -0.5)/m12$se, lower.tail = m12$estimate < 0.5) * 2





##################################
cohort_classes$pcsi[,1] == pcsi$ID

pcsi_classes= cohort_classes$pcsi[,2][which(cohort_classes$pcsi[,1] %in% pcsi$ID)]
cohort_classes$pcsi[,1][which(cohort_classes$pcsi[,1] %in% pcsi$ID)]==pcsi$ID
pcsi_grp=ifelse(as.numeric(as.character(pcsi$OS)) >=365,1,0)

tcga$ID=gsub( "_",".",tcga$ID)
tcga_classes= cohort_classes$tcga[,2][which(cohort_classes$tcga[,1] %in% tcga$ID)]
cohort_classes$tcga[,1][which(cohort_classes$tcga[,1] %in% tcga$ID)]==tcga$ID
tcga_grp=ifelse(as.numeric(as.character(tcga$OS)) >=365,1,0)


ouh_classes= cohort_classes$ouh[,2][which(cohort_classes$ouh[,1] %in% ouh$ID)]
cohort_classes$ouh[,1][which(cohort_classes$ouh[,1] %in% ouh$ID)]==ouh$ID
ouh_grp=ifelse(as.numeric(as.character(ouh$OS)) >=365,1,0)

icgc_arr_classes= cohort_classes$icgc_arr[,2][which(cohort_classes$icgc_arr[,1] %in% icgc_arr$ID)]
cohort_classes$icgc_arr[,1][which(cohort_classes$icgc_arr[,1] %in% icgc_arr$ID)]==icgc_arr$ID
icgc_arr_grp=ifelse(as.numeric(as.character(icgc_arr$OS)) >=365,1,0)


icgc_classes= cohort_classes$icgc[,2][which(cohort_classes$icgc[,1] %in% icgc$ID)]
cohort_classes$icgc[,1][which(cohort_classes$icgc[,1] %in% icgc$ID)]==icgc$ID
icgc_seq_grp=ifelse(as.numeric(as.character(icgc$OS)) >=365,1,0)



summary(glm(ouh_grp ~ Age + Sex + T_status+ N + M + Grade +pred_prob+ouh_classes, data=ouh, na.action = na.exclude))
summary(glm(pcsi_grp ~ Age + Sex + T_status+ N + M + Grade + pred_prob + pcsi_classes, data=pcsi, na.action = na.exclude))
summary(glm(tcga_grp ~ Age + Sex + T_status+ N + M + Grade + pred_prob +tcga_classes , data=tcga, na.action = na.exclude))
summary(glm(icgc_arr_grp ~ Age + Sex + T_status+ N + M + Grade + pred_prob + icgc_arr_classes, data=icgc_arr, na.action = na.exclude))
summary(glm(icgc_seq_grp ~ Age + Sex + T_status+ N + M + Grade + pred_prob + icgc_classes , data=icgc , na.action = na.exclude))










