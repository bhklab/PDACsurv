load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/clinical_features_ICGC_prob1.RData")
library(verification)
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
model1 <- lm(binary_grp ~ Age + Sex + T_status+ N + M + Grade, data=icgc , na.action = na.exclude)
model2 <- lm(binary_grp ~ pred_prob, data=icgc , na.action = na.exclude)
model3 <- lm(binary_grp ~ Age + Sex + T_status+ N + M + Grade +pred_prob, data=icgc , na.action = na.exclude)

summary(model1); anova(model1)
summary(model3); anova(model3)

####### Predicting probabilities in training cohort
icgc_cl_pred1=predict(model1,icgc, na.action = na.exclude)
icgc_cl_pred2=icgc$pred_prob[which(icgc$ID %in% names(icgc_cl_pred1))]
icgc_cl_pred3=predict(model3,icgc, na.action = na.exclude)

####### Predicting probabilities in validation cohorts

########## 1. PCSI

pcsi$Age=as.numeric(as.character(pcsi$Age))
pcsi$Sex = as.character(pcsi$Sex)
pcsi$binary_grp = as.numeric(as.character(pcsi$binary_grp))
pcsi$T_status = as.character(pcsi$T_status)


pcsi_cl_pred1=predict(model1, pcsi, na.action = na.exclude)
pcsi_cl_pred2=1-pcsi$pred_prob[which(pcsi$ID %in% names(pcsi_cl_pred1))]
pcsi_cl_pred3=predict(model3, pcsi, na.action = na.exclude)

########## 2. TCGA

tcga$Age=as.numeric(as.character(tcga$Age))
tcga$Sex = as.character(tcga$Sex)
tcga$binary_grp = as.numeric(as.character(tcga$binary_grp))
tcga$T_status = as.character(tcga$T_status)
tcga$pred_prob=as.numeric(as.character(tcga$pred_prob))
#tcga$Grade=as.character(tcga$Grade)

tcga_cl_pred1=predict(model1, tcga, na.action = na.exclude)
tcga_cl_pred2=1-tcga$pred_prob[which(tcga$ID %in% names(tcga_cl_pred1))]
tcga_cl_pred3=predict(model3, tcga, na.action = na.exclude)

########## 3. ICGC array

icgc_arr$Age=as.numeric(as.character(icgc_arr$Age))
icgc_arr$Sex = as.character(icgc_arr$Sex)
icgc_arr$binary_grp = as.numeric(as.character(icgc_arr$binary_grp))
icgc_arr$T_status = as.character(icgc_arr$T_status)


icgc_arr1=icgc_arr[order(icgc_arr$ID),]

icgc_arr_cl_pred1 = predict(model1, icgc_arr, na.action = na.exclude)
icgc_arr_cl_pred2 = 1-icgc_arr$pred_prob[which(icgc_arr1$ID %in% names(icgc_arr_cl_pred1))]
icgc_arr_cl_pred3 = predict(model3, icgc_arr, na.action = na.exclude)

########## 4.OUH

ouh$Age=as.numeric(as.character(ouh$Age))
ouh$Sex = as.character(ouh$Sex)
ouh$binary_grp = as.numeric(as.character(ouh$binary_grp))

ouh_cl_pred1 = predict(model1, ouh, na.action = na.exclude)
ouh_cl_pred2 = 1-ouh$pred_prob
ouh_cl_pred3 = predict(model3, ouh, na.action = na.exclude)


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

pcsi_roc3=reportROC(pcsi[names(pcsi_cl_pred3),]$binary_grp,pcsi_cl_pred3); pcsi_roc_se3=pcsi_roc3$AUC.SE
tcga_roc3=reportROC(tcga[names(tcga_cl_pred3),]$binary_grp,tcga_cl_pred3); tcga_roc_se3=tcga_roc3$AUC.SE
ouh_roc3=reportROC(ouh[names(ouh_cl_pred3),]$binary_grp,ouh_cl_pred3); ouh_roc_se3=ouh_roc3$AUC.SE
icgc_array_roc3 = reportROC(icgc_arr[names(icgc_arr_cl_pred3),]$binary_grp,icgc_arr_cl_pred3); icgc_array_roc_se3 = icgc_array_roc3$AUC.SE


pcsi_roc1_pval= roc.area(pcsi[names(pcsi_cl_pred1),]$binary_grp,pcsi_cl_pred1)$p.value  
tcga_roc1_pval= roc.area(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred1)$p.value  
ouh_roc1_pval= roc.area(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred1)$p.value  
icgc_array_roc1_pval= roc.area(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred1)$p.value  

pcsi_roc2_pval= roc.area(pcsi[names(pcsi_cl_pred1),]$binary_grp,1-pcsi_cl_pred2)$p.value  
tcga_roc2_pval= roc.area(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred2)$p.value  
ouh_roc2_pval= roc.area(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred2)$p.value  
icgc_array_roc2_pval= roc.area(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred2)$p.value  

pcsi_roc3_pval= roc.area(pcsi[names(pcsi_cl_pred3),]$binary_grp,pcsi_cl_pred3)$p.value  
tcga_roc3_pval= roc.area(tcga[names(tcga_cl_pred3),]$binary_grp,tcga_cl_pred3)$p.value  
ouh_roc3_pval= roc.area(ouh[names(ouh_cl_pred3),]$binary_grp,ouh_cl_pred3)$p.value  
icgc_array_roc3_pval= roc.area(icgc_arr[names(icgc_arr_cl_pred3),]$binary_grp,icgc_arr_cl_pred3)$p.value  


######## Plotting barplot
data <- structure(list("TCGA"= c(tcga_roc1$AUC, tcga_roc2$AUC,tcga_roc3$AUC ),
                       "PCSI" = c(pcsi_roc1$AUC, pcsi_roc2$AUC, pcsi_roc3$AUC), 
                       "ICGC-array" = c(icgc_array_roc1$AUC,icgc_array_roc2$AUC,icgc_array_roc3$AUC),
                       "OUH" = c(ouh_roc1$AUC,ouh_roc2$AUC, ouh_roc3$AUC)), 
                  .Names = c("TCGA","PCSI","ICGC-array","OUH"), class = "data.frame", row.names = c( "Clinicopathological model", "OS-TSP model", "Integrated clinical and OS-TSP model"))

attach(data)
print(data)

colours <- c( "lightcyan2",  "cyan3","gray" )
#pdf("/Users/vandanasandhu/Desktop/e.pdf")

barplot(as.matrix(data), main="", ylim= c(0, 0.8), ylab = "AUCs",  cex.main = 1.4, beside=TRUE, col=colours,border = 'NA',space=c(0.6,0.08,0.08,0.6,0.08,0.08,0.6,0.08,0.08,0.6,0.08,0.08))

########## Meta-estimate Calculations
model1_meta = combine.est(c(pcsi_roc1$AUC, tcga_roc1$AUC, ouh_roc1$AUC,icgc_array_roc1$AUC), c( pcsi_roc_se1, tcga_roc_se1, ouh_roc_se1,icgc_array_roc_se1),na.rm=TRUE,hetero=TRUE)$estimate
model1_seq = combine.est(c(pcsi_roc1$AUC, tcga_roc1$AUC), c( pcsi_roc_se1, tcga_roc_se1),na.rm=TRUE,hetero=FALSE)$estimate
model1_micro =combine.est(c( ouh_roc1$AUC,icgc_array_roc1$AUC), c( ouh_roc_se1,icgc_array_roc_se1),na.rm=FALSE,hetero=TRUE)$estimate

model2_meta = combine.est(c(pcsi_roc2$AUC, tcga_roc2$AUC, ouh_roc2$AUC,icgc_array_roc2$AUC), c( pcsi_roc_se2, tcga_roc_se2, ouh_roc_se2,icgc_array_roc_se2),na.rm=TRUE,hetero=TRUE)$estimate
model2_seq = combine.est(c(pcsi_roc2$AUC, tcga_roc2$AUC), c( pcsi_roc_se2, tcga_roc_se2),na.rm=TRUE,hetero=FALSE)$estimate
model2_micro =combine.est(c( ouh_roc2$AUC,icgc_array_roc2$AUC), c( ouh_roc_se2,icgc_array_roc_se2),hetero=FALSE,na.rm=TRUE)$estimate

model3_meta = combine.est(c(pcsi_roc3$AUC, tcga_roc3$AUC, ouh_roc3$AUC,icgc_array_roc3$AUC), c( pcsi_roc_se3, tcga_roc_se3, ouh_roc_se3,icgc_array_roc_se3),hetero=TRUE,na.rm=TRUE)$estimate
model3_seq = combine.est(c(pcsi_roc3$AUC, tcga_roc3$AUC), c( pcsi_roc_se3, tcga_roc_se3),na.rm=TRUE,hetero=FALSE)$estimate
model3_micro =combine.est(c( ouh_roc3$AUC,icgc_array_roc3$AUC), c( ouh_roc_se3,icgc_array_roc_se3),na.rm=TRUE,hetero=FALSE)$estimate

model1_meta_pval <- combine.test(p=c(pcsi_roc1_pval,tcga_roc1_pval, ouh_roc1_pval,icgc_array_roc1_pval),w=c(length(pcsi_cl_pred1),length(tcga_cl_pred1),length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),method="z.transform")
model1_seq_pval <- combine.test(p=c(pcsi_roc1_pval,tcga_roc1_pval),w=c(length(pcsi_cl_pred1),length(tcga_cl_pred1)),method="z.transform")
model1_micro_pval <- combine.test(p=c( ouh_roc1_pval,icgc_array_roc1_pval),w=c(length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),method="z.transform")

model2_meta_pval <- combine.test(p=c(pcsi_roc2_pval,tcga_roc2_pval, ouh_roc2_pval,icgc_array_roc2_pval),w=c(length(icgc_cl_pred2),length(tcga_cl_pred2),length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),method="z.transform")
model2_seq_pval <- combine.test(p=c(pcsi_roc2_pval,tcga_roc2_pval),w=c(length(icgc_cl_pred2),length(tcga_cl_pred2)),method="z.transform")
model2_micro_pval <- combine.test(p=c( ouh_roc2_pval,icgc_array_roc2_pval),w=c(length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),method="z.transform")

model3_meta_pval <- combine.test(p=c(pcsi_roc3_pval,tcga_roc3_pval, ouh_roc3_pval,icgc_array_roc3_pval),w=c(length(pcsi_cl_pred3),length(tcga_cl_pred3),length(ouh_cl_pred3),length(icgc_arr_cl_pred3)),method="z.transform")
model3_seq_pval <- combine.test(p=c(pcsi_roc3_pval,tcga_roc3_pval),w=c(length(pcsi_cl_pred3),length(tcga_cl_pred3)),method="z.transform")
model3_micro_pval <- combine.test(p=c( ouh_roc3_pval,icgc_array_roc3_pval),w=c(length(ouh_cl_pred3),length(icgc_arr_cl_pred3)),method="z.transform")

















