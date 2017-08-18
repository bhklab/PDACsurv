load("/Users/vandanasandhu/Desktop/os_predictor_project/clinical_features_seq_cohorts.RData")

############################## ###########################################

pcsi=data.frame(clinical_features$pcsi_clinical)
icgc=data.frame(clinical_features$icgc_clinical)
tcga=data.frame(clinical_features$tcga_clinical)
icgc_arr=data.frame(clinical_features$icgc_arr_clinical)
ouh=data.frame(clinical_features$ouh_clinical)

pcsi$Age=as.numeric(as.character(pcsi$Age))
pcsi$Sex = as.character(pcsi$Sex)
pcsi$binary_grp = as.numeric(as.character(pcsi$binary_grp))
pcsi$T_status=gsub(" ", "", pcsi$T_status, fixed = TRUE)


##### Model building comparisons
model1 <- lm(binary_grp ~ Age + Sex + T_status+ N + M + Grade, data=pcsi , na.action = na.exclude)
model2 <- lm(binary_grp ~ pred_prob, data=pcsi , na.action = na.exclude)
model3 <- lm(binary_grp ~ Age + Sex + T_status+ N + M + Grade +pred_prob, data=pcsi , na.action = na.exclude)


summary(model1); anova(model1)
summary(model2); anova(model2)
summary(model3); anova(model3)


####### Predicting probabilities in training cohort
pcsi_cl_pred1=predict(model1,pcsi, na.action = na.exclude)
pcsi_cl_pred2=predict(model2,pcsi, na.action = na.exclude)
pcsi_cl_pred3=predict(model3,pcsi, na.action = na.exclude)




####### Predicting probabilities in validation cohorts

########## 1. ICGC- sequencing

icgc$Age=as.numeric(as.character(icgc$Age))
icgc$Sex = as.character(icgc$Sex)
icgc$binary_grp = as.numeric(as.character(icgc$binary_grp))
icgc$T_status = as.character(icgc$T_status)

z1= which(icgc$T_status == "TX" )
z2= which(icgc$N == "NX" )
z3= which(icgc$G == "GX" )
z4= which(icgc$M == "M1" )
zz=c(z1,z2,z3,z4)

icgc=icgc[-zz,]
icgc_cl_pred1=predict(model1, icgc, na.action = na.exclude)
icgc_cl_pred2=predict(model2, icgc, na.action = na.exclude)
icgc_cl_pred3=predict(model3, icgc, na.action = na.exclude)

########## 2. TCGA

tcga$Age=as.numeric(as.character(tcga$Age))
tcga$Sex = as.character(tcga$Sex)
tcga$binary_grp = as.numeric(as.character(tcga$binary_grp))
tcga$T_status = as.character(tcga$T_status)
tcga$pred_prob=as.numeric(as.character(tcga$pred_prob))
#tcga$Grade=as.character(tcga$Grade)

tcga_cl_pred1=predict(model1, tcga, na.action = na.exclude)
tcga_cl_pred2=predict(model2, tcga, na.action = na.exclude)
tcga_cl_pred3=predict(model3, tcga, na.action = na.exclude)

########## 3. ICGC array

icgc_arr$Age=as.numeric(as.character(icgc_arr$Age))
icgc_arr$Sex = as.character(icgc_arr$Sex)
icgc_arr$binary_grp = as.numeric(as.character(icgc_arr$binary_grp))
icgc_arr$T_status = as.character(icgc_arr$T_status)

z1= which(icgc_arr$N == "NX" )
z2= which(icgc_arr$M == "M1" )
z3= which(icgc_arr$G == "GX" )
z4= which(icgc_arr$T_status == "TX" )

zz=c(z1,z2,z3,z4)

icgc_arr=icgc_arr[-zz,]
icgc_arr_cl_pred1 = predict(model1, icgc_arr, na.action = na.exclude)
icgc_arr_cl_pred2 = predict(model2, icgc_arr, na.action = na.exclude)
icgc_arr_cl_pred3 = predict(model3, icgc_arr, na.action = na.exclude)

########## 4.OUH

ouh$Age=as.numeric(as.character(ouh$Age))
ouh$Sex = as.character(ouh$Sex)
ouh$binary_grp = as.numeric(as.character(ouh$binary_grp))

ouh_cl_pred1 = predict(model1, ouh, na.action = na.exclude)
ouh_cl_pred2 = predict(model2, ouh, na.action = na.exclude)
ouh_cl_pred3 = predict(model3, ouh, na.action = na.exclude)


#MODEL 1: Clinical model AUC comparison ####################
###########################################################################
pcsi_roc1=reportROC(pcsi[names(pcsi_cl_pred1),]$binary_grp,pcsi_cl_pred1); pcsi_roc_se1=pcsi_roc1$AUC.SE
icgc_roc1=reportROC(icgc[names(icgc_cl_pred1),]$binary_grp,icgc_cl_pred1); icgc_roc_se1=icgc_roc1$AUC.SE
tcga_roc1=reportROC(tcga[names(tcga_cl_pred1),]$binary_grp,tcga_cl_pred1); tcga_roc_se1=tcga_roc1$AUC.SE
ouh_roc1=reportROC(ouh[names(ouh_cl_pred1),]$binary_grp,ouh_cl_pred1); ouh_roc_se1=ouh_roc1$AUC.SE
icgc_array_roc1 = reportROC(icgc_arr[names(icgc_arr_cl_pred1),]$binary_grp,icgc_arr_cl_pred1); icgc_array_roc_se1 = icgc_array_roc1$AUC.SE

pcsi_roc2=reportROC(pcsi[names(pcsi_cl_pred2),]$binary_grp,pcsi_cl_pred2); pcsi_roc_se2=pcsi_roc2$AUC.SE
icgc_roc2=reportROC(icgc[names(icgc_cl_pred2),]$binary_grp,icgc_cl_pred2); icgc_roc_se2=icgc_roc2$AUC.SE
tcga_roc2=reportROC(tcga[names(tcga_cl_pred2),]$binary_grp,tcga_cl_pred2); tcga_roc_se2=tcga_roc2$AUC.SE
ouh_roc2=reportROC(ouh[names(ouh_cl_pred2),]$binary_grp,ouh_cl_pred2); ouh_roc_se2=ouh_roc2$AUC.SE
icgc_array_roc2 = reportROC(icgc_arr[names(icgc_arr_cl_pred2),]$binary_grp,icgc_arr_cl_pred2);icgc_array_roc_se2 = icgc_array_roc2$AUC.SE

pcsi_roc3=reportROC(pcsi[names(pcsi_cl_pred3),]$binary_grp,pcsi_cl_pred3); pcsi_roc_se3=pcsi_roc3$AUC.SE
icgc_roc3=reportROC(icgc[names(icgc_cl_pred3),]$binary_grp,icgc_cl_pred3); icgc_roc_se3=icgc_roc3$AUC.SE
tcga_roc3=reportROC(tcga[names(tcga_cl_pred3),]$binary_grp,tcga_cl_pred3); tcga_roc_se3=tcga_roc3$AUC.SE
ouh_roc3=reportROC(ouh[names(ouh_cl_pred3),]$binary_grp,ouh_cl_pred3); ouh_roc_se3=ouh_roc3$AUC.SE
icgc_array_roc3 = reportROC(icgc_arr[names(icgc_arr_cl_pred3),]$binary_grp,icgc_arr_cl_pred3); icgc_array_roc_se3 = icgc_array_roc3$AUC.SE


data <- structure(list("TCGA"= c(tcga_roc1$AUC, tcga_roc2$AUC,tcga_roc3$AUC ),
                       "ICGC-sequencing" = c(icgc_roc1$AUC, icgc_roc2$AUC, icgc_roc3$AUC), 
                       "ICGC-array" = c(icgc_array_roc1$AUC,icgc_array_roc2$AUC,icgc_array_roc3$AUC),
                       "GSE60980" = c(ouh_roc1$AUC,ouh_roc2$AUC, ouh_roc3$AUC)), 
                  .Names = c("TCGA","ICGC-sequencing","ICGC-array","GSE60980"), class = "data.frame", row.names = c( "Clinicopathological model", "OS-TSP model", "Integrated clinical and OS-TSP model"))

attach(data)
print(data)

colours <- c(  "black", "yellow", "gray")
pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure4a.pdf")

barplot(as.matrix(data), main="", ylab = "AUCs", å, cex.main = 1.4, beside=TRUE, col=colours,border = 'NA')
dev.off()






model1_meta= combine.est(c(icgc_roc1$AUC, tcga_roc1$AUC, ouh_roc1$AUC,icgc_array_roc1$AUC), c( icgc_roc_se1, tcga_roc_se1, ouh_roc_se1,icgc_array_roc_se1),na.rm=TRUE)$estimate
model1_seq= combine.est(c(icgc_roc1$AUC, tcga_roc1$AUC), c( icgc_roc_se1, tcga_roc_se1),na.rm=TRUE)$estimate
model1_micro=combine.est(c( ouh_roc1$AUC,icgc_array_roc1$AUC), c( ouh_roc_se1,icgc_array_roc_se1),na.rm=TRUE)$estimate
model2_meta= combine.est(c(icgc_roc2$AUC, tcga_roc2$AUC, ouh_roc2$AUC,icgc_array_roc2$AUC), c( icgc_roc_se2, tcga_roc_se2, ouh_roc_se2,icgc_array_roc_se2),na.rm=TRUE)$estimate
model2_seq= combine.est(c(icgc_roc2$AUC, tcga_roc2$AUC), c( icgc_roc_se2, tcga_roc_se2),na.rm=TRUE)$estimate
model2_micro=combine.est(c( ouh_roc2$AUC,icgc_array_roc2$AUC), c( ouh_roc_se2,icgc_array_roc_se2),na.rm=TRUE)$estimate
model3_meta= combine.est(c(icgc_roc3$AUC, tcga_roc3$AUC, ouh_roc3$AUC,icgc_array_roc3$AUC), c( icgc_roc_se3, tcga_roc_se3, ouh_roc_se3,icgc_array_roc_se3),na.rm=TRUE)$estimate
model3_seq= combine.est(c(icgc_roc3$AUC, tcga_roc3$AUC), c( icgc_roc_se3, tcga_roc_se3),na.rm=TRUE)$estimate
model3_micro=combine.est(c( ouh_roc3$AUC,icgc_array_roc3$AUC), c( ouh_roc_se3,icgc_array_roc_se3),na.rm=TRUE)$estimate


