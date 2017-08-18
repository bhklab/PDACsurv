load("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Data_R_objects/RS_updated_model_ba.RData")
load("/Users/vandanasandhu/Documents/Projects/Project2_MetadataSubtyping/Survival_cohorts/OS_Predictor_Aug_15_2017.RData")

##################  ##################  ##################  ##################  
##################  Predict function for validation in independent cohorts

predict_ktsp = function(val_mat, val_grp){
  val_pred <- list()
  
  for(i in 1: length(selected_model) ){
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), selected_model[[i]])  
  }
  
  list_z1=vector()
  freq_early=vector()
  i=1
  
  for (i in 1: nrow(val_mat)){
    for (k in 1: length(selected_model) ){
      list_z1=append(list_z1,as.numeric(val_pred[[k]][i]))
    }
    
    freq_early[i] =length(list_z1[list_z1==1])/length(selected_model) 
    list_z1=vector()
  }
  ret_list=list(predicted_probabilities=freq_early)
  return(ret_list)
}

#### Validation cohorts for AUC estinates 
###################### ICGC
icgc_cohort=val_coh$ICGC_seq

g1=which(as.numeric(as.character(icgc_cohort$OS))<=365 &  as.numeric(as.character(icgc_cohort$OS_Status))==1);g2=which(as.numeric(as.character(icgc_cohort$OS))>365)
g_ind=sort(c(g1,g2))
icgc_cohort=icgc_cohort[g_ind,]

icgc_mat<-data.matrix(sapply(icgc_cohort[1:nrow(icgc_cohort) ,4:ncol(icgc_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(icgc_mat)=rownames(icgc_cohort); icgc_grp=ifelse(as.numeric(as.character(icgc_cohort$OS))>=365,1,0)

###################### TCGA
tcga_cohort=val_coh$TCGA

g1=which(as.numeric(as.character(tcga_cohort$OS))<=365 &  as.numeric(as.character(tcga_cohort$OS_Status))==1);g2=which(as.numeric(as.character(tcga_cohort$OS))>365)
g_ind=sort(c(g1,g2))
tcga_cohort=tcga_cohort[g_ind,]

tcga_mat<-data.matrix(sapply(tcga_cohort[1:nrow(tcga_cohort) ,4:ncol(tcga_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(tcga_mat)=rownames(tcga_cohort); tcga_grp=ifelse(as.numeric(as.character(tcga_cohort$OS))>=365,1,0)

###################### Moffitt
moff_cohort=val_coh$Moffitt

g1=which(as.numeric(as.character(moff_cohort$OS))<=365 &  as.numeric(as.character(moff_cohort$OS_Status))==1);g2=which(as.numeric(as.character(moff_cohort$OS))>365)
g_ind=sort(c(g1,g2))
moff_cohort=moff_cohort[g_ind,]

moff_mat<-data.matrix(sapply(moff_cohort[1:nrow(moff_cohort) ,4:ncol(moff_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(moff_mat)=rownames(moff_cohort); moff_grp=ifelse(as.numeric(as.character(moff_cohort$OS))>=365,1,0)

###################### OUH
ouh_cohort=val_coh$OUH

g1=which(as.numeric(as.character(ouh_cohort$OS))<=365 &  as.numeric(as.character(ouh_cohort$OS_Status))==1);g2=which(as.numeric(as.character(ouh_cohort$OS))>365)
g_ind=sort(c(g1,g2))
ouh_cohort=ouh_cohort[g_ind,]

ouh_mat<-data.matrix(sapply(ouh_cohort[1:nrow(ouh_cohort) ,4:ncol(ouh_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(ouh_mat)=rownames(ouh_cohort); ouh_grp=ifelse(as.numeric(as.character(ouh_cohort$OS))>=365,1,0)

###################### zhang
zhang_cohort=val_coh$Zhang

g1=which(as.numeric(as.character(zhang_cohort$OS))<=365 &  as.numeric(as.character(zhang_cohort$OS_Status))==1);g2=which(as.numeric(as.character(zhang_cohort$OS))>365)
g_ind=sort(c(g1,g2))
zhang_cohort=zhang_cohort[g_ind,]

zhang_mat<-data.matrix(sapply(zhang_cohort[1:nrow(zhang_cohort) ,4:ncol(zhang_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(zhang_mat)=rownames(zhang_cohort); zhang_grp=ifelse(as.numeric(as.character(zhang_cohort$OS))>=365,1,0)

###################### Winter
winter_cohort=val_coh$Winter

g1=which(as.numeric(as.character(winter_cohort$OS))<=365 &  as.numeric(as.character(winter_cohort$OS_Status))==1);g2=which(as.numeric(as.character(winter_cohort$OS))>365)
g_ind=sort(c(g1,g2))
winter_cohort=winter_cohort[g_ind,]

winter_mat<-data.matrix(sapply(winter_cohort[1:nrow(winter_cohort) ,4:ncol(winter_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(winter_mat)=rownames(winter_cohort); winter_grp=ifelse(as.numeric(as.character(winter_cohort$OS))>=365,1,0)

###################### ICGC_array
icgc_array_cohort=val_coh$ICGC_array

g1=which(as.numeric(as.character(icgc_array_cohort$OS))<=365 &  as.numeric(as.character(icgc_array_cohort$OS_Status))==1);g2=which(as.numeric(as.character(icgc_array_cohort$OS))>365)
g_ind=sort(c(g1,g2))
icgc_array_cohort=icgc_array_cohort[g_ind,]

icgc_array_mat<-data.matrix(sapply(icgc_array_cohort[1:nrow(icgc_array_cohort) ,4:ncol(icgc_array_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(icgc_array_mat)=rownames(icgc_array_cohort); icgc_array_grp=ifelse(as.numeric(as.character(icgc_array_cohort$OS))>=365,1,0)

###################### pcsi
pcsi_cohort=val_coh$PCSI

g1=which(as.numeric(as.character(pcsi_cohort$OS))<=365 &  as.numeric(as.character(pcsi_cohort$OS_Status))==1);g2=which(as.numeric(as.character(pcsi_cohort$OS))>365)
g_ind=sort(c(g1,g2))
pcsi_cohort=pcsi_cohort[g_ind,]

pcsi_mat<-data.matrix(sapply(pcsi_cohort[1:nrow(pcsi_cohort) ,4:ncol(pcsi_cohort)], function(xx) as.numeric(as.character(xx))))
rownames(pcsi_mat)=rownames(pcsi_cohort); pcsi_grp=ifelse(as.numeric(as.character(pcsi_cohort$OS))>=365,1,0)



######## Predicting prediction probabilities using TSP-Ensembl model

pcsi_list=predict_ktsp(pcsi_mat, pcsi_grp)
icgc_list=predict_ktsp(icgc_mat, icgc_grp)
tcga_list=predict_ktsp(tcga_mat, tcga_grp)
moff_list=predict_ktsp(moff_mat, moff_grp)
ouh_list=predict_ktsp(ouh_mat, ouh_grp)
zhang_list=predict_ktsp(zhang_mat, zhang_grp)
winter_list=predict_ktsp(winter_mat, winter_grp)
icgc_array_list=predict_ktsp(icgc_array_mat, icgc_array_grp)


######################
###################### Plotting ROC CURVES

pcsi_roc=reportROC(pcsi_grp,pcsi_list[[1]]); pcsi_roc_se=pcsi_roc$AUC.SE
icgc_roc=reportROC(icgc_grp,icgc_list[[1]]); icgc_roc_se=icgc_roc$AUC.SE
tcga_roc=reportROC(tcga_grp,tcga_list[[1]]); tcga_roc_se=tcga_roc$AUC.SE
moff_roc=reportROC(moff_grp,moff_list[[1]]); moff_roc_se=moff_roc$AUC.SE
zhang_roc=reportROC(zhang_grp,zhang_list[[1]]); zhang_roc_se=zhang_roc$AUC.SE
winter_roc=reportROC(winter_grp,winter_list[[1]]); winter_roc_se=winter_roc$AUC.SE
ouh_roc=reportROC(ouh_grp,ouh_list[[1]]); ouh_roc_se=ouh_roc$AUC.SE
icgc_array_roc=reportROC(icgc_array_grp,icgc_array_list[[1]]); icgc_array_roc_se=icgc_array_roc$AUC.SE

library(pROC)
pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure2a.pdf")
plot(roc(icgc_grp,icgc_list[[1]]),lwd=4, col="#fb9a99",lty=1)
plot(roc(tcga_grp,tcga_list[[1]]),lwd=4, col="chartreuse3",add=TRUE,lty=1)
plot(roc(moff_grp,moff_list[[1]]),lwd=4, col="darkgoldenrod1",add=TRUE,lty=3)
plot(roc(zhang_grp,zhang_list[[1]]),lwd=4, col="wheat4",add=TRUE,lty=3)
plot(roc(winter_grp,winter_list[[1]]),lwd=4, col="cornflowerblue",add=TRUE,lty=3)
plot(roc(ouh_grp,ouh_list[[1]]),lwd=4, col="mediumorchid2",add=TRUE,lty=3)
plot(roc(icgc_array_grp,icgc_array_list[[1]]),lwd=4, col="turquoise3",add=TRUE,lty=3)

legend("bottomright",legend=c(paste("TCGA :",round(tcga_roc$AUC,digits=2),sep=" "),
                              paste("ICGC-sequencing :",round(icgc_roc$AUC,digits=2),sep=" "),
                              paste("ICGC-array :",round(icgc_array_roc$AUC,digits=2),sep=" "),
                              paste("GSE71729 :",round(moff_roc$AUC,digits=2),sep=" "),
                              paste("GSE28735 :",round(zhang_roc$AUC,digits=2),sep=" "),
                              paste("E-MEXP-2780 :",round(winter_roc$AUC,digits=2),sep=" "),
                              paste("GSE60980 :",round(ouh_roc$AUC,digits=2),sep=" ")),
                             
       fill=c("chartreuse3","#fb9a99","turquoise3","darkgoldenrod1","wheat4","cornflowerblue","mediumorchid2"),y.intersp = 1, cex=0.9,bty = "n")


meta_auc = combine.est(c(icgc_roc$AUC, tcga_roc$AUC, moff_roc$AUC, zhang_roc$AUC, winter_roc$AUC, ouh_roc$AUC,icgc_array_roc$AUC), c( icgc_roc_se, tcga_roc_se, moff_roc_se,zhang_roc_se, winter_roc_se, ouh_roc_se,icgc_array_roc_se),na.rm=TRUE)$estimate
seq_auc = combine.est(c(icgc_roc$AUC, tcga_roc$AUC), c( icgc_roc_se, tcga_roc_se),na.rm=TRUE)$estimate
micro_auc = combine.est(c( moff_roc$AUC, zhang_roc$AUC, winter_roc$AUC, ouh_roc$AUC,icgc_array_roc$AUC), c(  moff_roc_se,zhang_roc_se, winter_roc_se, ouh_roc_se,icgc_array_roc_se),na.rm=TRUE)$estimate

#######################Probability distribution boxplots

df=cbind(pred=c(pcsi_list[[1]],icgc_list[[1]],tcga_list[[1]],moff_list[[1]],zhang_list[[1]],winter_list[[1]],ouh_list[[1]],icgc_array_list[[1]]),grp=c(pcsi_grp, icgc_grp, tcga_grp, moff_grp, zhang_grp, winter_grp, ouh_grp,icgc_array_grp), coh=c(rep("1.PCSI",length(pcsi_grp)),rep("2.ICGC",length(icgc_grp)),rep("3.TCGA",length(tcga_grp)),rep("4.Moffitt",length(moff_grp)),rep("5.Zhang",length(zhang_grp)),rep("6.Winter",length(winter_grp)),rep("7.OUH",length(ouh_grp)),rep("8.ICGC_array",length(icgc_array_grp))))

df.m=as.data.frame(df)

df.m$pred=as.numeric(as.character(df.m$pred))
require(ggplot2)
ggplot(data = df.m, aes(x=coh, y=pred)) + geom_boxplot(aes(fill=grp)) + geom_point(aes(y=pred, group=grp), position = position_dodge(width=0.75))


################################### Validation cohorts for D-index and C-index, including censored samples as well
###################### ICGC
icgc_cohort=val_coh$ICGC_seq
xx=icgc_cohort
g1=which(as.numeric(as.character(icgc_cohort$OS))<=365 ); g2=which(as.numeric(as.character(icgc_cohort$OS))>365)
g_ind=sort(c(g1,g2))
icgc_cohort=icgc_cohort[g_ind,]

icgc_mat<-data.matrix(sapply(matrix(icgc_cohort[1:nrow(icgc_cohort) ,4:ncol(icgc_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(icgc_mat)=icgc_cohort[,1]
colnames(icgc_mat)=colnames(icgc_cohort)[4:ncol(xx)]
icgc_grp=ifelse(as.numeric(as.character(icgc_cohort$OS))>=365,1,0)

###################### tcga
tcga_cohort=val_coh$TCGA

xx=tcga_cohort
g1=which(as.numeric(as.character(tcga_cohort$OS))<=365 ); g2=which(as.numeric(as.character(tcga_cohort$OS))>365)
g_ind=sort(c(g1,g2))
tcga_cohort=tcga_cohort[g_ind,]

tcga_mat<-data.matrix(sapply(matrix(tcga_cohort[1:nrow(tcga_cohort) ,4:ncol(tcga_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(tcga_mat)=tcga_cohort[,1]
colnames(tcga_mat)=colnames(tcga_cohort)[4:ncol(xx)]
tcga_grp=ifelse(as.numeric(as.character(tcga_cohort$OS))>=365,1,0)

###################### Moffitt
moff_cohort=val_coh$Moffitt

xx=moff_cohort
g1=which(as.numeric(as.character(moff_cohort$OS))<=365 ); g2=which(as.numeric(as.character(moff_cohort$OS))>365)
g_ind=sort(c(g1,g2))
moff_cohort=moff_cohort[g_ind,]

moff_mat<-data.matrix(sapply(matrix(moff_cohort[1:nrow(moff_cohort) ,4:ncol(moff_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(moff_mat)=moff_cohort[,1]
colnames(moff_mat)=colnames(moff_cohort)[4:ncol(xx)]
moff_grp=ifelse(as.numeric(as.character(moff_cohort$OS))>=365,1,0)

###################### OUH
ouh_cohort=val_coh$OUH
xx=ouh_cohort
g1=which(as.numeric(as.character(ouh_cohort$OS))<=365 ); g2=which(as.numeric(as.character(ouh_cohort$OS))>365)
g_ind=sort(c(g1,g2))
ouh_cohort=ouh_cohort[g_ind,]

ouh_mat<-data.matrix(sapply(matrix(ouh_cohort[1:nrow(ouh_cohort) ,4:ncol(ouh_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(ouh_mat)=ouh_cohort[,1]
colnames(ouh_mat)=colnames(ouh_cohort)[4:ncol(xx)]
ouh_grp=ifelse(as.numeric(as.character(ouh_cohort$OS))>=365,1,0)

###################### zhang
zhang_cohort=val_coh$Zhang
#zhang_cohort=tran_val_coh[which(tran_val_coh$cohort=="Zhang"),]
xx=zhang_cohort
g1=which(as.numeric(as.character(zhang_cohort$OS))<=365 ); g2=which(as.numeric(as.character(zhang_cohort$OS))>365)
g_ind=sort(c(g1,g2))
zhang_cohort=zhang_cohort[g_ind,]

zhang_mat<-data.matrix(sapply(matrix(zhang_cohort[1:nrow(zhang_cohort) ,4:ncol(zhang_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(zhang_mat)=zhang_cohort[,1]
colnames(zhang_mat)=colnames(zhang_cohort)[4:ncol(xx)]
zhang_grp=ifelse(as.numeric(as.character(zhang_cohort$OS))>=365,1,0)

###################### Winter
winter_cohort=val_coh$Winter
xx=winter_cohort
g1=which(as.numeric(as.character(winter_cohort$OS))<=365 ); g2=which(as.numeric(as.character(winter_cohort$OS))>365)
g_ind=sort(c(g1,g2))
winter_cohort=winter_cohort[g_ind,]

winter_mat<-data.matrix(sapply(matrix(winter_cohort[1:nrow(winter_cohort) ,4:ncol(winter_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(winter_mat)=winter_cohort[,1]
colnames(winter_mat)=colnames(winter_cohort)[4:ncol(xx)]
winter_grp=ifelse(as.numeric(as.character(winter_cohort$OS))>=365,1,0)

###################### ICGC_array
icgc_array_cohort=val_coh$ICGC_array
xx=icgc_array_cohort
g1=which(as.numeric(as.character(icgc_array_cohort$OS))<=365 ); g2=which(as.numeric(as.character(icgc_array_cohort$OS))>365)
g_ind=sort(c(g1,g2))
icgc_array_cohort=icgc_array_cohort[g_ind,]

icgc_array_mat<-data.matrix(sapply(matrix(icgc_array_cohort[1:nrow(icgc_array_cohort) ,4:ncol(icgc_array_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(icgc_array_mat)=icgc_array_cohort[,1]
colnames(icgc_array_mat)=colnames(icgc_array_cohort)[4:ncol(xx)]
icgc_array_grp=ifelse(as.numeric(as.character(icgc_array_cohort$OS))>=365,1,0)


######## Predicting probabilities for all samples in validation cohort
icgc_list=predict_ktsp(icgc_mat, icgc_grp)
tcga_list=predict_ktsp(tcga_mat, tcga_grp)
moff_list=predict_ktsp(moff_mat, moff_grp)
ouh_list=predict_ktsp(ouh_mat, ouh_grp)
zhang_list=predict_ktsp(zhang_mat, zhang_grp)
winter_list=predict_ktsp(winter_mat, winter_grp)
icgc_array_list=predict_ktsp(icgc_array_mat, icgc_array_grp)



### Dindex estimate calculation

dindex_ouh <- D.index(x=ouh_list[[1]], surv.time=as.numeric(as.character(ouh_cohort$OS)), surv.event=as.numeric(as.character(ouh_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_icgc <- D.index(x=icgc_list[[1]], surv.time=as.numeric(as.character(icgc_cohort$OS)), surv.event=as.numeric(as.character(icgc_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_tcga <- D.index(x=tcga_list[[1]], surv.time=as.numeric(as.character(tcga_cohort$OS)), surv.event=as.numeric(as.character(tcga_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_winter <- D.index(x=winter_list[[1]], surv.time=as.numeric(as.character(winter_cohort$OS)), surv.event=as.numeric(as.character(winter_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_zhang <- D.index(x=zhang_list[[1]], surv.time=as.numeric(as.character(zhang_cohort$OS)), surv.event=as.numeric(as.character(zhang_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_icgc_array <- D.index(x=icgc_array_list[[1]], surv.time=as.numeric(as.character(icgc_array_cohort$OS)), surv.event=as.numeric(as.character(icgc_array_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_moff <- D.index(x=moff_list[[1]], surv.time=as.numeric(as.character(moff_cohort$OS)), surv.event=as.numeric(as.character(moff_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");

### Concordance index calculation
con_ouh <- concordance.index(x=ouh_list[[1]], surv.time=as.numeric(as.character(ouh_cohort$OS)), surv.event=as.numeric(as.character(ouh_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_icgc <- concordance.index(x=icgc_list[[1]], surv.time=as.numeric(as.character(icgc_cohort$OS)), surv.event=as.numeric(as.character(icgc_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_tcga <- concordance.index(x=tcga_list[[1]], surv.time=as.numeric(as.character(tcga_cohort$OS)), surv.event=as.numeric(as.character(tcga_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_winter <- concordance.index(x=winter_list[[1]], surv.time=as.numeric(as.character(winter_cohort$OS)), surv.event=as.numeric(as.character(winter_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_zhang <- concordance.index(x=zhang_list[[1]], surv.time=as.numeric(as.character(zhang_cohort$OS)), surv.event=as.numeric(as.character(zhang_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_moff <- concordance.index(x=moff_list[[1]], surv.time=as.numeric(as.character(moff_cohort$OS)), surv.event=as.numeric(as.character(moff_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_icgc_array <- concordance.index(x=icgc_array_list[[1]], surv.time=as.numeric(as.character(icgc_array_cohort$OS)), surv.event=as.numeric(as.character(icgc_array_cohort$OS_Status)), na.rm=TRUE, method="noether");

####################################### Calculating meta-estimates of D-index and Concordance index


###  Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR OVERALL DATA
dindex_meta <- combine.est(c( dindex_icgc$d.index, dindex_tcga$d.index, dindex_ouh$d.index,dindex_winter$d.index,dindex_zhang$d.index,dindex_moff$d.index,dindex_icgc_array$d.index),c(dindex_icgc$se, dindex_tcga$se, dindex_ouh$se,dindex_winter$se,dindex_zhang$se,dindex_moff$se,dindex_icgc_array$se),na.rm = TRUE)
dindex_meta_lower <- dindex_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta$se
dindex_meta_upper <- dindex_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta$se
dindex_meta_pval <- combine.test(p=c(dindex_icgc$p.value, dindex_tcga$p.value, dindex_ouh$p.value,dindex_winter$p.value,dindex_zhang$p.value,dindex_moff$p.value,dindex_icgc_array$p.value),w=c(length(icgc_list[[1]]), length(tcga_list[[1]]),length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),length(moff_list[[1]]),length(icgc_array_list[[1]])),hetero = FALSE,method="z.transform")


con_meta <- combine.est(c( con_icgc$c.index, con_tcga$c.index, con_ouh$c.index,con_winter$c.index,con_zhang$c.index,con_moff$c.index,con_icgc_array$c.index),c(con_icgc$se, con_tcga$se, con_ouh$se,con_winter$se,con_zhang$se,con_moff$se,con_icgc_array$se),na.rm = TRUE)
con_meta_lower <- con_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta$se
con_meta_upper <- con_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta$se
con_meta_pval <- combine.test(p=c(con_icgc$p.value, con_tcga$p.value, con_ouh$p.value,con_winter$p.value,con_zhang$p.value,con_moff$p.value,con_icgc_array$p.value),w=c(length(icgc_list[[1]]), length(tcga_list[[1]]),length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),length(moff_list[[1]]),length(icgc_array_list[[1]])),hetero = FALSE,method="z.transform")


### Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR sequencing cohort

dindex_seq <- combine.est(c( dindex_icgc$d.index, dindex_tcga$d.index),c(dindex_icgc$se, dindex_tcga$se),na.rm = TRUE)
dindex_seq_lower <- dindex_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq$se
dindex_seq_upper <- dindex_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq$se
dindex_seq_pval <- combine.test(p=c(dindex_icgc$p.value, dindex_tcga$p.value),w=c(length(icgc_list[[1]]), length(tcga_list[[1]])),hetero = FALSE,method="z.transform")

con_seq <- combine.est(c( con_icgc$c.index, con_tcga$c.index),c(con_icgc$se, con_tcga$se),na.rm = TRUE)
con_seq_lower <- con_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq$se
con_seq_upper <- con_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq$se
con_seq_pval <- combine.test(p=c(con_icgc$p.value, con_tcga$p.value),w=c(length(icgc_list[[1]]), length(tcga_list[[1]])),hetero = FALSE,method="z.transform")

### Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR microarray cohort

dindex_micro <- combine.est(c( dindex_ouh$d.index,dindex_winter$d.index,dindex_zhang$d.index,dindex_moff$d.index,dindex_icgc_array$d.index),c( dindex_ouh$se,dindex_winter$se,dindex_zhang$se,dindex_moff$se,dindex_icgc_array$se),na.rm = TRUE)
dindex_micro_lower <- dindex_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro$se
dindex_micro_upper <- dindex_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro$se
dindex_micro_pval <- combine.test(p=c(dindex_ouh$p.value,dindex_winter$p.value,dindex_zhang$p.value,dindex_moff$p.value,dindex_icgc_array$p.value),w=c(length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),length(moff_list[[1]]),length(icgc_array_list[[1]])),hetero = FALSE,method="z.transform")


con_micro <- combine.est(c(  con_ouh$c.index,con_winter$c.index,con_zhang$c.index,con_moff$c.index,con_icgc_array$c.index),c( con_ouh$se,con_winter$se,con_zhang$se,con_moff$se,con_icgc_array$se),na.rm = TRUE)
con_micro_lower <- con_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro$se
con_micro_upper <- con_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro$se
con_micro_pval <- combine.test(p=c(con_ouh$p.value,con_winter$p.value,con_zhang$p.value,con_moff$p.value,con_icgc_array$p.value),w=c(length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),length(moff_list[[1]]),length(icgc_array_list[[1]])),hetero = FALSE,method="z.transform")



##### Plotting Forestplot of  D index ###############
r.mean <- c( log2(dindex_tcga$d.index),log2(dindex_icgc$d.index), log2(dindex_icgc_array$d.index), log2(dindex_moff$d.index),log2(dindex_zhang$d.index),log2(dindex_winter$d.index), log2(dindex_ouh$d.index),  log2(dindex_seq$estimate),log2(dindex_micro$estimate),log2(dindex_meta$estimate))
r.lower <- c( log2(dindex_tcga$lower), log2(dindex_icgc$lower), log2(dindex_icgc_array$lower), log2(dindex_moff$lower),  log2(dindex_zhang$lower), log2(dindex_winter$lower), log2(dindex_ouh$lower),log2(dindex_seq_lower),log2(dindex_micro_lower), log2(dindex_meta_lower))
r.upper <- c( log2(dindex_tcga$upper),log2(dindex_icgc$upper),log2(dindex_icgc_array$upper),  log2(dindex_moff$upper),  log2(dindex_zhang$upper), log2(dindex_winter$upper), log2(dindex_ouh$upper),log2(dindex_seq_upper), log2(dindex_micro_upper),log2(dindex_meta_upper))
r.pval <- round(c(dindex_tcga$p.value, dindex_icgc$p.value, dindex_icgc_array$p.value,dindex_moff$p.value,dindex_zhang$p.value,   dindex_winter$p.value, dindex_ouh$p.value,dindex_seq_pval,dindex_micro_pval, dindex_meta_pval),2)

r.pval1 <- c(sprintf("%.1E", dindex_tcga$p.value),sprintf("%.1E", dindex_icgc$p.value),sprintf("%.1E", dindex_icgc_array$p.value),sprintf("%.1E", dindex_moff$p.value), sprintf("%.1E", dindex_zhang$p.value), sprintf("%.1E", dindex_winter$p.value), sprintf("%.1E", dindex_ouh$p.value),sprintf("%.1E", dindex_seq_pval),sprintf("%.1E", dindex_micro_pval),sprintf("%.1E", dindex_meta_pval))

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","ICGC-sequencing","ICGC-array","GSE71729","GSE28735","E-MEXP-2780","GSE60980","Sequencing","Microarray","Overall")

data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts",rownames(t)),
  c("P values",r.pval1))
pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure3a.pdf")

forestplot(tabletext2,data2,xlab="Log2 D-index",is.summary=c(TRUE,rep(FALSE,7),TRUE, TRUE, TRUE),clip=c(-3,3.0),cex=9,col = fpColors(lines="royalblue", box="darkblue", summary ="darkred",text="black"),title="D-Index",zero=0,graphwidth=unit(2, "inches"),  align=c("l"))
dev.off()

##### Plotting Forestplot of Concordance index ###############

r.mean <- c(con_tcga$c.index, con_icgc$c.index, con_icgc_array$c.index,  con_moff$c.index,con_zhang$c.index,  con_winter$c.index, con_ouh$c.index, con_seq$estimate,con_micro$estimate,con_meta$estimate)
r.lower <- c( con_tcga$lower,  con_icgc$lower, con_icgc_array$lower,con_moff$lower,con_zhang$lower, con_winter$lower,con_ouh$lower,  con_seq_lower,con_micro_lower, con_meta_lower)
r.upper <- c( con_tcga$upper, con_icgc$upper,  con_icgc_array$upper,con_moff$upper,  con_zhang$upper,con_winter$upper,con_ouh$upper, con_seq_upper, con_micro_upper, con_meta_upper)
r.pval <- round(c( con_tcga$p.value,con_icgc$p.value, con_icgc_array$p.value, con_moff$p.value, con_zhang$p.value, con_winter$p.value, con_ouh$p.value, con_seq_pval,con_micro_pval,con_meta_pval),2)
r.pval1 <- c(sprintf("%.1E", con_tcga$p.value),sprintf("%.1E", con_icgc$p.value),sprintf("%.1E", con_icgc_array$p.value), sprintf("%.1E", con_moff$p.value),sprintf("%.1E", con_zhang$p.value),sprintf("%.1E", con_winter$p.value), sprintf("%.1E", con_ouh$p.value),sprintf("%.1E", con_seq_pval),sprintf("%.1E", con_micro_pval),sprintf("%.1E", con_meta_pval))

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","ICGC-sequencing","ICGC-array","GSE71729","GSE28735","E-MEXP-2780","GSE60980","Sequencing","Microarray","Overall")


data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -11L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts",rownames(t)),
  c("P values",r.pval1))
pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure3b.pdf")


forestplot(tabletext2,data2,xlab="Concordance index",is.summary=c(TRUE,rep(FALSE,7),TRUE, TRUE, TRUE), clip=c(0,1),col = fpColors(lines="royalblue", box="darkblue", summary ="darkred"),title="Concordance-index",zero=0.5,graphwidth=unit(2, "inches"),align=c("l"))
dev.off()

#txt_gp =  fpTxtGp(label = gpar(fontfamily = "Verdana"),ticks = gpar(cex=0.8)), 

#####
