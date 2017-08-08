load("/Users/vandanasandhu/Documents/Projects/Project2_MetadataSubtyping/Survival_cohorts/OS_Predictor_updated_PDAC.RData")
ss=read.table("/Users/vandanasandhu/Documents/Projects/Project2_MetadataSubtyping/Subtyping-PCSI-OUH/icgc_array_survival.txt",sep="\t",header=T)


tran_val_coh=data.frame(t(val_coh))

##############  ICGC-sequencing data  ######################################################################
icgc_cohort=tran_val_coh[which(tran_val_coh$cohort %in% c("ICGC")),]

g1=which(as.numeric(as.character(icgc_cohort$OS))<=365  &  as.numeric(as.character(icgc_cohort$OS_Status))==1); g2=which(as.numeric(as.character(icgc_cohort$OS))>365)
g_ind=sort(c(g1,g2))

icgc_cohort=icgc_cohort[g_ind,]

icgc_mat <- data.matrix(sapply(matrix(icgc_cohort[1:nrow(icgc_cohort) ,4:ncol(icgc_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(icgc_mat)=rownames(icgc_cohort)
colnames(icgc_mat)=colnames(icgc_cohort)[4:ncol(icgc_cohort)]
icgc_grp=ifelse(as.numeric(as.character(icgc_cohort$OS))>=365,1,0)



##############  ICGC-array data  ######################################################################
icgc_array_cohort=tran_val_coh[which(tran_val_coh$cohort=="ICGC_array"),]

g1=which(as.numeric(as.character(icgc_array_cohort$OS))<=365 &  as.numeric(as.character(icgc_array_cohort$OS_Status))==1); g2=which(as.numeric(as.character(icgc_array_cohort$OS))>365)
g_ind=sort(c(g1,g2))

icgc_array_cohort=icgc_array_cohort[g_ind,]

icgc_array_mat <- data.matrix(sapply(matrix(icgc_array_cohort[1:nrow(icgc_array_cohort) ,4:ncol(icgc_array_cohort)]), function(xx) as.numeric(as.character(xx))))
 
rownames(icgc_array_mat)=as.character(ss[,2][g_ind])
colnames(icgc_array_mat)=colnames(icgc_array_cohort)[4:ncol(icgc_array_cohort)]

icgc_array_grp=ifelse(as.numeric(as.character(icgc_array_cohort$OS))>=365,1,0)


############# Predicting the probabilities using TSP-Ensembl model

icgc_arr_list=predict_ktsp(icgc_array_mat, icgc_array_grp)
icgc_seq_list=predict_ktsp(icgc_mat, icgc_grp)


### Finding common samples between ICGC-array and ICGC-sequencing

z1=which(rownames(icgc_array_mat) %in% rownames(icgc_mat))
z2=which(rownames(icgc_mat) %in% rownames(icgc_array_mat))

pr_arr=icgc_arr_list[[1]][z1]
pr_seq=icgc_seq_list[[1]][z2]
icgc_arr_grp=icgc_array_grp[z1]
icgc_seq_grp=icgc_grp[z2]

A=matrix(cbind(rownames(icgc_array_mat)[z1],pr_arr),ncol=2,nrow=82)
A=A[order(A[,1]),]
B=matrix(cbind(rownames(icgc_mat)[z2],pr_seq),ncol=2,nrow=82)
B=B[order(B[,1]),]

#### Correlation analysis and Scatter plots
Z=cbind(A,B)
plot(as.numeric(Z[,2]),as.numeric(Z[,4]),col="#6a3d9a",pch=19, xlab="ICGC array", ylab="ICGC-seq")
cor.test(as.numeric(Z[,2]),as.numeric(Z[,4]),method = "spearman")


################################### ROC curves for ICGC-array and ICGC-sequencing ##################################


icgc_array_roc=reportROC(icgc_array_grp,as.numeric(icgc_arr_list[[1]]))
icgc_array_roc_se=icgc_array_roc$AUC.SE

icgc_seq_roc=reportROC(icgc_grp,icgc_seq_list[[1]])
icgc_seq_roc_se=icgc_seq_roc$AUC.SE

plot(roc(icgc_array_grp,icgc_arr_list[[1]]),lwd=2,col="blue",lty=2)
plot(roc(icgc_grp,icgc_seq_list[[1]]),lwd=2, col="red",add=TRUE,lty=2)
legend("bottomright",legend=c(paste("ICGC-ARRAYS",round(icgc_array_roc$AUC,digits=2),sep=" "),
                              paste("ICGC-SEQ",round(icgc_seq_roc$AUC,digits=2),sep=" ")),
       
       fill=c("Blue","Red"),y.intersp = 0.7, cex=0.9)


