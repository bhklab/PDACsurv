load("/Users/vandanasandhu/Documents/Projects/Project2_MetadataSubtyping/Survival_cohorts/OS_Predictor_Aug_15_2017.RData")

##############  ICGC-sequencing data  ######################################################################
icgc_cohort=val_coh$ICGC_seq

g1=which(as.numeric(as.character(icgc_cohort$OS))<=365  &  as.numeric(as.character(icgc_cohort$OS_Status))==1); g2=which(as.numeric(as.character(icgc_cohort$OS))>365)
g_ind=sort(c(g1,g2))

icgc_cohort=icgc_cohort[g_ind,]

icgc_mat <- data.matrix(sapply(matrix(icgc_cohort[1:nrow(icgc_cohort) ,4:ncol(icgc_cohort)]), function(xx) as.numeric(as.character(xx))))
rownames(icgc_mat)=rownames(icgc_cohort)
colnames(icgc_mat)=colnames(icgc_cohort)[4:ncol(icgc_cohort)]
icgc_grp=ifelse(as.numeric(as.character(icgc_cohort$OS))>=365,1,0)



##############  ICGC-array data  ######################################################################
icgc_array_cohort=val_coh$ICGC_array_all
g1=which(as.numeric(as.character(icgc_array_cohort$OS))<=365 &  as.numeric(as.character(icgc_array_cohort$OS_Status))==1); g2=which(as.numeric(as.character(icgc_array_cohort$OS))>365)
g_ind=sort(c(g1,g2))

icgc_array_cohort=icgc_array_cohort[g_ind,]

icgc_array_mat <- data.matrix(sapply(matrix(icgc_array_cohort[1:nrow(icgc_array_cohort) ,4:ncol(icgc_array_cohort)]), function(xx) as.numeric(as.character(xx))))
 
rownames(icgc_array_mat)=rownames(icgc_array_cohort)
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
pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure2c.pdf")

Z=cbind(A,B)
plot(as.numeric(Z[,2]),as.numeric(Z[,4]),col="#6a3d9a",pch=19, xlab="ICGC array", ylab="ICGC-sequencing")
cor.test(as.numeric(Z[,2]),as.numeric(Z[,4]),method = "spearman")
legend("bottomright",legend=c( paste("Spearman correlation R:",0.9,sep=" "),
                               paste("P value:",2.2e-16,sep=" ")), bty='n')

################################### ROC curves for ICGC-array and ICGC-sequencing ##################################


icgc_array_roc=reportROC(icgc_arr_grp,pr_arr)
icgc_array_roc_se=icgc_array_roc$AUC.SE

icgc_seq_roc=reportROC(icgc_seq_grp,pr_seq)
icgc_seq_roc_se=icgc_seq_roc$AUC.SE

pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure2b.pdf")

plot(roc(icgc_array_grp,icgc_arr_list[[1]]),lwd=4,col="turquoise3",lty=3)
plot(roc(icgc_grp,icgc_seq_list[[1]]),lwd=4, col="#fb9a99",add=TRUE,lty=1)
legend("bottomright",legend=c( paste("ICGC-sequencing",round(icgc_seq_roc$AUC,digits=2),sep=" "),
                              paste("ICGC-array",round(icgc_array_roc$AUC,digits=2),sep=" ")),
       fill=c("#fb9a99","turquoise3"),y.intersp = 1, cex=1,bty = "n")




