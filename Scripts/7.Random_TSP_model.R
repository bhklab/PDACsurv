require(switchBox)
library(vcdExtra)
library(caret)
library(forestplot)
library("ktspair")
library(pROC)
library(survcomp)
library(survival)
args(SWAP.Train.KTSP)
require(data.table)
library(reportROC)

load("/Users/vandanasandhu/Documents/Projects/Project2_MetadataSubtyping/Survival_cohorts/OS_Predictor_Aug_15_2017.RData")
load("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Data_R_objects/random_classifier.RData")


###################### PCSI cohort as training cohort
######################

pcsi_cohort=val_coh$PCSI

g1=which(as.numeric(as.character(pcsi_cohort$OS))<=365 &  as.numeric(as.character(pcsi_cohort$OS_Status))==1); g2=which(as.numeric(as.character(pcsi_cohort$OS))>365)
g_ind=sort(c(g1,g2))

xx=pcsi_cohort[g_ind,]
pcsi_cohort=xx

pcsi_mat <-data.matrix(sapply(xx[1:nrow(xx) ,4:ncol(xx)], function(xx) as.numeric(as.character(xx))))
rownames(pcsi_mat)=xx[,1][1:nrow(pcsi_mat)]
pcsi_grp=ifelse(as.numeric(as.character(pcsi_cohort$OS))>=365,1,0)

#########################################################
# Generating 1000 random TSP models 


pred <- list()
sel_pred <- list()
count=0;
b_acc <- vector()
i=1
model <- list()
models_no=1000
count=1
selected_model=list()
set.seed(1987)
sel_b_acc=list()

for(i in 1:models_no){
  #set.seed(1)
  x5 <-sample(which(pcsi_grp==0), 30, replace=F)     # Selecting random 30 samples from group 1 
  y5 <-sample(which(pcsi_grp==1), 30, replace=F)     # Selecting random 30 samples from group 2
  
  x1=pcsi_mat[c(x5,y5),]                
  
  y_index=c(x5, y5)                   # Selecting the classes of re-sampled samples
  y1=pcsi_grp[y_index]
  
  
  ### k-TSP models
  
  zzz=paste('classifier',i,sep="") 
  model[[i]]<- SWAP.KTSP.Train(t(x1), as.factor(y1) )
  
  selected_model[[count]]=model[[i]]
  selected_model[[count]]$TSPs[,1]=sample(colnames(pcsi_mat),length(selected_model[[count]]$TSPs[,1]))
  selected_model[[count]]$TSPs[,2]=sample(colnames(pcsi_mat),length(selected_model[[count]]$TSPs[,1]))
  
  count=count+1;
  print(i)
}


################ Predicting probabilities function using random classifier model

predict_ktsp = function(val_mat, val_grp){
  val_pred <- list()
  val_pred_freq<- list()
  auc=vector()
  auc_se=vector()
  
  for(i in 1: length(model) ){
    
    
    
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), selected_model[[i]])  
    
    a=reportROC(val_grp, as.numeric(as.character(val_pred[[i]])))
    auc[i]=a$AUC
    auc_se[i]=a$SE
    
  }
  
  ret_list=list(val_pred,auc,auc_se)
  return(ret_list)  
}



######## Predicting probabliliets using random model for all the validation cohorts

pcsi_list=predict_ktsp(pcsi_mat, pcsi_grp)
icgc_list=predict_ktsp(icgc_mat, icgc_grp)
tcga_list=predict_ktsp(tcga_mat, tcga_grp)
moff_list=predict_ktsp(moff_mat, moff_grp)
ouh_list=predict_ktsp(ouh_mat, ouh_grp)
zhang_list=predict_ktsp(zhang_mat, zhang_grp)
winter_list=predict_ktsp(winter_mat, winter_grp)
icgc_arr_list=predict_ktsp(icgc_array_mat,icgc_array_grp )
###################### Calculating meta-estimates for all the 1000 models using all the cohorts

meta_auc=list()
seq_auc=list()
microarray_auc=list()

for( i in 1:1000){
  meta_auc[[i]] = combine.est(c(icgc_list[[2]][i], tcga_list[[2]][i], moff_list[[2]][i], zhang_list[[2]][i],winter_list[[2]][i], ouh_list[[2]][i], icgc_arr_list[[2]][i]), c( icgc_list[[3]][i], tcga_list[[3]][i], moff_list[[3]][i], zhang_list[[3]][i], winter_list[[3]][i], ouh_list[[3]][i],icgc_arr_list[[3]][i]),na.rm=TRUE)$estimate
  seq_auc[[i]] = combine.est(c(icgc_list[[2]][i], tcga_list[[2]][i]), c( icgc_list[[3]][i], tcga_list[[3]][i]),na.rm=TRUE)$estimate
  microarray_auc[[i]] = combine.est(c( moff_list[[2]][i], zhang_list[[2]][i],winter_list[[2]][i], ouh_list[[2]][i],icgc_arr_list[[2]][i]), c( moff_list[[3]][i], zhang_list[[3]][i], winter_list[[3]][i], ouh_list[[3]][i],icgc_arr_list[[3]][i]),na.rm=TRUE)$estimate
  
}

######## Plotting the density plot 

pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure6b.pdf")
myData <-data.frame(Overall=unlist(meta_auc), Sequencing =unlist(seq_auc), Microarray= unlist(microarray_auc)) 
data=melt(myData)
colnames(data)=c("Platforms","AUCs")
ggplot(data,aes(x=AUCs, fill=Platforms)) + geom_density(alpha = 0.3, position = "stack") +  geom_vline(aes(xintercept=0.72), color="pink", linetype="dashed", size=1)+ geom_vline(aes(xintercept=0.77), color="green", linetype="dashed", size=1) + geom_vline(aes(xintercept=0.68), color="lightblue", linetype="dashed", size=1)

