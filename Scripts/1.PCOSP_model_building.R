
########## Required libraries

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
library(verification)

#### Load datasets
load("training_cohorts.RData")

icgc_seq_cohort = training_cohort$icgc_seq_cohort
icgc_array_cohort = training_cohort$icgc_array_cohort

rownames(icgc_array_cohort) == rownames(icgc_seq_cohort)

###### Excluding samples censored before 1-yr
g1=which(as.numeric(as.character(icgc_seq_cohort$OS))<=365 &  as.numeric(as.character(icgc_seq_cohort$OS_Status))==1);g2=which(as.numeric(as.character(icgc_seq_cohort$OS))>365)
g_ind=sort(c(g1,g2))

icgc_seq_cohort=icgc_seq_cohort[g_ind,]
icgc_array_cohort=icgc_array_cohort[g_ind,]

merge_common <- rbind(icgc_seq_cohort,icgc_array_cohort)      ### Merged common ICGC seq and array data



######################
###################### Training the model on ICGC seq/array common samples cohort
## Classes for training

xx=merge_common
xmat<-xx[1:nrow(xx) ,1:(ncol(xx)-2)]                       ## Removing survival data columns
merge_common_mat <- data.matrix(sapply(xmat, function(xx) as.numeric(as.character(xx))))
rownames(merge_common_mat)=rownames(merge_common)
merge_common_grp=ifelse(as.numeric(as.character(xx$OS))>=365,1,0)

#######################################################
######################## Generating 1000 TSP models

pred <- list()
sel_pred <- list()
count=0;
b_acc <- vector()
F1 <- vector()
i=1
model <- list()
models_no=1000
count=1
selected_model=list()
set.seed(1987)
sel_b_acc=list()

for(i in 1:1000){
  
  x5 <-sample(which(merge_common_grp==0), 40, replace=F)     # Selecting random 40 samples from Low survival group
  y5 <-sample(which(merge_common_grp==1), 40, replace=F)     # Selecting random 40 samples from High survival group
  
  x1=merge_common_mat[c(x5,y5),]                
  y_index=c(x5, y5)                                  # Selecting the classes of re-sampled samples
  y1=merge_common_grp[y_index]
  
  
  ### Building k-TSP models
  zzz=paste('classifier',i,sep="") 
  model[[i]]<- SWAP.KTSP.Train(t(x1), as.factor(y1) )
  
  z=setdiff(1:164,c(x5,y5))                              ### Finding test samples excluded in training set
  test=merge_common_mat[z,]   
  test_grp=merge_common_grp[z]
  
  ### Testing the model on out of bag samples
  pred[[i]] <- SWAP.KTSP.Classify(t(test), model[[i]])   ### Predicting the classes of test set
  
  cc=confusionMatrix(pred[[i]], test_grp,  mode = "prec_recall")
  b_acc[i]=as.numeric(cc$byClass)[11]
  F1[i]=as.numeric(cc$byClass)[7]
  print(i)
}

selected_model = model[which(b_acc> 0.60)]
length(selected_model)
save(selected_model, file="/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/PCOSP.RData")


