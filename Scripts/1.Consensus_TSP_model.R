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

load("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Data_R_objects/RS_updated_model_ba.RData")
load("/Users/vandanasandhu/Documents/Projects/Project2_MetadataSubtyping/Survival_cohorts/OS_Predictor_Aug_15_2017.RData")


######################
###################### Training the model on PCSI cohort

pcsi_cohort=val_coh$PCSI

## Classes for training
g1=which(as.numeric(as.character(pcsi_cohort$OS))<=365 &  as.numeric(as.character(pcsi_cohort$OS_Status))==1)
g2=which(as.numeric(as.character(pcsi_cohort$OS))>365)
g_ind=sort(c(g1,g2))



xx=pcsi_cohort[g_ind,]
xmat<-xx[1:nrow(xx) ,4:ncol(xx)]
pcsi_mat <- data.matrix(sapply(xmat, function(xx) as.numeric(as.character(xx))))
rownames(pcsi_mat)=xx[,1][1:nrow(xmat)]
pcsi_grp=ifelse(as.numeric(as.character(xx$OS))>=365,1,0)

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

for(i in 1:models_no){

  x5 <-sample(which(pcsi_grp==0), 30, replace=F)     # Selecting random 30 samples from group 1 
  y5 <-sample(which(pcsi_grp==1), 30, replace=F)     # Selecting random 30 samples from group 2
  
  x1=pcsi_mat[c(x5,y5),]                
  
  y_index=c(x5, y5)                                  # Selecting the classes of re-sampled samples
  y1=pcsi_grp[y_index]
  
  
  ### Identifying KTSPs
  
  zzz=paste('classifier',i,sep="") 
  model[[i]]<- SWAP.KTSP.Train(t(x1), as.factor(y1) )
 
  z=setdiff(1:113,c(x5,y5))                              ### Finding test samples excluded in training set
  test=pcsi_mat[z,]   
  test_grp=pcsi_grp[z]
  
 
  ### Predicting on the test samples
  pred[[i]] <- SWAP.KTSP.Classify(t(test), model[[i]])   ### Predicting the classes of test set
  cc=confusionMatrix(pred[[i]], test_grp)
  b_acc[i]=as.numeric(cc$byClass)[11]
  print(i)
  
  #### Filtering model with balance accuracy more than 0.6
 
  if(b_acc[i] > 0.6){
    
    selected_model[[count]]=model[[i]]
    sel_pred[[count]] <- SWAP.KTSP.Classify(t(test), selected_model[[count]])   
    sel_b_acc[[count]] <- b_acc[i]
    count=count+1
  }
  
} 


save(selected_model,"consensus_model.RData")
