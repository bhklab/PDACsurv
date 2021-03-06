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
install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")
library("easyGgplot2")

load("training_cohorts.RData")
load("PDAC_Expression_dataset.RData")

icgc_seq_cohort = training_cohort$icgc_seq_cohort
icgc_array_cohort = training_cohort$icgc_array_cohort

rownames(icgc_array_cohort) == rownames(icgc_seq_cohort)

################## Excluding samples censored before 1-yr
g1=which(as.numeric(as.character(icgc_seq_cohort$OS))<=365 &  as.numeric(as.character(icgc_seq_cohort$OS_Status))==1);g2=which(as.numeric(as.character(icgc_seq_cohort$OS))>365)
g_ind=sort(c(g1,g2))

icgc_seq_cohort=icgc_seq_cohort[g_ind,]
icgc_array_cohort=icgc_array_cohort[g_ind,]

merge_common <- rbind(icgc_seq_cohort,icgc_array_cohort)      ### Merged common ICGC seq and array data

###################### Training the model on ICGC seq/array common samples cohort
## Classes for training
xx=merge_common
xmat<-xx[1:nrow(xx) ,1:(ncol(xx)-2)]                       ## Removing survival data columns
merge_common_mat <- data.matrix(sapply(xmat, function(xx) as.numeric(as.character(xx))))
rownames(merge_common_mat)=rownames(merge_common)
merge_common_grp=ifelse(as.numeric(as.character(xx$OS))>=365,1,0)

#########################################################
# Generating 1000 random models by re-shuffling the labels


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
  
  x5 <-sample(which(merge_common_grp==0), 40, replace=F)     # Selecting random 30 samples from group 1 
  y5 <-sample(which(merge_common_grp==1), 40, replace=F)     # Selecting random 30 samples from group 2
  
  x1=merge_common_mat[c(x5,y5),]                
  
  y_index=c(x5, y5)                                     # Selecting the classes of re-sampled samples
  shuffle_merge_grp=sample(merge_common_grp)            # Shuffling the labels
  y1= shuffle_merge_grp[y_index]
  
  
  ### Identifying KTSP models
  
  zzz=paste('classifier',i,sep="") 
  model[[i]]<- SWAP.KTSP.Train(t(x1), as.factor(y1) )
  print(i);
  
  
  z=setdiff(1:164,c(x5,y5))                              ### Finding test samples excluded in training set
  test=merge_common_mat[z,]   
  test_grp=merge_common_grp[z]
  
  
  ### Predicting on the test samples
  pred[[i]] <- SWAP.KTSP.Classify(t(test), model[[i]])   ### Predicting the classes of test set
  
  cc=confusionMatrix(pred[[i]], test_grp,  mode = "prec_recall")
  b_acc[i]=as.numeric(cc$byClass)[11]
}

reshuffle_model=model
save(reshuffle_model,file="/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/Random_Reshuffle_model.RData")


############# Random Reshuffle model function  #####################################
load("Random_Reshuffle_model.RData")
source("Validation_cohorts_formatting.R")

predict_ktsp = function(val_mat, val_grp){
  val_pred <- list()
  val_pred_freq<- list()
  auc=vector()
  auc_se=vector()
  for(i in 1: length(reshuffle_model) ){
    
    
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), reshuffle_model[[i]])  
    
    a=reportROC(val_grp, as.numeric(as.character(val_pred[[i]])),plot=FALSE)
    auc[i]=a$AUC
    auc_se[i]=a$AUC.SE
    
  }
  
  ret_list=list(val_pred,auc,auc_se)
  return(ret_list)  
}

############# Function call for independent cohorts #######################################
pcsi_mat=pcsi_mat[g_pcsi,]
chen_mat=chen_mat[g_chen,]
kirby_mat=kirby_mat[g_kirby,]
collisson_mat=collisson_mat[g_coll,]
ouh_mat=ouh_mat[g_ouh,]
winter_mat=winter_mat[g_winter,]
tcga_mat=tcga_mat[g_tcga,]
unc_mat=unc_mat[g_unc,]
zhang_mat=zhang_mat[g_zhang,]
icgc_array_mat=icgc_array_mat[g_icgc_arr,]

pcsi_list=predict_ktsp(pcsi_mat, pcsi_grp)
chen_list=predict_ktsp(chen_mat, chen_grp)
kirby_list=predict_ktsp(kirby_mat, kirby_grp)
collisson_list=predict_ktsp(collisson_mat, collisson_grp)
ouh_list=predict_ktsp(ouh_mat, ouh_grp)
winter_list=predict_ktsp(winter_mat, winter_grp)
tcga_list=predict_ktsp(tcga_mat, tcga_grp)
unc_list=predict_ktsp(unc_mat, unc_grp)
zhang_list=predict_ktsp(zhang_mat, zhang_grp)
icgc_arr_list=predict_ktsp(icgc_array_mat,icgc_array_grp )

###################### Calculating meta-estimates for all the 1000 models using all the cohorts

meta_auc=list()
seq_auc=list()
microarray_auc=list()

for( i in 1:1000){
  meta_auc[[i]] = combine.est(c(pcsi_list[[2]][i], tcga_list[[2]][i], unc_list[[2]][i], zhang_list[[2]][i],winter_list[[2]][i], ouh_list[[2]][i], icgc_arr_list[[2]][i], chen_list[[2]][i], kirby_list[[2]][i], collisson_list[[2]][i]), c( pcsi_list[[3]][i], tcga_list[[3]][i], unc_list[[3]][i], zhang_list[[3]][i], winter_list[[3]][i], ouh_list[[3]][i],icgc_arr_list[[3]][i], chen_list[[3]][i], kirby_list[[3]][i], collisson_list[[3]][i]),hetero=TRUE,na.rm=TRUE)$estimate
  seq_auc[[i]] = combine.est(c(pcsi_list[[2]][i], tcga_list[[2]][i], kirby_list[[2]][i]), c( pcsi_list[[3]][i], tcga_list[[3]][i],kirby_list[[3]][i]),na.rm=TRUE,hetero=TRUE)$estimate
  microarray_auc[[i]] = combine.est(c( unc_list[[2]][i], zhang_list[[2]][i],winter_list[[2]][i], ouh_list[[2]][i],icgc_arr_list[[2]][i], chen_list[[2]][i], collisson_list[[2]][i]), c( unc_list[[3]][i], zhang_list[[3]][i], winter_list[[3]][i], ouh_list[[3]][i],icgc_arr_list[[3]][i], chen_list[[3]][i], collisson_list[[3]][i]),na.rm=TRUE, hetero=TRUE)$estimate
  
}


######## Plotting the density plot 

platforms=c( rep("1. Array-based", 1000),  rep("2. Sequencing", 1000), rep("3. Overall",1000))
BAC=c(unlist(microarray_auc), unlist(seq_auc), unlist(meta_auc))
dd=data.frame(platforms=platforms, BAC=BAC)
vline.dat =data.frame(platforms=levels(dd$platforms), v1= c(0.69, 0.72, 0.70))

ggplot2.density(data=dd, xName='BAC', 
                groupName='platforms', legendPosition="top",
                faceting=TRUE, facetingVarNames="platforms",removePanelGrid=TRUE, removePanelBorder=TRUE, showLegend=FALSE, backgroundColor="white",
                fillGroupDensity = TRUE, colorGroupDensityLine = TRUE, 
                xTickLabelFont=c(6, "plain", "black"), yTickLabelFont=c(6, "plain", "black")) + 
  labs(x = "Balanced Accuracy", y="Density")+
  ggtitle("Random reshuffling of labels") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline( data=vline.dat, aes(xintercept=v1, color= platforms), linetype="dashed",size=1.5) 








