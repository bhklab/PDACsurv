######### Load libraries

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

#### Load PCOSP model and Validation cohort
setwd("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/")
load("PCOSP.RData")
source("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/Scripts/PCOSP_score_estimation.R")
source("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/Scripts/Validation_cohorts_formatting.R")

######## PCOSP score estimations for all validation cohort
pcsi_list=pcosp_prob(pcsi_mat)
icgc_list=pcosp_prob(icgc_mat)
tcga_list=pcosp_prob(tcga_mat)
icgc_array_list=pcosp_prob(icgc_array_mat)
ouh_list=pcosp_prob(ouh_mat)
zhang_list=pcosp_prob(zhang_mat)
winter_list=pcosp_prob(winter_mat)
unc_list=pcosp_prob(unc_mat)
collisson_list=pcosp_prob(collisson_mat)
chen_list=pcosp_prob(chen_mat)
kirby_list=pcosp_prob(kirby_mat)

######################
### Dindex estimate calculation

dindex_ouh <- D.index(x=ouh_list[[1]], surv.time=as.numeric(as.character(ouh_cohort$OS)), 
              surv.event=as.numeric(as.character(ouh_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_icgc <- D.index(x=icgc_list[[1]], surv.time=as.numeric(as.character(icgc_cohort$OS)), 
               surv.event=as.numeric(as.character(icgc_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_tcga <- D.index(x=tcga_list[[1]], surv.time=as.numeric(as.character(tcga_cohort$OS)), 
               surv.event=as.numeric(as.character(tcga_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_winter <- D.index(x=winter_list[[1]], surv.time=as.numeric(as.character(winter_cohort$OS)),
                surv.event=as.numeric(as.character(winter_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_zhang <- D.index(x=zhang_list[[1]], surv.time=as.numeric(as.character(zhang_cohort$OS)), 
                surv.event=as.numeric(as.character(zhang_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_icgc_array <- D.index(x=icgc_array_list[[1]], surv.time=as.numeric(as.character(icgc_array_cohort$OS)),
                    surv.event=as.numeric(as.character(icgc_array_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_unc <- D.index(x=unc_list[[1]], surv.time=as.numeric(as.character(unc_cohort$OS)), 
              surv.event=as.numeric(as.character(unc_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_collisson <- D.index(x=collisson_list[[1]], surv.time=as.numeric(as.character(collisson_cohort$OS)), 
                    surv.event=as.numeric(as.character(collisson_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_chen <- D.index(x=chen_list[[1]], surv.time=as.numeric(as.character(chen_cohort$OS)), 
               surv.event=as.numeric(as.character(chen_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_kirby <- D.index(x=kirby_list[[1]], surv.time=as.numeric(as.character(kirby_cohort$OS)),
                surv.event=as.numeric(as.character(kirby_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
dindex_pcsi <- D.index(x=pcsi_list[[1]], surv.time=as.numeric(as.character(pcsi_cohort$OS)), 
               surv.event=as.numeric(as.character(pcsi_cohort$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");


########### Concordance index calculation
con_ouh <- concordance.index(x=ouh_list[[1]], surv.time=as.numeric(as.character(ouh_cohort$OS)), 
                             surv.event=as.numeric(as.character(ouh_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_icgc <- concordance.index(x=icgc_list[[1]], surv.time=as.numeric(as.character(icgc_cohort$OS)), 
                              surv.event=as.numeric(as.character(icgc_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_tcga <- concordance.index(x=tcga_list[[1]], surv.time=as.numeric(as.character(tcga_cohort$OS)), 
                              surv.event=as.numeric(as.character(tcga_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_winter <- concordance.index(x=winter_list[[1]], surv.time=as.numeric(as.character(winter_cohort$OS)), 
                                surv.event=as.numeric(as.character(winter_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_zhang <- concordance.index(x=zhang_list[[1]], surv.time=as.numeric(as.character(zhang_cohort$OS)), 
                               surv.event=as.numeric(as.character(zhang_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_unc <- concordance.index(x=unc_list[[1]], surv.time=as.numeric(as.character(unc_cohort$OS)), 
                             surv.event=as.numeric(as.character(unc_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_icgc_array <- concordance.index(x=icgc_array_list[[1]], surv.time=as.numeric(as.character(icgc_array_cohort$OS)), 
                                    surv.event=as.numeric(as.character(icgc_array_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_collisson<- concordance.index(x=collisson_list[[1]], surv.time=as.numeric(as.character(collisson_cohort$OS)), 
                                  surv.event=as.numeric(as.character(collisson_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_chen<- concordance.index(x=chen_list[[1]], surv.time=as.numeric(as.character(chen_cohort$OS)), 
                             surv.event=as.numeric(as.character(chen_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_kirby<- concordance.index(x=kirby_list[[1]], surv.time=as.numeric(as.character(kirby_cohort$OS)), 
                              surv.event=as.numeric(as.character(kirby_cohort$OS_Status)), na.rm=TRUE, method="noether");
con_pcsi<- concordance.index(x=pcsi_list[[1]], surv.time=as.numeric(as.character(pcsi_cohort$OS)), 
                             surv.event=as.numeric(as.character(pcsi_cohort$OS_Status)), na.rm=TRUE, method="noether");

####################################### Calculating meta-estimates of D-index and Concordance index


###  Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR OVERALL DATA
dindex_meta <- combine.est(c( dindex_pcsi$d.index, dindex_tcga$d.index, dindex_kirby$d.index, 
                              dindex_ouh$d.index, dindex_winter$d.index,dindex_zhang$d.index,
                              dindex_unc$d.index,dindex_icgc_array$d.index,dindex_collisson$d.index, dindex_chen$d.index),
                            c(dindex_pcsi$se, dindex_tcga$se, dindex_kirby$se, 
                              dindex_ouh$se,dindex_winter$se, dindex_zhang$se,
                              dindex_unc$se,dindex_icgc_array$se,dindex_collisson$se, 
                              dindex_chen$se),na.rm = TRUE,hetero = TRUE)
dindex_meta_lower <- dindex_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta$se
dindex_meta_upper <- dindex_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta$se
# #dindex_meta_pval <- combine.test(p=c(dindex_pcsi$p.value, dindex_tcga$p.value, dindex_kirby$p.value, 
#                                      dindex_ouh$p.value,dindex_winter$p.value,dindex_zhang$p.value,
#                                      dindex_unc$p.value,dindex_icgc_array$p.value, dindex_collisson$p.value,dindex_chen$p.value),
#                                  
#                                  w=c(length(pcsi_list[[1]]), length(tcga_list[[1]]),length(kirby_list[[1]]),
#                                      length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),
#                                      length(unc_list[[1]]),length(icgc_array_list[[1]]), length(collisson_list[[1]]), 
#                                      length(chen_list[[1]])),hetero = FALSE,method="z.transform")


dindex_meta_pval= 2*pnorm(-abs(log(dindex_meta$estimate)/dindex_meta$se))

#pnorm((dindex_meta$estimate - 0.5)/dindex_meta$se, lower.tail = dindex_meta$estimate < 0.5) * 2


con_meta <- combine.est(c( con_pcsi$c.index, con_tcga$c.index, con_kirby$c.index, 
                           con_ouh$c.index,con_winter$c.index,con_zhang$c.index,
                           con_unc$c.index,con_icgc_array$c.index, con_collisson$c.index,  con_chen$c.index),
                        
                        c(con_pcsi$se, con_tcga$se, con_kirby$se,
                          con_ouh$se, con_winter$se,con_zhang$se,
                          con_unc$se,con_icgc_array$se,con_collisson$se,  con_chen$se),na.rm = TRUE,hetero = TRUE)

con_meta_lower <- con_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta$se
con_meta_upper <- con_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta$se
# con_meta_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value, con_kirby$p.value,
#                                   con_ouh$p.value, con_winter$p.value, con_zhang$p.value,
#                                   con_unc$p.value, con_icgc_array$p.value, con_collisson$p.value, con_chen$p.value),
#                               
#                               w=c(length(pcsi_list[[1]]), length(tcga_list[[1]]),length(kirby_list[[1]]), 
#                                   length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),
#                                   length(unc_list[[1]]),length(icgc_array_list[[1]]), length(collisson_list[[1]]), 
#                                   length(chen_list[[1]])),method="z.transform")
# 
con_meta_pval=pnorm((con_meta$estimate - 0.5)/con_meta$se, lower.tail = con_meta$estimate < 0.5) * 2


### Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR sequencing cohort

dindex_seq <- combine.est(c( dindex_pcsi$d.index, dindex_tcga$d.index, dindex_kirby$d.index),
                          c(dindex_pcsi$se, dindex_tcga$se, dindex_kirby$se),na.rm = TRUE, hetero = FALSE) ## Since the platform is same
dindex_seq_lower <- dindex_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq$se
dindex_seq_upper <- dindex_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq$se
# dindex_seq_pval <- combine.test(p=c(dindex_pcsi$p.value, dindex_tcga$p.value, dindex_kirby$p.value),
#                                 w=c(length(pcsi_list[[1]]), length(tcga_list[[1]]), length(kirby_list[[1]])),method="z.transform")

dindex_seq_pval<-2*pnorm(-abs(log(dindex_seq$estimate)/dindex_seq$se))


#pnorm((dindex_seq$estimate - 0.5)/dindex_seq$se, lower.tail = dindex_seq$estimate < 0.5) * 2


con_seq <- combine.est(c( con_pcsi$c.index, con_tcga$c.index, con_kirby$c.index),
                       c(con_pcsi$se, con_tcga$se, con_kirby$se),na.rm = TRUE,  hetero = FALSE)## Since the platform is same
con_seq_lower <- con_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq$se
con_seq_upper <- con_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq$se
# con_seq_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value, con_kirby$p.value),
#                              w=c(length(pcsi_list[[1]]), length(tcga_list[[1]]), length(kirby_list[[1]])),method="z.transform")
con_seq_pval<-pnorm((con_seq$estimate - 0.5)/con_seq$se, lower.tail = con_seq$estimate < 0.5) * 2

### Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR microarray cohort

dindex_micro <- combine.est(c( dindex_ouh$d.index,dindex_winter$d.index,dindex_zhang$d.index,
                               dindex_unc$d.index,dindex_icgc_array$d.index, dindex_collisson$d.index, dindex_chen$d.index),
                            c( dindex_ouh$se,dindex_winter$se,dindex_zhang$se,
                               dindex_unc$se,dindex_icgc_array$se, dindex_collisson$se, dindex_chen$se),na.rm = TRUE, hetero = FALSE)

dindex_micro_lower <- dindex_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro$se
dindex_micro_upper <- dindex_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro$se
# dindex_micro_pval <- combine.test(p=c(dindex_ouh$p.value,dindex_winter$p.value,dindex_zhang$p.value,
#                                       dindex_unc$p.value,dindex_icgc_array$p.value, dindex_collisson$p.value, dindex_chen$p.value),
#                                   w=c(length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),
#                                       length(unc_list[[1]]),length(icgc_array_list[[1]]), length(collisson_list[[1]]), 
#                                       length(chen_list[[1]])),hetero = FALSE,method="z.transform")

dindex_micro_pval<- 2*pnorm(-abs(log(dindex_micro$estimate)/dindex_micro$se))

#pnorm((dindex_micro$estimate - 0.5)/dindex_micro$se, lower.tail = dindex_micro$estimate < 0.5) * 2
  
con_micro <- combine.est(c(  con_ouh$c.index,con_winter$c.index,con_zhang$c.index,
                             con_unc$c.index,con_icgc_array$c.index, con_collisson$c.index,  con_chen$c.index),
                         c( con_ouh$se,con_winter$se,con_zhang$se,
                            con_unc$se,con_icgc_array$se, con_collisson$se, con_chen$se),na.rm = TRUE, hetero = FALSE)
con_micro_lower <- con_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro$se
con_micro_upper <- con_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro$se
# con_micro_pval <- combine.test(p=c(con_ouh$p.value,con_winter$p.value,con_zhang$p.value,con_unc$p.value,
#                                    con_icgc_array$p.value, con_collisson$p.value, con_chen$p.value),
#                                w=c(length(ouh_list[[1]]),length(winter_list[[1]]),length(zhang_list[[1]]),length(unc_list[[1]]),
#                                    length(icgc_array_list[[1]]), length(collisson_list[[1]]),  length(chen_list[[1]])),
#                                hetero = FALSE,method="z.transform")

con_micro_pval<-pnorm((con_micro$estimate - 0.5)/con_micro$se, lower.tail = con_micro$estimate < 0.5) * 2



##### Plotting Forestplot of  D index ###############
r.mean <- c( log2(dindex_tcga$d.index), log2(dindex_pcsi$d.index), log2(dindex_kirby$d.index),  
             log2(dindex_icgc_array$d.index), log2(dindex_unc$d.index),log2(dindex_chen$d.index),
             log2(dindex_ouh$d.index),  log2(dindex_zhang$d.index), log2(dindex_winter$d.index),  log2(dindex_collisson$d.index),  
             log2(dindex_seq$estimate),log2(dindex_micro$estimate),log2(dindex_meta$estimate))

r.lower <- c( log2(dindex_tcga$lower), log2(dindex_pcsi$lower), log2(dindex_kirby$lower),
              log2(dindex_icgc_array$lower), log2(dindex_unc$lower),  log2(dindex_chen$lower), 
              log2(dindex_ouh$lower),log2(dindex_zhang$lower), log2(dindex_winter$lower),  log2(dindex_collisson$lower), 
              log2(dindex_seq_lower),log2(dindex_micro_lower), log2(dindex_meta_lower))

r.upper <- c( log2(dindex_tcga$upper), log2(dindex_pcsi$upper), log2(dindex_kirby$upper),
              log2(dindex_icgc_array$upper),  log2(dindex_unc$upper),  log2(dindex_chen$upper),
              log2(dindex_ouh$upper), log2(dindex_zhang$upper), log2(dindex_winter$upper),log2(dindex_collisson$upper), 
              log2(dindex_seq_upper), log2(dindex_micro_upper),log2(dindex_meta_upper))

r.pval <- round(c(dindex_tcga$p.value, dindex_pcsi$p.value, dindex_kirby$p.value,
                  dindex_icgc_array$p.value,dindex_unc$p.value, dindex_chen$p.value,
                  dindex_ouh$p.value, dindex_zhang$p.value, dindex_winter$p.value,dindex_collisson$p.value,   
                  dindex_seq_pval,dindex_micro_pval, dindex_meta_pval),2)

r.pval1 <- c(sprintf("%.1E", dindex_tcga$p.value),sprintf("%.1E", dindex_pcsi$p.value), sprintf("%.1E", dindex_kirby$p.value),
             sprintf("%.1E", dindex_icgc_array$p.value),sprintf("%.1E", dindex_unc$p.value),  sprintf("%.1E", dindex_chen$p.value), 
             sprintf("%.1E", dindex_ouh$p.value), sprintf("%.1E", dindex_zhang$p.value),  sprintf("%.1E", dindex_winter$p.value), 
             sprintf("%.1E", dindex_collisson$p.value), 
             sprintf("%.1E", dindex_seq_pval),sprintf("%.1E", dindex_micro_pval),sprintf("%.1E", dindex_meta_pval))

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","PCSI","Kirby","ICGC-array","UNC","Chen","OUH","Zhang","Winter","Collisson","Sequencing","Microarray","Overall")

data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -14L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts",rownames(t)),
  c("P values",r.pval1))
pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Figures/New_Figure_DINDEX1.pdf")
length_seq= length(tcga_list[[1]]) + length(pcsi_list[[1]])+ length(kirby_list[[1]])
length_micro = length(icgc_array_list[[1]])+length(unc_list[[1]])+ length(chen_list[[1]])+
  length(ouh_list[[1]])+length(zhang_list[[1]])+length(winter_list[[1]])+length(collisson_list[[1]])
length_c=  c(NA,c(length(tcga_list[[1]]), length(pcsi_list[[1]]),length(kirby_list[[1]]),
                  length(icgc_array_list[[1]]),length(unc_list[[1]]), length(chen_list[[1]]),
                  length(ouh_list[[1]]),length(zhang_list[[1]]),length(winter_list[[1]]),length(collisson_list[[1]]),
                  length_seq, length_micro, length_seq +length_micro)/1000 )
fn <- local({
  i = 0
  
  b_clrs =  c(c("#FF7F00","#FF7F00","#FF7F00"),c("#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4"))
  l_clrs =  c(c("#FF7F00","#FF7F00","#FF7F00"),c("#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4"))
  #s_clrs =c(rep("red",10),"green","pink","yellow","orange")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    #fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fn1 <- local({
  i = 0
  
  s_clrs =c("#FF7F00","#1F78B4","grey57")
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})




forestplot(tabletext2,data2,xlab="Log2 HR",is.summary=c(TRUE,rep(FALSE,10),TRUE, TRUE, TRUE),clip=c(-1,2.5),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  xlab  = gpar(fontfamily = "Helvetica", cex = 1)),
           col = fpColors(text="black"),title=" ",new_page = FALSE,
           fn.ci_norm = fn,  fn.ci_sum = fn1, zero=0,graphwidth=unit(2, "inches"),  align=c("l"), pch=16,boxsize = length_c+0.2)
dev.off()

##### Plotting Forestplot of Concordance index ###############

r.mean <- c(con_tcga$c.index,  con_pcsi$c.index, con_kirby$c.index,  con_icgc_array$c.index,  
            con_unc$c.index, con_chen$c.index, con_ouh$c.index, con_zhang$c.index, 
            con_winter$c.index,  con_collisson$c.index,  con_seq$estimate,con_micro$estimate,con_meta$estimate)

r.lower <- c( con_tcga$lower, con_pcsi$lower,   con_kirby$lower, con_icgc_array$lower,
              con_unc$lower,con_chen$lower,con_ouh$lower, con_zhang$lower,   con_winter$lower,
              con_collisson$lower, con_seq_lower,con_micro_lower, con_meta_lower)

r.upper <- c( con_tcga$upper,  con_pcsi$upper, con_kirby$upper,   con_icgc_array$upper,
              con_unc$upper,  con_chen$upper,con_ouh$upper,con_zhang$upper,  
              con_winter$upper,   con_collisson$upper,con_seq_upper, con_micro_upper, con_meta_upper)

r.pval <- round(c( con_tcga$p.value,  con_pcsi$p.value, con_kirby$p.value, con_icgc_array$p.value, 
                   con_unc$p.value,  con_chen$p.value, con_ouh$p.value, con_zhang$p.value,
                   con_winter$p.value,con_collisson$p.value,  con_seq_pval,con_micro_pval,con_meta_pval),2)

r.pval1 <- c(sprintf("%.1E", con_tcga$p.value), sprintf("%.1E", con_pcsi$p.value), sprintf("%.1E", con_kirby$p.value), 
             sprintf("%.1E", con_icgc_array$p.value), sprintf("%.1E", con_unc$p.value), sprintf("%.1E",  con_chen$p.value),
             sprintf("%.1E", con_ouh$p.value), sprintf("%.1E", con_zhang$p.value), sprintf("%.1E", con_winter$p.value),  
             sprintf("%.1E", con_collisson$p.value),sprintf("%.1E", con_seq_pval),sprintf("%.1E", con_micro_pval),
             sprintf("%.1E", con_meta_pval))

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","PCSI","Kirby", "ICGC-array","UNC","Chen","OUH", "Zhang","Winter","Collisson","Sequencing","Microarray","Overall")


data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -14L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts",rownames(t)),
  c("P values",r.pval1))
pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Figures/New_Figure_CONCORDANCE1.pdf")
fn <- local({
  i = 0
  
  b_clrs =  c(c("#FF7F00","#FF7F00","#FF7F00"),c("#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4"))
  l_clrs =  c(c("#FF7F00","#FF7F00","#FF7F00"),c("#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4","#1F78B4"))
  #s_clrs =c(rep("red",10),"green","pink","yellow","orange")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    #fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fn1 <- local({
  i = 0
  
  s_clrs =c("#FF7F00","#1F78B4","grey57")
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

forestplot(tabletext2,data2,xlab="C-index",is.summary=c(TRUE,rep(FALSE,10),TRUE, TRUE, TRUE), clip=c(0.3,0.8), 
           txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.9),  xlab  = gpar(fontfamily = "Helvetica", cex = 1)), 
           col = fpColors(text="black"), fn.ci_norm = fn,  fn.ci_sum = fn1,title="",zero=0.5,graphwidth=unit(2, "inches"),align=c("l"),new_page = FALSE,
           boxsize = length_c+0.2)
dev.off()

#txt_gp =  fpTxtGp(label = gpar(fontfamily = "Verdana"),ticks = gpar(cex=0.8)), 

#####




###################### Plotting ROC CURVES

pcsi_roc=roc(pcsi_grp,pcsi_list[[1]][g_pcsi])$auc[1]
pcsi_roc_se=reportROC(pcsi_grp,pcsi_list[[1]][g_pcsi], plot = FALSE)$AUC.SE

tcga_roc=roc(tcga_grp,tcga_list[[1]][g_tcga])$auc[1]
tcga_roc_se=reportROC(tcga_grp,tcga_list[[1]][g_tcga], plot = FALSE)$AUC.SE

unc_roc=roc(unc_grp,unc_list[[1]][g_unc])$auc[1] 
unc_roc_se=reportROC(unc_grp,unc_list[[1]][g_unc], plot = FALSE)$AUC.SE

zhang_roc=roc(zhang_grp,zhang_list[[1]][g_zhang])$auc[1]
zhang_roc_se=reportROC(zhang_grp,zhang_list[[1]][g_zhang], plot = FALSE)$AUC.SE

winter_roc=roc(winter_grp,winter_list[[1]][g_winter])$auc[1]
winter_roc_se=reportROC(winter_grp,winter_list[[1]][g_winter], plot = FALSE)$AUC.SE

ouh_roc=roc(ouh_grp,ouh_list[[1]][g_ouh])$auc[1]
ouh_roc_se=reportROC(ouh_grp,ouh_list[[1]][g_ouh], plot = FALSE)$AUC.SE

icgc_array_roc=roc(icgc_array_grp,icgc_array_list[[1]][g_icgc_arr])$auc[1]
icgc_array_roc_se=reportROC(icgc_array_grp,icgc_array_list[[1]][g_icgc_arr], plot = FALSE)$AUC.SE

collisson_roc=roc(collisson_grp,collisson_list[[1]][g_coll])$auc[1]
collisson_roc_se=reportROC(collisson_grp,collisson_list[[1]][g_coll], plot = FALSE)$AUC.SE

chen_roc=roc(chen_grp,chen_list[[1]][g_chen])$auc[1]
chen_roc_se=reportROC(chen_grp,chen_list[[1]][g_chen], plot = FALSE)$AUC.SE

kirby_roc=roc(kirby_grp,kirby_list[[1]][g_kirby])$auc[1]
kirby_roc_se=reportROC(kirby_grp,kirby_list[[1]][g_kirby], plot = FALSE)$AUC.SE

meta_auc = combine.est(c(pcsi_roc, tcga_roc, kirby_roc, unc_roc,
                         zhang_roc, winter_roc, ouh_roc,icgc_array_roc, 
                         collisson_roc, chen_roc),
                       c( pcsi_roc_se, tcga_roc_se,kirby_roc_se, unc_roc_se,zhang_roc_se, 
                          winter_roc_se, ouh_roc_se,icgc_array_roc_se, collisson_roc_se, 
                          chen_roc_se),na.rm=TRUE,hetero = TRUE)$estimate

seq_auc = combine.est(c(pcsi_roc, tcga_roc, kirby_roc), 
                      c( pcsi_roc_se, tcga_roc_se,kirby_roc_se),na.rm=TRUE,hetero =TRUE)$estimate

micro_auc = combine.est(c( unc_roc, zhang_roc, winter_roc, ouh_roc,icgc_array_roc,collisson_roc,  chen_roc), 
                        c(  unc_roc_se,zhang_roc_se, winter_roc_se, ouh_roc_se,icgc_array_roc_se, collisson_roc_se,  chen_roc_se),
                        na.rm=TRUE,hetero = TRUE)$estimate


pcsi_pval = roc.area(pcsi_grp,1-pcsi_list[[1]][g_pcsi])$p.value  
tcga_pval = roc.area(tcga_grp,1-tcga_list[[1]][g_tcga])$p.value  
ouh_pval = roc.area(ouh_grp,1-ouh_list[[1]][g_ouh])$p.value  
unc_pval =roc.area(unc_grp,1-unc_list[[1]][g_unc])$p.value  
zhang_pval =roc.area(zhang_grp,1-zhang_list[[1]][g_zhang])$p.value
winter_pval =roc.area(winter_grp,1-winter_list[[1]][g_winter])$p.value
icgc_array_pval = roc.area(icgc_array_grp,1-icgc_array_list[[1]][g_icgc_arr])$p.value  
collisson_pval=roc.area(collisson_grp,1-collisson_list[[1]][g_coll])$p.value
chen_pval=roc.area(chen_grp,1-chen_list[[1]][g_chen])$p.value
kirby_pval=roc.area(kirby_grp,1-kirby_list[[1]][g_kirby])$p.value

m1= combine.est(c(pcsi_roc, tcga_roc, kirby_roc, unc_roc,
                  zhang_roc, winter_roc, ouh_roc,icgc_array_roc, 
                  collisson_roc, chen_roc),
                c( pcsi_roc_se, tcga_roc_se,kirby_roc_se, unc_roc_se,zhang_roc_se, 
                   winter_roc_se, ouh_roc_se,icgc_array_roc_se, collisson_roc_se, 
                   chen_roc_se),na.rm=TRUE,hetero = TRUE)
m2=combine.est(c(pcsi_roc, tcga_roc, kirby_roc), 
               c( pcsi_roc_se, tcga_roc_se,kirby_roc_se),na.rm=TRUE,hetero =TRUE)
m3= combine.est(c( unc_roc, zhang_roc, winter_roc, ouh_roc,icgc_array_roc,collisson_roc,  chen_roc), 
                c(  unc_roc_se,zhang_roc_se, winter_roc_se, ouh_roc_se,icgc_array_roc_se, collisson_roc_se,  chen_roc_se),
                na.rm=TRUE,hetero = TRUE)

meta_pval= pnorm((m1$estimate - 0.5)/m1$se, lower.tail = m1$estimate < 0.5) * 2
seq_pval=pnorm((m2$estimate -0.5)/m2$se, lower.tail = m2$estimate < 0.5) * 2
micro_pval=pnorm((m3$estimate -0.5)/m3$se, lower.tail = m3$estimate < 0.5) * 2
# 
# meta_pval = combine.test(p=c(pcsi_pval,tcga_pval, ouh_pval, unc_pval, zhang_pval,winter_pval,icgc_array_pval,
#                              collisson_pval,chen_pval,kirby_pval ),
#                          w= c( length(pcsi_list[[1]][g_icgc_seq]), length(tcga_list[[1]][g_tcga]),  length(ouh_list[[1]][g_ouh]),
#                                length(unc_list[[1]][g_unc]), length(zhang_list[[1]][g_zhang]), length(winter_list[[1]][g_winter]),  
#                                length(icgc_array_list[[1]][g_icgc_arr]), length(collisson_list[[1]][g_coll]), length(chen_list[[1]][g_chen]),
#                                length(kirby_list[[1]][g_kirby])),na.rm=TRUE,method="z.transform")
# 
# seq_pval = combine.test(p=c(pcsi_pval,tcga_pval, kirby_pval), w= c( length(pcsi_list[[1]][g_icgc_seq]), length(tcga_list[[1]][g_tcga]), 
#                                                                     length(kirby_list[[1]][g_kirby])),na.rm=TRUE)
# 
# micro_pval = combine.test(p=c(ouh_pval, unc_pval, zhang_pval,winter_pval,icgc_array_pval,collisson_pval,chen_pval),
#                           w= c( length(ouh_list[[1]][g_ouh]),length(unc_list[[1]][g_unc]), length(zhang_list[[1]][g_zhang]), 
#                                 length(winter_list[[1]][g_winter]), length(icgc_array_list[[1]][g_icgc_arr]), 
#                                 length(collisson_list[[1]][g_coll]),length(chen_list[[1]][g_chen])),na.rm=TRUE,method="z.transform")


library(pROC)


pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Figures/New_Figure_AUC.pdf")
plot(roc(tcga_grp,tcga_list[[1]][g_tcga]),lwd=4, col="chartreuse3",lty=1)
plot(roc(kirby_grp,kirby_list[[1]][g_kirby]),lwd=4, col="magenta",add=TRUE,lty=1)
plot(roc(pcsi_grp,pcsi_list[[1]][g_pcsi]),lwd=4, col="#fb9a99",add=TRUE,lty=1)
plot(roc(icgc_array_grp,icgc_array_list[[1]][g_icgc_arr]),lwd=4, col="turquoise3",add=TRUE,lty=3)
plot(roc(unc_grp,unc_list[[1]][g_unc]),lwd=4, col="darkgoldenrod1",add=TRUE,lty=3)
plot(roc(zhang_grp,zhang_list[[1]][g_zhang]),lwd=4, col="wheat4",add=TRUE,lty=3)
plot(roc(chen_grp,chen_list[[1]][g_chen]),lwd=4, col="green",add=TRUE,lty=3)
plot(roc(collisson_grp,collisson_list[[1]][g_coll]),lwd=4, col="red",add=TRUE,lty=3)
plot(roc(winter_grp,winter_list[[1]][g_winter]),lwd=4, col="cornflowerblue",add=TRUE,lty=3)
plot(roc(ouh_grp,ouh_list[[1]][g_ouh]),lwd=4, col="mediumorchid2",add=TRUE,lty=3)

legend("bottomright",legend=c(paste("TCGA: ",round(tcga_roc,digits=2)," (P = ", sprintf("%.1E", tcga_pval),")", sep=""),
                              paste("Kirby: ",round(kirby_roc,digits=2)," (P = ", sprintf("%.1E", kirby_pval),")", sep=""),
                              paste("PCSI: ",round(pcsi_roc,digits=2)," (P = ", sprintf("%.1E", pcsi_pval),")", sep=""),
                              paste("ICGC-array: ",round(icgc_array_roc,digits=2)," (P = ", sprintf("%.1E", icgc_array_pval),")", sep=""),
                              paste("UNC: ",round(unc_roc,digits=2)," (P = ", sprintf("%.1E", unc_pval),")", sep=""),
                              paste("Zhang: ",round(zhang_roc,digits=2)," (P = ", sprintf("%.1E", zhang_pval),")", sep=""),
                              paste("Chen: ",round(chen_roc,digits=2)," (P = ", sprintf("%.1E", chen_pval),")", sep=""),
                              paste("Collisson: ",round(collisson_roc,digits=2)," (P = ", sprintf("%.1E", collisson_pval),")", sep=""),
                              paste("Winter: ",round(winter_roc,digits=2)," (P = ", sprintf("%.1E", winter_pval),")", sep=""),
                              paste("OUH: ",round(ouh_roc,digits=2)," (P = ", sprintf("%.1E", ouh_pval),")", sep="")),
       
       
       fill=c( "chartreuse3","magenta","#fb9a99","turquoise3","darkgoldenrod1","wheat4","green","red","cornflowerblue","mediumorchid2"),y.intersp = 1, cex=0.9,bty = "n")

dev.off()

