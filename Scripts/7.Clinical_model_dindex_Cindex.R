
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/Clinical_models1.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/clinical_features_included_censored_As_well_new_PROB.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/RData/cohorts_subtypes_Average_data.RData")

######3 Including  censored data for C-index and D-INDEX caluculation and comprasion across the models.
model1=model1

pcsi=data.frame(clinical_features$pcsi_clinical)
icgc=data.frame(clinical_features$icgc_clinical)
tcga=data.frame(clinical_features$tcga_clinical)
icgc_arr=data.frame(clinical_features$icgc_arr_clinical)
ouh=data.frame(clinical_features$ouh_clinical)

pcsi$Age=as.numeric(as.character(pcsi$Age))
pcsi$Sex = as.character(pcsi$Sex)
pcsi$T_status=gsub(" ", "", pcsi$T_status, fixed = TRUE)
pcsi$pred_prob=as.numeric(as.character(pcsi$pred_prob))

pcsi_cl_pred1=1-predict(model1,pcsi, na.action = na.exclude,type="response")
pcsi_cl_pred2= pcsi$pred_prob[which(pcsi$ID %in% names(pcsi_cl_pred1))]


######### Validation cohorts including censored samples as well
########## ICGC

icgc$Age=as.numeric(as.character(icgc$Age))
icgc$Sex = as.character(icgc$Sex)
icgc$T_status = as.character(icgc$T_status)

icgc_cl_pred1=1-predict(model1, icgc, na.action = na.exclude,type="response")
icgc_cl_pred2= icgc$pred_prob[which(icgc$ID %in% names(icgc_cl_pred1))]

########## TCGA

tcga$Age=as.numeric(as.character(tcga$Age))
tcga$Sex = as.character(tcga$Sex)
tcga$T_status = as.character(tcga$T_status)
tcga$pred_prob=as.numeric(as.character(tcga$pred_prob))

tcga_cl_pred1=1-predict(model1, tcga, na.action = na.exclude,type="response")
tcga_cl_pred2= tcga$pred_prob[which(tcga$ID %in% names(tcga_cl_pred1))]

########## ICGC array

icgc_arr$Age=as.numeric(as.character(icgc_arr$Age))
icgc_arr$Sex = as.character(icgc_arr$Sex)
icgc_arr$T_status = as.character(icgc_arr$T_status)
icgc_arr$pred_prob=as.numeric(as.character(icgc_arr$pred_prob))

icgc_arr_cl_pred1=1-predict(model1, icgc_arr, na.action = na.exclude,type="response")
icgc_arr_cl_pred2= icgc_arr$pred_prob[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))]

########## OUH

ouh$Age=as.numeric(as.character(ouh$Age))
ouh$Sex = as.character(ouh$Sex)

ouh_cl_pred1 = 1-predict(model1, ouh, na.action = na.exclude,type="response")
ouh_cl_pred2= ouh$pred_prob

###########################################################################
###########################################################################
###########################################################################
### Concordance indices and Dindex calcualtion 

## Clinical model
dindex_ouh <- D.index(x=ouh_cl_pred1, surv.time=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS)), surv.event=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_ouh <- concordance.index(x=ouh_cl_pred1, surv.time=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS)), surv.event=as.numeric(as.character(ouh[names(ouh_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_pcsi <- D.index(x=pcsi_cl_pred1, surv.time=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS)), surv.event=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_pcsi <- concordance.index(x=pcsi_cl_pred1, surv.time=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS)), surv.event=as.numeric(as.character(pcsi[names(pcsi_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_tcga <- D.index(x=tcga_cl_pred1, surv.time=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS)), surv.event=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_tcga <- concordance.index(x=tcga_cl_pred1, surv.time=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS)), surv.event=as.numeric(as.character(tcga[names(tcga_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");
dindex_icgc_array <- D.index(x=icgc_arr_cl_pred1 , surv.time=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS)), surv.event=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_icgc_array <- concordance.index(x=icgc_arr_cl_pred1 , surv.time=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS)), surv.event=as.numeric(as.character(icgc_arr[names(icgc_arr_cl_pred1),]$OS_Status)), na.rm=TRUE, method="noether");

##PCOSP
dindex_ouh1 <- D.index(x=ouh_cl_pred2, surv.time=as.numeric(as.character(ouh$OS)), surv.event=as.numeric(as.character(ouh$OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_ouh1 <- concordance.index(x=ouh_cl_pred2, surv.time=as.numeric(as.character(ouh$OS)), surv.event=as.numeric(as.character(ouh$OS_Status)), na.rm=TRUE, method="noether");
dindex_pcsi1 <- D.index(x=pcsi_cl_pred2, surv.time=as.numeric(as.character(pcsi$OS[which(pcsi$ID %in% names(pcsi_cl_pred1))])), surv.event=as.numeric(as.character(pcsi$OS_Status[which(pcsi$ID %in% names(pcsi_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_pcsi1 <- concordance.index(x=pcsi_cl_pred2, surv.time=as.numeric(as.character(pcsi$OS[which(pcsi$ID %in% names(pcsi_cl_pred1))])), surv.event=as.numeric(as.character(pcsi$OS_Status[which(pcsi$ID %in% names(pcsi_cl_pred1))])), na.rm=TRUE, method="noether");
dindex_tcga1 <- D.index(x=tcga_cl_pred2, surv.time=as.numeric(as.character(tcga$OS[which(tcga$ID %in% names(tcga_cl_pred1))])), surv.event=as.numeric(as.character(tcga$OS_Status[which(tcga$ID %in% names(tcga_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_tcga1 <- concordance.index(x=tcga_cl_pred2, surv.time=as.numeric(as.character(tcga$OS[which(tcga$ID %in% names(tcga_cl_pred1))])), surv.event=as.numeric(as.character(tcga$OS_Status[which(tcga$ID %in% names(tcga_cl_pred1))])), na.rm=TRUE, method="noether");
dindex_icgc_array1 <- D.index(x=icgc_arr_cl_pred2 , surv.time=as.numeric(as.character(icgc_arr$OS[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), surv.event=as.numeric(as.character(icgc_arr$OS_Status[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
con_icgc_array1 <- concordance.index(x=icgc_arr_cl_pred2 , surv.time=as.numeric(as.character(icgc_arr$OS[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), surv.event=as.numeric(as.character(icgc_arr$OS_Status[which(icgc_arr$ID %in% names(icgc_arr_cl_pred1))])), na.rm=TRUE, method="noether");


####################################### Meta estimates calculations
## Clinical model
dindex_meta <- combine.est(c( dindex_pcsi$d.index, dindex_tcga$d.index, dindex_ouh$d.index,dindex_icgc_array$d.index),
                           c(dindex_pcsi$se, dindex_tcga$se, dindex_ouh$se,dindex_icgc_array$se),na.rm = TRUE,hetero =  TRUE)

con_meta <- combine.est(c( con_pcsi$c.index, con_tcga$c.index, con_ouh$c.index,con_icgc_array$c.index),
                        c(con_pcsi$se, con_tcga$se, con_ouh$se,con_icgc_array$se),na.rm = TRUE,hetero =  TRUE)

dindex_seq <- combine.est(c( dindex_pcsi$d.index, dindex_tcga$d.index),c(dindex_pcsi$se, dindex_tcga$se),na.rm = TRUE)
con_seq <- combine.est(c( con_pcsi$c.index, con_tcga$c.index),c(con_pcsi$se, con_tcga$se),na.rm = TRUE)

dindex_micro <- combine.est(c( dindex_ouh$d.index,dindex_icgc_array$d.index),c( dindex_ouh$se,dindex_icgc_array$se),na.rm = TRUE)
con_micro <- combine.est(c(  con_ouh$c.index,con_icgc_array$c.index),c( con_ouh$se,con_icgc_array$se),na.rm = TRUE)

dindex_meta_lower <- dindex_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta$se
con_meta_lower <- con_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta$se
dindex_seq_lower <- dindex_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq$se
con_seq_lower <- con_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq$se
dindex_micro_lower <- dindex_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro$se
con_micro_lower <- con_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro$se
dindex_meta_upper <- dindex_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta$se
con_meta_upper <- con_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta$se
dindex_seq_upper <- dindex_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq$se
con_seq_upper <- con_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq$se
dindex_micro_upper <- dindex_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro$se
con_micro_upper <- con_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro$se
# 
# dindex_meta_pval <- combine.test(p=c(dindex_pcsi$p.value, dindex_tcga$p.value, dindex_ouh$p.value,dindex_icgc_array$p.value),
#                                  w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1),length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
# dindex_micro_pval <- combine.test(p=c(dindex_ouh$p.value,dindex_icgc_array$p.value),
#                                   w=c(length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
# dindex_seq_pval <- combine.test(p=c(dindex_pcsi$p.value, dindex_tcga$p.value),
#                                 w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1)),hetero = FALSE,method="z.transform")
# con_meta_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value, con_ouh$p.value,con_icgc_array$p.value),
#                               w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1),length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
# con_micro_pval <- combine.test(p=c(con_ouh$p.value,con_icgc_array$p.value),w=c(length(ouh_cl_pred1),length(icgc_arr_cl_pred1)),hetero = FALSE,method="z.transform")
# con_seq_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value),w=c(length(pcsi_cl_pred1), length(tcga_cl_pred1)),hetero = FALSE,method="z.transform")


dindex_meta_pval <- 2*pnorm(-abs(log(dindex_meta$estimate)/dindex_meta$se))
#pnorm((dindex_meta$estimate -0.5)/dindex_meta$se, lower.tail = dindex_meta$estimate < 0.5) * 2
dindex_micro_pval <- 2*pnorm(-abs(log(dindex_micro$estimate)/dindex_micro$se))
dindex_seq_pval <- 2*pnorm(-abs(log(dindex_seq$estimate)/dindex_seq$se))

con_meta_pval <- pnorm((con_meta$estimate -0.5)/con_meta$se, lower.tail = con_meta$estimate < 0.5) * 2
con_micro_pval <- pnorm((con_micro$estimate -0.5)/con_micro$se, lower.tail = con_micro$estimate < 0.5) * 2
con_seq_pval <- pnorm((con_seq$estimate -0.5)/con_seq$se, lower.tail = con_seq$estimate < 0.5) * 2





## PCOSP
dindex_meta1 <- combine.est(c( dindex_pcsi1$d.index, dindex_tcga1$d.index, dindex_ouh1$d.index,dindex_icgc_array1$d.index),c(dindex_pcsi1$se, dindex_tcga1$se, dindex_ouh1$se,dindex_icgc_array1$se),na.rm = TRUE, hetero =  TRUE)
con_meta1 <- combine.est(c(con_pcsi1$c.index, con_tcga1$c.index, con_ouh1$c.index, con_icgc_array1$c.index),c(con_pcsi1$se, con_tcga1$se, con_ouh1$se,con_icgc_array1$se),na.rm = TRUE,hetero =  TRUE)
dindex_seq1 <- combine.est(c( dindex_pcsi1$d.index, dindex_tcga1$d.index),c(dindex_pcsi1$se, dindex_tcga1$se),na.rm = TRUE)
con_seq1 <- combine.est(c( con_pcsi1$c.index, con_tcga1$c.index),c(con_pcsi1$se, con_tcga1$se),na.rm = TRUE)
dindex_micro1 <- combine.est(c( dindex_ouh1$d.index,dindex_icgc_array1$d.index),c( dindex_ouh1$se,dindex_icgc_array1$se),na.rm = TRUE)
con_micro1 <- combine.est(c(  con_ouh1$c.index,con_icgc_array1$c.index),c( con_ouh1$se,con_icgc_array1$se),na.rm = TRUE)
dindex_meta_lower1 <- dindex_meta1$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta1$se
con_meta_lower1 <- con_meta1$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta1$se
dindex_seq_lower1 <- dindex_seq1$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq1$se
con_seq_lower1 <- con_seq1$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq1$se
dindex_micro_lower1 <- dindex_micro1$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro1$se
con_micro_lower1 <- con_micro1$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro1$se
dindex_meta_upper1 <- dindex_meta1$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta1$se
con_meta_upper1 <- con_meta1$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta1$se
dindex_seq_upper1 <- dindex_seq1$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq1$se
con_seq_upper1 <- con_seq1$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq1$se
dindex_micro_upper1 <- dindex_micro1$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro1$se
con_micro_upper1 <- con_micro1$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro1$se
# 
# dindex_meta_pval1 <- combine.test(p=c(dindex_pcsi1$p.value, dindex_tcga1$p.value, dindex_ouh1$p.value,dindex_icgc_array1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2),length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
# dindex_micro_pval1 <- combine.test(p=c(dindex_ouh1$p.value,dindex_icgc_array1$p.value),w=c(length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
# dindex_seq_pval1 <- combine.test(p=c(dindex_pcsi1$p.value, dindex_tcga1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2)),hetero = FALSE,method="z.transform")
# con_meta_pval1 <- combine.test(p=c(con_pcsi1$p.value, con_tcga1$p.value, con_ouh1$p.value,con_icgc_array1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2),length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
# con_micro_pval1 <- combine.test(p=c(con_ouh1$p.value,con_icgc_array1$p.value),w=c(length(ouh_cl_pred2),length(icgc_arr_cl_pred2)),hetero = FALSE,method="z.transform")
# con_seq_pval1 <- combine.test(p=c(con_pcsi1$p.value, con_tcga1$p.value),w=c(length(pcsi_cl_pred2), length(tcga_cl_pred2)),hetero = FALSE,method="z.transform")

dindex_meta_pval1 <- 2*pnorm(-abs(log(dindex_meta1$estimate)/dindex_meta1$se))
#pnorm((dindex_meta$estimate -0.5)/dindex_meta$se, lower.tail = dindex_meta$estimate < 0.5) * 2
dindex_micro_pval1 <- 2*pnorm(-abs(log(dindex_micro1$estimate)/dindex_micro1$se))
dindex_seq_pval1 <- 2*pnorm(-abs(log(dindex_seq1$estimate)/dindex_seq1$se))

con_meta_pval1 <- pnorm((con_meta1$estimate -0.5)/con_meta1$se, lower.tail = con_meta1$estimate < 0.5) * 2
con_micro_pval1 <- pnorm((con_micro1$estimate -0.5)/con_micro1$se, lower.tail = con_micro1$estimate < 0.5) * 2
con_seq_pval1 <- pnorm((con_seq1$estimate -0.5)/con_seq1$se, lower.tail = con_seq1$estimate < 0.5) * 2


############# Plotting Concordance index and comparison across models
#pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Figures/Clinical_c_INDEX1.pdf")
r.mean <- c( NA,con_tcga$c.index,con_tcga1$c.index, NA,NA,
             con_pcsi$c.index,  con_pcsi1$c.index, NA,NA,
             con_icgc_array$c.index,con_icgc_array1$c.index,NA,NA,
             con_ouh$c.index,con_ouh1$c.index, NA,NA,
             con_seq$estimate,con_seq1$estimate,NA,NA,
           con_micro$estimate,con_micro1$estimate, NA,NA,
           con_meta$estimate,con_meta1$estimate
            )
           
r.lower<- c( NA,con_tcga$lower,con_tcga1$lower, NA,NA,
             con_pcsi$lower,  con_pcsi1$lower, NA,NA,
             con_icgc_array$lower,con_icgc_array1$lower,NA,NA,
             con_ouh$lower,con_ouh1$lower, NA,NA,
             con_seq_lower,con_seq_lower1,NA,NA,
             con_micro_lower,con_micro_lower1, NA,NA,
             con_meta_lower,con_meta_lower1
)
r.upper<- c(  NA,con_tcga$upper,con_tcga1$upper, NA,NA,
             con_pcsi$upper,  con_pcsi1$upper, NA,NA,
             con_icgc_array$upper,con_icgc_array1$upper,NA,NA,
             con_ouh$upper,con_ouh1$upper, NA,NA,
             con_seq_upper,con_seq_upper1,NA,NA,
             con_micro_upper,con_micro_upper1, NA,NA,
             con_meta_upper,con_meta_upper1
)


r.pval<- c( NA,con_tcga$p.value,con_tcga1$p.value, NA,NA,
            con_pcsi$p.value,  con_pcsi1$p.value, NA,NA,
            con_icgc_array$p.value,con_icgc_array1$p.value,NA,NA,
            con_ouh$p.value,con_ouh1$p.value, NA,NA,
            con_seq_pval,con_seq_pval1,NA,NA,
            con_micro_pval,con_micro_pval1, NA,NA,
            con_meta_pval,con_meta_pval1
)

r.pval1<- c( NA,sprintf("%.1E", con_tcga$p.value) ,sprintf("%.1E", con_tcga1$p.value), NA,NA,
             sprintf("%.1E", con_pcsi$p.value),  sprintf("%.1E", con_pcsi1$p.value), NA,NA,
             sprintf("%.1E", con_icgc_array$p.value), sprintf("%.1E", con_icgc_array1$p.value),NA,NA,
             sprintf("%.1E", con_ouh$p.value), sprintf("%.1E", con_ouh1$p.value), NA,NA,
             sprintf("%.1E", con_seq_pval) ,sprintf("%.1E", con_seq_pval1),NA,NA,
             sprintf("%.1E", con_micro_pval) ,sprintf("%.1E", con_micro_pval1), NA,NA,
             sprintf("%.1E", con_meta_pval) ,sprintf("%.1E", con_meta_pval1)
)



#pdf("/Users/vandanasandhu/Desktop/c.pdf")
t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","Clinical model", "PCOSP", "",
                  "PCSI","Clinical model", "PCOSP", "",
                  "ICGC-array", "Clinical model", "PCOSP", "",
                  "OUH", "Clinical model", "PCOSP", "",
                  "Sequencing", "Clinical model", "PCOSP", "",
                  "Microarray","Clinical model", "PCOSP","", 
                  "Overall", "Clinical model", "PCOSP" )
data2 <-  
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA,  -28L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts",rownames(t)),
  c("P values",r.pval1))

seq=length(tcga_cl_pred1) +  length(pcsi_cl_pred1)
arr=length(icgc_arr_cl_pred1)+ length(ouh_cl_pred1)

#pdf("/Users/vandanasandhu/Desktop/c.pdf")

fn <- local({
  i = 0
  
  b_clrs =  c("palevioletred1","darkgrey","palevioletred1","darkgrey", "palevioletred1","darkgrey","palevioletred1","darkgrey")
  l_clrs =    c("palevioletred1","darkgrey","palevioletred1","darkgrey", "palevioletred1","darkgrey","palevioletred1","darkgrey")
  #s_clrs =c(rep("palevioletred1",10),"green","pink","darkgrey","orange")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    #fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fn1 <- local({
  i = 0
  
  s_clrs =c("palevioletred1","darkgrey","palevioletred1","darkgrey","black","black")
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})


forestplot(tabletext2,data2,xlab="C-index",is.summary=c(TRUE, TRUE, FALSE, FALSE, 
                                                                  TRUE, TRUE, FALSE, FALSE, 
                                                                  TRUE, TRUE, FALSE, FALSE, 
                                                                  TRUE, TRUE, FALSE, FALSE, 
                                                                  TRUE, TRUE,TRUE, TRUE,
                                                                  TRUE, TRUE, TRUE, TRUE,
                                                                  TRUE, TRUE, TRUE, TRUE), clip=c(0,3.0),cex=10,
           fn.ci_norm = fn,  fn.ci_sum = fn1,zero=0.5,graphwidth=unit(2, "inches"), align=c("l"), new_page = FALSE,txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  
            xlab  = gpar(fontfamily = "Helvetica", cex = 1)), col = fpColors(text="black"))
dev.off()

############# Plotting D-index and comparison across models

pdf("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Figures/Clinical_dindex.pdf")

r.mean <- c( NA,log2(dindex_tcga$d.index),log2(dindex_tcga1$d.index), NA,NA,
             log2(dindex_pcsi$d.index),  log2(dindex_pcsi1$d.index), NA,NA,
             log2(dindex_icgc_array$d.index),log2(dindex_icgc_array1$d.index),NA,NA,
             log2(dindex_ouh$d.index),log2(dindex_ouh1$d.index), NA,NA,
             log2(dindex_seq$estimate),log2(dindex_seq1$estimate),NA,NA,
             log2(dindex_micro$estimate),log2(dindex_micro1$estimate), NA,NA,
             log2(dindex_meta$estimate),log2(dindex_meta1$estimate)
)

r.lower<- c( NA,log2(dindex_tcga$lower),log2(dindex_tcga1$lower), NA,NA,
              log2(dindex_pcsi$lower),  log2(dindex_pcsi1$lower), NA,NA,
              log2(dindex_icgc_array$lower),log2(dindex_icgc_array1$lower),NA,NA,
              log2(dindex_ouh$lower),log2(dindex_ouh1$lower), NA,NA,
              log2(dindex_seq_lower),log2(dindex_seq_lower1),NA,NA,
              log2(dindex_micro_lower),log2(dindex_micro_lower1), NA,NA,
              log2(dindex_meta_lower),log2(dindex_meta_lower1)
)
r.upper<- c(  NA,log2(dindex_tcga$upper),log2(dindex_tcga1$upper), NA,NA,
               log2(dindex_pcsi$upper),  log2(dindex_pcsi1$upper), NA,NA,
               log2(dindex_icgc_array$upper),log2(dindex_icgc_array1$upper),NA,NA,
               log2(dindex_ouh$upper),log2(dindex_ouh1$upper), NA,NA,
               log2(dindex_seq_upper),log2(dindex_seq_upper1),NA,NA,
               log2(dindex_micro_upper),log2(dindex_micro_upper1), NA,NA,
               log2(dindex_meta_upper),log2(dindex_meta_upper1)
)



r.pval<- c( NA,dindex_tcga$p.value,dindex_tcga1$p.value, NA,NA,
            dindex_pcsi$p.value,  dindex_pcsi1$p.value, NA,NA,
            dindex_icgc_array$p.value,dindex_icgc_array1$p.value,NA,NA,
            dindex_ouh$p.value,dindex_ouh1$p.value, NA,NA,
            dindex_seq_pval,dindex_seq_pval1,NA,NA,
            dindex_micro_pval,dindex_micro_pval1, NA,NA,
            dindex_meta_pval,dindex_meta_pval1
)

r.pval1<- c( NA,sprintf("%.1E", dindex_tcga$p.value) ,sprintf("%.1E", dindex_tcga1$p.value), NA,NA,
             sprintf("%.1E", dindex_pcsi$p.value),  sprintf("%.1E", dindex_pcsi1$p.value), NA,NA,
             sprintf("%.1E", dindex_icgc_array$p.value), sprintf("%.1E", dindex_icgc_array1$p.value),NA,NA,
             sprintf("%.1E", dindex_ouh$p.value), sprintf("%.1E", dindex_ouh1$p.value), NA,NA,
             sprintf("%.1E", dindex_seq_pval) ,sprintf("%.1E", dindex_seq_pval1),NA,NA,
             sprintf("%.1E", dindex_micro_pval) ,sprintf("%.1E", dindex_micro_pval1), NA,NA,
             sprintf("%.1E", dindex_meta_pval) ,sprintf("%.1E", dindex_meta_pval1)
)

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("TCGA","Clinical model", "PCOSP", "",
                  "PCSI","Clinical model", "PCOSP", "",
                  "ICGC-array", "Clinical model", "PCOSP", "",
                  "OUH", "Clinical model", "PCOSP", "",
                  "Sequencing", "Clinical model", "PCOSP", "",
                  "Microarray","Clinical model", "PCOSP","", 
                  "Overall", "Clinical model", "PCOSP" )
data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA,  -28L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Cohorts",rownames(t)),
  c("P values",r.pval1))

seq=length(tcga_cl_pred1) +  length(pcsi_cl_pred1)
arr=length(icgc_arr_cl_pred1)+ length(ouh_cl_pred1)

fn <- local({
  i = 0
  
  
  b_clrs =  c("palevioletred1","darkgrey","palevioletred1","darkgrey", "palevioletred1","darkgrey","palevioletred1","darkgrey")
  l_clrs =    c("palevioletred1","darkgrey","palevioletred1","darkgrey", "palevioletred1","darkgrey","palevioletred1","darkgrey")
  #s_clrs =c(rep("red",10),"green","pink","yellow","orange")
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
    #fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

fn1 <- local({
  i = 0
  
  s_clrs =c("palevioletred1","darkgrey","palevioletred1","darkgrey","black","black")
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

       
forestplot(tabletext2,data2,xlab="Log2 HR",is.summary=c(TRUE, TRUE, FALSE, FALSE, 
                                                             TRUE, TRUE, FALSE, FALSE, 
                                                             TRUE, TRUE, FALSE, FALSE, 
                                                             TRUE, TRUE, FALSE, FALSE, 
                                                             TRUE, TRUE,TRUE, TRUE,
                                                             TRUE, TRUE, TRUE, TRUE,
                                                             TRUE, TRUE, TRUE, TRUE),clip=c(-2,3.0),cex=9,
           fn.ci_norm = fn,  fn.ci_sum = fn1,zero=0,graphwidth=unit(2, "inches"), align=c("l"), new_page = FALSE,txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  
                                                                                                                                  xlab  = gpar(fontfamily = "Helvetica", cex = 1)), col = fpColors(text="black"))



dev.off()

###########################

pcosp_clinical_cindex = cindex.comp.meta(list.cindex1 = list(con_pcsi1, con_tcga1, con_ouh1, con_icgc_array1),
                                  list.cindex2 = list(con_pcsi, con_tcga, con_ouh, con_icgc_array))
pcosp_clinical_cindex

pcosp_clinical_dindex = dindex.comp.meta(list.dindex1 = list(dindex_pcsi1, dindex_tcga1, dindex_ouh1, dindex_icgc_array1),
                                  list.dindex2 = list(dindex_pcsi, dindex_tcga, dindex_ouh, dindex_icgc_array))
pcosp_clinical_dindex 
