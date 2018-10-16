
predict_meta_estimates = function(pcsi_list, tcga_list, kirby_list,icgc_array_list,unc_list,winter_list,  collisson_list, ouh_list,  
                                  pcsi_OS, pcsi_OS_Status,  
                                  tcga_OS, tcga_OS_Status,
                                  kirby_OS, kirby_OS_Status,  
                                  icgc_array_OS, icgc_array_OS_Status,  
                                  unc_OS, unc_OS_Status,
                                  winter_OS, winter_OS_Status,  
                                  collisson_OS, collisson_OS_Status,
                                  ouh_OS, ouh_OS_Status ){
  
  dindex_ouh <- D.index(x=ouh_list, surv.time=as.numeric(as.character(ouh_OS)), surv.event=as.numeric(as.character(ouh_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_tcga <- D.index(x=tcga_list, surv.time=as.numeric(as.character(tcga_OS)), surv.event=as.numeric(as.character(tcga_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_winter <- D.index(x=winter_list, surv.time=as.numeric(as.character(winter_OS)), surv.event=as.numeric(as.character(winter_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_icgc_array <- D.index(x=icgc_array_list, surv.time=as.numeric(as.character(icgc_array_OS)), surv.event=as.numeric(as.character(icgc_array_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_unc <- D.index(x=unc_list, surv.time=as.numeric(as.character(unc_OS)), surv.event=as.numeric(as.character(unc_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_collisson <- D.index(x=collisson_list, surv.time=as.numeric(as.character(collisson_OS)), surv.event=as.numeric(as.character(collisson_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_kirby <- D.index(x=kirby_list, surv.time=as.numeric(as.character(kirby_OS)), surv.event=as.numeric(as.character(kirby_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  dindex_pcsi <- D.index(x=pcsi_list, surv.time=as.numeric(as.character(pcsi_OS)), surv.event=as.numeric(as.character(pcsi_OS_Status)), na.rm=TRUE, alpha = 0.05, method.test = "logrank");
  
  
  ########### Concordance index calculation
  con_ouh <- concordance.index(x=ouh_list, surv.time=as.numeric(as.character(ouh_OS)), surv.event=as.numeric(as.character(ouh_OS_Status)), na.rm=TRUE, method="noether");
  con_tcga <- concordance.index(x=tcga_list, surv.time=as.numeric(as.character(tcga_OS)), surv.event=as.numeric(as.character(tcga_OS_Status)), na.rm=TRUE, method="noether");
  con_winter <- concordance.index(x=winter_list, surv.time=as.numeric(as.character(winter_OS)), surv.event=as.numeric(as.character(winter_OS_Status)), na.rm=TRUE, method="noether");
  con_unc <- concordance.index(x=unc_list, surv.time=as.numeric(as.character(unc_OS)), surv.event=as.numeric(as.character(unc_OS_Status)), na.rm=TRUE, method="noether");
  con_icgc_array <- concordance.index(x=icgc_array_list, surv.time=as.numeric(as.character(icgc_array_OS)), surv.event=as.numeric(as.character(icgc_array_OS_Status)), na.rm=TRUE, method="noether");
  con_collisson<- concordance.index(x=collisson_list, surv.time=as.numeric(as.character(collisson_OS)), surv.event=as.numeric(as.character(collisson_OS_Status)), na.rm=TRUE, method="noether");
  con_kirby<- concordance.index(x=kirby_list, surv.time=as.numeric(as.character(kirby_OS)), surv.event=as.numeric(as.character(kirby_OS_Status)), na.rm=TRUE, method="noether");
  con_pcsi<- concordance.index(x=pcsi_list, surv.time=as.numeric(as.character(pcsi_OS)), surv.event=as.numeric(as.character(pcsi_OS_Status)), na.rm=TRUE, method="noether");
  
  ####################################### Calculating meta-estimates of D-index and Concordance index
  ###  Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR OVERALL DATA
 dindex_meta <- combine.est(c( dindex_pcsi$d.index,  dindex_tcga$d.index, dindex_kirby$d.index, dindex_ouh$d.index,dindex_winter$d.index,dindex_unc$d.index,dindex_icgc_array$d.index, dindex_collisson$d.index),
                             c(dindex_pcsi$se, dindex_tcga$se, dindex_kirby$se, dindex_ouh$se,dindex_winter$se,dindex_unc$se,dindex_icgc_array$se,dindex_collisson$se),na.rm = TRUE,hetero = TRUE)
  dindex_meta_lower <- dindex_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_meta$se
  dindex_meta_upper <- dindex_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_meta$se
 # dindex_meta_pval <- combine.test(p=c(dindex_pcsi$p.value,  dindex_tcga$p.value, dindex_kirby$p.value, dindex_ouh$p.value,dindex_winter$p.value,dindex_unc$p.value,dindex_icgc_array$p.value, dindex_collisson$p.value),
  #                                 w=c(length(pcsi_list), length(tcga_list),length(kirby_list),length(ouh_list),length(winter_list),length(unc_list),length(icgc_array_list), length(collisson_list)),hetero = FALSE,method="z.transform")
  
  
  dindex_meta_pval <- 2*pnorm(-abs(log(dindex_meta$estimate)/dindex_meta$se))

  con_meta <- combine.est(c(con_pcsi$c.index, con_tcga$c.index, con_kirby$c.index, con_ouh$c.index,con_winter$c.index,con_unc$c.index,con_icgc_array$c.index, con_collisson$c.index),
                          c(con_pcsi$se,  con_tcga$se, con_kirby$se,con_ouh$se,con_winter$se,con_unc$se,con_icgc_array$se, con_collisson$se),na.rm = TRUE,hetero = TRUE)
  con_meta_lower <- con_meta$estimate + qnorm(0.025, lower.tail=TRUE) *  con_meta$se
  con_meta_upper <- con_meta$estimate + qnorm(0.025, lower.tail=FALSE) *  con_meta$se
  #con_meta_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value, con_kirby$p.value,con_ouh$p.value,con_winter$p.value,con_unc$p.value,con_icgc_array$p.value, con_collisson$p.value),
   #                             w=c(length(pcsi_list), length(tcga_list),length(kirby_list), length(ouh_list),length(winter_list),length(unc_list),length(icgc_array_list), length(collisson_list) ),method="z.transform")
    
    con_meta_pval <- pnorm((con_meta$estimate -0.5)/con_meta$se, lower.tail = con_meta$estimate < 0.5) * 2

  
  ### Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR sequencing cohort

  dindex_seq <- combine.est(c( dindex_pcsi$d.index,dindex_tcga$d.index, dindex_kirby$d.index),c(dindex_pcsi$se,  dindex_tcga$se, dindex_kirby$se),na.rm = TRUE, hetero = TRUE) 
  dindex_seq_lower <- dindex_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_seq$se
  dindex_seq_upper <- dindex_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_seq$se
  #dindex_seq_pval <- combine.test(p=c(dindex_pcsi$p.value,dindex_tcga$p.value, dindex_kirby$p.value),
   #                               w=c(length(pcsi_list),length(tcga_list), length(kirby_list)),method="z.transform")
  
  dindex_seq_pval<- 2*pnorm(-abs(log(dindex_seq$estimate)/dindex_seq$se))

  con_seq <- combine.est(c( con_pcsi$c.index, con_tcga$c.index, con_kirby$c.index),c(con_pcsi$se,  con_tcga$se, con_kirby$se),na.rm = TRUE,  hetero = TRUE)
  con_seq_lower <- con_seq$estimate + qnorm(0.025, lower.tail=TRUE) *  con_seq$se
  con_seq_upper <- con_seq$estimate + qnorm(0.025, lower.tail=FALSE) *  con_seq$se
 # con_seq_pval <- combine.test(p=c(con_pcsi$p.value, con_tcga$p.value, con_kirby$p.value),
  #                             w=c(length(pcsi_list), length(tcga_list), length(kirby_list)),method="z.transform")
  
    con_seq_pval<- pnorm((con_seq$estimate -0.5)/con_seq$se, lower.tail = con_seq$estimate < 0.5) * 2
  
  ### Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR microarray cohort
  
  dindex_micro <- combine.est(c( dindex_ouh$d.index,dindex_winter$d.index,dindex_unc$d.index,dindex_icgc_array$d.index, dindex_collisson$d.index),c( dindex_ouh$se,dindex_winter$se,dindex_unc$se,dindex_icgc_array$se, dindex_collisson$se),na.rm = TRUE, hetero = TRUE)
  dindex_micro_lower <- dindex_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  dindex_micro$se
  dindex_micro_upper <- dindex_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  dindex_micro$se
  #dindex_micro_pval <- combine.test(p=c(dindex_ouh$p.value,dindex_winter$p.value,dindex_unc$p.value,dindex_icgc_array$p.value, dindex_collisson$p.value),
   #                                 w=c(length(ouh_list),length(winter_list),length(unc_list),length(icgc_array_list), length(collisson_list)),hetero = FALSE,method="z.transform")
  
dindex_micro_pval<- 2*pnorm(-abs(log(dindex_micro$estimate)/dindex_micro$se))
  
  con_micro <- combine.est(c(  con_ouh$c.index,con_winter$c.index,con_unc$c.index,con_icgc_array$c.index, con_collisson$c.index),c( con_ouh$se,con_winter$se,con_unc$se,con_icgc_array$se, con_collisson$se),na.rm = TRUE, hetero = TRUE)
  con_micro_lower <- con_micro$estimate + qnorm(0.025, lower.tail=TRUE) *  con_micro$se
  con_micro_upper <- con_micro$estimate + qnorm(0.025, lower.tail=FALSE) *  con_micro$se
  #con_micro_pval <- combine.test(p=c(con_ouh$p.value,con_winter$p.value,con_unc$p.value,con_icgc_array$p.value, con_collisson$p.value),
   #                              w=c(length(ouh_list),length(winter_list),length(unc_list),length(icgc_array_list), length(collisson_list)),hetero = FALSE,method="z.transform")
  
      con_micro_pval<- pnorm((con_micro$estimate -0.5)/con_micro$se, lower.tail = con_micro$estimate < 0.5) * 2

  
  results = list(dindex_meta=dindex_meta, dindex_meta_lower=dindex_meta_lower, dindex_meta_upper=dindex_meta_upper, dindex_meta_pval=dindex_meta_pval,
                 dindex_seq=dindex_seq,dindex_seq_lower= dindex_seq_lower, dindex_seq_upper=dindex_seq_upper,dindex_seq_pval= dindex_seq_pval,
                 dindex_micro=dindex_micro,dindex_micro_lower= dindex_micro_lower, dindex_micro_upper=dindex_micro_upper,dindex_micro_pval= dindex_micro_pval,
                 con_meta=con_meta,con_meta_lower= con_meta_lower, con_meta_upper=con_meta_upper,con_meta_pval= con_meta_pval,
                 con_seq=con_seq,con_seq_lower= con_seq_lower, con_seq_upper=con_seq_upper,con_seq_pval= con_seq_pval,
                 con_micro=con_micro,con_micro_lower= con_micro_lower, con_micro_upper=con_micro_upper,con_micro_pval= con_micro_pval,
                 
                 con_pcsi=con_pcsi, con_tcga=con_tcga, con_ouh=con_ouh, con_collisson= con_collisson, con_icgc_array= con_icgc_array, con_kirby= con_kirby, con_unc=con_unc, con_winter=con_winter,
                 dindex_pcsi=dindex_pcsi, dindex_tcga=dindex_tcga, dindex_ouh=dindex_ouh, dindex_collisson= dindex_collisson, dindex_icgc_array= dindex_icgc_array, dindex_kirby= dindex_kirby, dindex_unc=dindex_unc, dindex_winter=dindex_winter)
  
  return(results) 
  
}

