#' Calculate the
#'
#' @param modelProbabilities
#' @param clinicalFeatures
#' @param seqCohorts
#' @param
#'
caclulateModelDandCindex <- function(modelProbabilities, clinicalFeatures,
                                     seqCohorts, model=1)
  {

  namesClinical <- lapply(modelProbabilities[[model]],
                          function(cohort)
                            names(cohort$clinical))

  # Subset clinical cohorts to only sample with clinical predictions
  cFeatures <- structure(lapply(seq_along(namesClinical),
                      function(i, namesClinical, clinicalFeatures)
                        clinicalFeatures[[i]][clinicalFeatures[[i]]$ID %in% namesClinical[[i]],],
                      namesClinical=namesClinical,
                      clinicalFeatures=clinicalFeatures),
                      .Names=names(clinicalFeatures))

  clinicalProbs <- lapply(modelProbabilities[[model]], function(cohort) cohort$clinical)
  PCOSPprobs <-lapply(modelProbabilities[[model]], function(cohort) cohort$PCOSP)

  clinicalStats <- constructMetaEstimatesDF(clinicalProbs, cFeatures, seqCohorts)
  PCOSPStats <- contructMetaEstimatesDF(PCOSPprobs, cFeatures, seqCohorts)

  return(list(
    "clinical"=clinicalStats,
    "PCOSP"=PCOSPstats
  ))
}


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
