#' Validate the Meta Estimate Calculation
#'
##TODO:: HEEWON - what does this function do?
#'
#' The function calculates D-index and Concordance index for independent cohorts
#'   and calculates the meta-estimae of C-index and D-index. The forestplot of
#'   C-index and D-index was plotted for all the cohorts.
#'
#' @example
#' # Load validation data
#' data(validationCohorts)
#' data(selectedModels)
#'
#TODO:: Set nthread to 1
#' # Validate the meta-estimate calculation
#' validationStats <- validateMetaEstimateCalculation(validationCohorts, selectedModels,
#'     seqCohorts=c("PCSI", "TCGA", "Kirby"), nthread=15)
#'
#' @param validationCohorts A named \code{list} of validation cohorts.
#' @param selectedModels A \code{list} of selected models from the
#'     `buildPCOSPmodels` function.
#' @param seqCohorts A \code{character} vector specifying which cohorts
#'     are from sequencing data. It is assumed that non-sequencing cohorts
#'     are from micro-array cohorts.
#' @param nthread A \code{numeric} vector with an integer number of thread to
#'     parallelize across.
#' @param saveDir A \code{character} vector specifying the directory in which
#'     to save results. If missing the function returns the data.
#'
#' @return A named \code{list} of statistical meta-estimates for the combined,
#'   sequencing and microarray cohorts supplied in validationCohorts. Top level
#'   list labels indicate which cohorts, with the subsequent list level
#'   indicating the type of statistic (i.e., dindex vs concordance index)
#'
#' @importFrom survcomp D.index concordance.index combine.est
#' @importFrom forestplot forestplot
#' @importFrom reportROC reportROC
#' @importFrom pROC roc
#' @importFrom verification roc.area
#' @export
calculateValidationStats <- function(validationCohorts, selectedModels, seqCohorts,
                                nthread, saveDir) {

  PCOSPscoreList <- calculatePCOSPscores(validationCohorts, selectedModels, nthread)

  ## Dindex estimate calculation
  DindexList <- .estimateDindex(PCOSPscoreList, validationCohorts)

  ## Concordance index calculation
  concordanceIndexList <- .estimateConcordanceIndex(PCOSPscoreList,
                                                    validationCohorts)

  ## Calculating meta-estimates of D-index and Concordance index

  ##  Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR OVERALL DATA
  combinedStats <- metaEstimateStats(DindexList, concordanceIndexList,
                                     hetero=TRUE)

  ## Determine which cohorts are from sequencing data
  isSeq <- names(validationCohorts) %in% seqCohorts

  ## Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR sequencing cohort
  sequencingStats <- metaEstimateStats(DindexList[isSeq],
                                       concordanceIndexList[isSeq],
                                       hetero=FALSE)

  ## Meta-estimate of d-INDEX AND CONCORDANCE INDEX FOR microarray cohort
  arrayStats <- metaEstimateStats(DindexList[!isSeq],
                                  concordanceIndexList[!isSeq],
                                  hetero=FALSE)

  # Extract statistics for Dindex and concordanceIndex into a list of data.frames
  list(
    "dIndex"=data.frame(
      "mean"=c(vapply(DindexList, function(x) log2(x$d.index), FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) log2(x$metaEstimate$dIndex$estimate), FUN.VALUE=numeric(1))),
      "lower"=c(vapply(DindexList, function(x) log2(x$lower), FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) log2(x$lowerTail$dIndex), FUN.VALUE=numeric(1))),
      "upper"=c(vapply(DindexList, function(x) log2(x$upper), FUN.VALUE=numeric(1)),
                 vapply(list(sequencingStats, arrayStats, combinedStats),
                        function(x) log2(x$upperTail$dIndex), FUN.VALUE=numeric(1))),
      "pval"=c(vapply(DindexList, function(x) x$p.value, FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) x$pValue$dIndex, FUN.VALUE=numeric(1))),
      row.names=c(names(DindexList), "Sequencing", "Microarray", "Overall")
    ),
    "cIndex"=data.frame(
      "mean"=c(vapply(concordanceIndexList, function(x) x$c.index, FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) x$metaEstimate$cIndex$estimate, FUN.VALUE=numeric(1))),
      "lower"=c(vapply(concordanceIndexList, function(x) x$lower, FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) x$lowerTail$cIndex, FUN.VALUE=numeric(1))),
      "upper"=c(vapply(concordanceIndexList, function(x) x$upper, FUN.VALUE=numeric(1)),
                vapply(list(sequencingStats, arrayStats, combinedStats),
                       function(x) x$upperTail$cIndex, FUN.VALUE=numeric(1))),
      "pval"=c(vapply(concordanceIndexList, function(x) x$p.value, FUN.VALUE=numeric(1)),
               vapply(list(sequencingStats, arrayStats, combinedStats),
                      function(x) x$pValue$cIndex, FUN.VALUE=numeric(1))),
      row.names=c(names(concordanceIndexList), "Sequencing", "Microarray", "Overall")
    ),
    "PCOSPscores"=PCOSPscoreList,
    "isSequencing"=isSeq
  )
}

#' Calculate the cohort-wise PCOSP score from a list of validation cohorts
#'
#' @param validationCohorts A \code{list} of validation cohorts
#' @param selectedModels A \code{list} of selected models, as returned by the
#'   `buildPCOSPmodels` function.
#' @param nthread
#'
#' @return A \code{list} of PCOSP scores for each validation cohorot
#'
#' @export
calculatePCOSPscores <- function(validationCohorts, selectedModels, nthread) {

  formattedValCohorts <- formatValidationCohorts(validationCohorts)

  # Temporily change number of cores to parallelize over
  opts <- options()
  options("mc.cores"=nthread)
  on.exit(options(opts))

  ## Extract matrices from the cohort list
  cohortMatrixList <- lapply(formattedValCohorts, function(cohort) cohort$mat)
  ##TODO:: This is not used here, do we need it?
  cohortGroupList <- lapply(formattedValCohorts, function(cohort) cohort$grp)

  # Estimate PCOSP scores from cohort matrixes
  PCOSPscoreList <- bplapply(cohortMatrixList,
                             function(cohortMat, selectedModels)
                               estimatePCOSPprob(cohortMat, selectedModels),
                             selectedModels=selectedModels)
  return(PCOSPscoreList)
}

#'
#'
#'
#'
metaEstimateStats <- function(DindexList, concordanceIndexList, hetero) {
  DindexMetaEstimate <- .metaEstimateDindex(DindexList, hetero)
  concordanceIndexMetaEstimate <-
    .metaEstimateConcordanceIndex(concordanceIndexList, hetero)
  stats <- .zipLists(DindexMetaEstimate,
                     concordanceIndexMetaEstimate,
                     sublistNames=c("dIndex", "cIndex"))
  return(stats)
}

#' Estimate the D-index for a list of validation cohorts
#'
#'
#'
#' @importFrom survcomp D.index
.estimateDindex <- function(PCOSPscoreList, validationCohorts) {
  structure(lapply(seq_along(PCOSPscoreList), function(i) {
    cohort <- validationCohorts[[i]]
    PCOSPscore <- PCOSPscoreList[[i]]
    D.index(x=PCOSPscore$predicted_probabilities,
            surv.time=as.numeric.factor(cohort$OS),
            surv.event=as.numeric.factor(cohort$OS_Status),
            na.rm=TRUE,
            alpha=0.05,
            method.test="logrank")
  }), .Names=names(PCOSPscoreList))
}

#' Estimate the concordance index for a list of validation cohorts
#'
#' @param PSCOSPscoreList A named \code{list} of PCOSPscores as calculated
#'   with `estimatePCOSPscore`
#' @param validationCohorts A named \code{list} of validation cohorts
#'
#' @return
#'
#' @importFrom survcomp concordance.index
.estimateConcordanceIndex <- function(PCOSPscoreList, validationCohorts) {
  structure(lapply(seq_along(validationCohorts), function(i) {
    cohort <- validationCohorts[[i]]
    PCOSPscore <- PCOSPscoreList[[i]]
    concordance.index(x=PCOSPscore$predicted_probabilities,
            surv.time=as.numeric.factor(cohort$OS),
            surv.event=as.numeric.factor(cohort$OS_Status),
            method="noether")
  }), .Names=names(validationCohorts))
}

#' Meta-estimate the Dindex for a list of validation cohorts
#'
#'
#'
#' @importFrom survcomp combine.est
.metaEstimateDindex <- function(DindexList, hetero) {

  Dindexes <- vapply(DindexList, function(cohort) cohort$d.index,
                     FUN.VALUE=numeric(1))
  DindexSEs <- vapply(DindexList, function(cohort) cohort$se,
                      FUN.VALUE=numeric(1))

  DindexMetaEstimate <- combine.est(Dindexes, DindexSEs, na.rm=TRUE, hetero=hetero)

  ##TODO:: Define a metaestimate S4 object? Can then dispatch plots on it
  return(list(
    "metaEstimate"= DindexMetaEstimate,
    "lowerTail"=DindexMetaEstimate$estimate + qnorm(0.025, lower.tail=TRUE) *
      DindexMetaEstimate$se,
    "upperTail"=DindexMetaEstimate$estimate + qnorm(0.025, lower.tail=FALSE) *
      DindexMetaEstimate$se,
    "se"=DindexMetaEstimate$se,
    "pValue"=2*pnorm(-abs(log(DindexMetaEstimate$estimate)/DindexMetaEstimate$se)),
    "cohortNames"=names(DindexList)
  ))
}

#' Meta-estimate the concordance index for
#'
#' @importFrom survcomp concordance.index
.metaEstimateConcordanceIndex <- function(concordanceIndexList, hetero) {

  conIndexes <- vapply(concordanceIndexList, function(cohort) cohort$c.index,
                     FUN.VALUE=numeric(1))
  conIndexSEs <- vapply(concordanceIndexList, function(cohort) cohort$se,
                      FUN.VALUE=numeric(1))

  conIndexMetaEstimate <- combine.est(conIndexes, conIndexSEs, na.rm=TRUE, hetero=hetero)

  ##TODO:: Define a metaestimate S4 object? Can then dispatch plots on it
  return(list(
    "metaEstimate"= conIndexMetaEstimate,
    "lowerTail"=conIndexMetaEstimate$estimate + qnorm(0.025, lower.tail=TRUE) *
      conIndexMetaEstimate$se,
    "upperTail"=conIndexMetaEstimate$estimate + qnorm(0.025, lower.tail=FALSE) *
      conIndexMetaEstimate$se,
    "se"=conIndexMetaEstimate$se,
    "pValue"=2*pnorm((conIndexMetaEstimate$estimate - 0.5)/conIndexMetaEstimate$se,
                     lower.tail= conIndexMetaEstimate$estimate < 0.5),
    "cohortNames"=names(concordanceIndexList)
  ))
}

##TODO:: Should this go in utilties?
#' Merge n lists element-wise into a list of sublists n items long
#'
#' Take n lists and combine them element-wise into a single list with each sublist
#'   containing the corresponding elements of the original lists. Assigns
#'   sublistNames as the names of the new zipped sublists.
#'
#' @param ... Two or more \code{list}s to zip together element-wise. New list
#'     is named using first list.
#' @param sublistNames A \code{character} vector n items long, specifying the
#'     labels for the zipped items
#'
#' @keywords internal
.zipLists <- function(..., sublistNames) {
  # Error checking
  if (length(list(...)) != length(sublistNames))
    stop("Must be as many names as lists!")

  # Merge lists into sublists element-wise
  zipped <- mapply(list, ..., SIMPLIFY=FALSE)

  # Name each zipped item in the sublist
  zipped <- lapply(zipped, function(stat) structure(stat, .Names=sublistNames))
  return(zipped)
}














































































































.plottingFunctions <- function(refactorMe) {
  ## Plotting Forestplot of  D index
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





  # Plotting Forestplot of Concordance index --------------------------------

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
  # Plotting ROC CURVES -----------------------------------------------------

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

  ouh_roc=roc(ouh_grp, ouh_list[[1]][g_ouh])$auc[1]
  ouh_roc_se=reportROC(ouh_grp, ouh_list[[1]][g_ouh], plot = FALSE)$AUC.SE

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
}
