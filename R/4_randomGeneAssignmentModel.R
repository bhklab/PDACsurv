#' Fit a random gene assignment model?
#'
##TODO:: HEEWON - documentation
#'
#' @examples
#' data(training_cohorts)
#' randomGeneAssignmentModel(training_cohorts)
#'
#' @param data A \code{} containing stuff
#' @param saveDir A \code{character} vector containing the path to the desired
#'   save directory
#'
#' @return A \code{}
#'
#' @importFrom switchBox SWAP.KTSP.Classify
#' @importFrom reportROC reportROC
#' @export
randomGeneAssignmentModel <- function(data, saveDir) {

  icgc_seq_cohort <- data$icgc_seq_cohort
  icgc_array_cohort <- data$icgc_array_cohort

  # Excluding samples censored before 1-yr ----------------------------------
  merge_common <- mergeCommonData(icgc_seq_cohort,
                                  icgc_array_cohort)

  # Training the model on ICGC seq/array common samples cohort --------------

  # Generating 1000 random models by re-shuffling the labels
  random_gene_model <- reshuffleRandomModels(merge_common)

  if (!missing(saveDir)) {
    saveRDS(random_gene_model,
            file=file.path(saveDir, 'randomGeneModel.rds'),
            compress='bzip2')
  }

  # Predicting probabliliets using random model for all the validation cohorts

  validationCohorts <- formatValidationCohorts(PDCAexpressionData)

  ktspPredictionList <-

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
    meta_auc[[i]] = combine.est(c(pcsi_list[[2]][i],
                                  tcga_list[[2]][i],
                                  unc_list[[2]][i],
                                  zhang_list[[2]][i],
                                  winter_list[[2]][i],
                                  ouh_list[[2]][i],
                                  icgc_arr_list[[2]][i], chen_list[[2]][i], kirby_list[[2]][i], collisson_list[[2]][i]), c( pcsi_list[[3]][i], tcga_list[[3]][i], unc_list[[3]][i], zhang_list[[3]][i], winter_list[[3]][i], ouh_list[[3]][i],icgc_arr_list[[3]][i], chen_list[[3]][i], kirby_list[[3]][i], collisson_list[[3]][i]),na.rm=TRUE)$estimate
    seq_auc[[i]] = combine.est(c(pcsi_list[[2]][i], tcga_list[[2]][i], kirby_list[[2]][i]), c( pcsi_list[[3]][i], tcga_list[[3]][i],kirby_list[[3]][i]),na.rm=TRUE)$estimate
    microarray_auc[[i]] = combine.est(c( unc_list[[2]][i], zhang_list[[2]][i],winter_list[[2]][i], ouh_list[[2]][i],icgc_arr_list[[2]][i], chen_list[[2]][i], collisson_list[[2]][i]), c( unc_list[[3]][i], zhang_list[[3]][i], winter_list[[3]][i], ouh_list[[3]][i],icgc_arr_list[[3]][i], chen_list[[3]][i], collisson_list[[3]][i]),na.rm=TRUE)$estimate

  }

}

.predict_ktsp <- function(list){
  val_pred <- list()
  val_pred_freq<- list()
  auc=vector()
  auc_se=vector()

  for(i in 1: length(random_gene_model) ){
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), random_gene_model[[i]])

    a <- reportROC(val_grp, as.numeric(as.character(val_pred[[i]])),plot = FALSE)
    auc[i]=a$AUC
    auc_se[i]=a$AUC.SE

  }

  ret_list=list(val_pred, auc, auc_se)
  return(ret_list)
}
