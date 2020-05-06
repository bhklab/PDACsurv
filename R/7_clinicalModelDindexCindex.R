#' Calculate the
#'
#' @param modelProbabilities A
#' @param clinicalFeatures A
#' @param seqCohorts A
#' @param model A
#'
#' @export
calculateModelDandCindex <- function(modelProbabilities, clinicalFeatures,
                                     seqCohorts, model=1) {

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

  ##TODO:: Determine why we invert the probabiltiies here and if the names
  ##    are still correct?
  clinicalProbs <- lapply(modelProbabilities[[model]], function(cohort) 1 - cohort$clinical)
  PCOSPprobs <-lapply(modelProbabilities[[model]], function(cohort) 1 - cohort$PCOSP)

  clinicalStats <- constructMetaEstimatesDF(clinicalProbs, cFeatures, seqCohorts)
  PCOSPstats <- constructMetaEstimatesDF(PCOSPprobs, cFeatures, seqCohorts)

  return(list(
    "clinical"=clinicalStats,
    "PCOSP"=PCOSPstats
  ))
}

###########################

# pcosp_clinical_cindex = cindex.comp.meta(list.cindex1 = list(con_pcsi1, con_tcga1, con_ouh1, con_icgc_array1),
#                                   list.cindex2 = list(con_pcsi, con_tcga, con_ouh, con_icgc_array))
# pcosp_clinical_cindex
#
# pcosp_clinical_dindex = dindex.comp.meta(list.dindex1 = list(dindex_pcsi1, dindex_tcga1, dindex_ouh1, dindex_icgc_array1),
#                                   list.dindex2 = list(dindex_pcsi, dindex_tcga, dindex_ouh, dindex_icgc_array))
# pcosp_clinical_dindex
