#' Build PCOSP models with the group labels randomly shuffled
#'
## TODO:: HEEWON Document this function
#'
#' @param trainingCohorts A named \code{list} of training cohorts for which
#'     to fit and select PCOSP models.
#' @param numModels A code{numeric} vector with an integer number of models
#'     to run.
#' @param saveDir \code{character} A path to a directory to save the model. If you
#'   exclude this the function will return the model object instead.
#' @param nthread \code{integer} The number of threads to parallelize across
#' @param original A \code{logical} vector, if true this function calls
#'    the deprecated (and slow) function used in the original PCOSP paper.
#'    This is inlcuded to ensure the paper is reproducible using this pacakge.
#'
#' @return \code{list} Either returns the model object or, is \code{saveDir} is
#'   specified it saves to disk instead and return the path
#'
#' @section Warning: This function uses random numbers; remember to
#'   \code{set.seed()} before running to ensure reproducible results
#'
#' @export
buildRandomGeneAssignmentModels <- function(trainingCohorts, numModels, nthread,
                                           saveDir, original) {

  # Set number of threads to parallelize over
  if (!missing(nthread)) {
    ops <- options()
    options("mc.cores"=nthread)
    on.exit(options(ops))
  }

  # Extract cohorts from trainingCohorts
  seqCohort <- trainingCohorts$icgc_seq_cohort
  arrayCohort <- trainingCohorts$icgc_array_cohort

  # Merged common ICGC seq and array trainingCohorts
  commonData <- mergeCommonData(seqCohort, arrayCohort)

  # Training the model on ICGC seq/array common samples cohort
  cohortMatrix <- convertCohortToMatrix(commonData)
  cohortMatrixGroups <- ifelse(as.numeric.factor(commonData$OS) >= 365, 1, 0)

  if (missing(original) || original == FALSE ) {
    selectedModels <- .generateRGAmodels(cohortMatrix, cohortMatrixGroups,
                                         numModels)
  } else {
    selectedModels <- .buildRGAmodels(cohortMatrix, cohortMatrixGroups, numModels)
  }

  # Save to disk or return
  if (!missing(saveDir)) {
    saveRDS(selectedModels,
            file=paste0(file.path(saveDir, 'RGAmodels'), '.rds'))
    return(paste0('Saved model to ', saveDir))
  } else {
    return(selectedModels)
  }
}

##TODO:: See if we can refactor part of this to be reused in reshuffleRandomModels
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Train
.generateRGAmodels <- function(cohortMatrix, cohortMatrixGroups, numModels) {

  trainingDataRowIdxs <- lapply(rep(40, numModels),
                                .randomSampleIndexShuffle,
                                labels=cohortMatrixGroups,
                                groups=sort(unique(cohortMatrixGroups)))
  system.time({
    trainedModels <- bplapply(trainingDataRowIdxs,
                              function(idx, data)
                                SWAP.KTSP.Train(t(data[idx, ]), levels(idx)),
                              data=cohortMatrix)
  })

  geneNames <- colnames(cohortMatrix)
  numGenes <- vapply(trainedModels, function(m) nrow(m$TSPs), FUN.VALUE=numeric(1))

  RGAmodels <- lapply(seq_along(trainedModels),
                               function(idx, models, n, genes) {
                                  models[[idx]]$TSPs[, 1] <- sample(geneNames, n[idx])
                                  models[[idx]]$TSPs[, 2] <- sample(geneNames, n[idx])
                                  return(models[[idx]])
                               },
                               models=trainedModels,
                               n=numGenes,
                               genes=geneNames)

  ##TODO:: Do we need this? Not selecting models based on balanced accuray?
  # testingDataRowIdxs <- lapply(trainingDataRowIdxs,
  #                              function(idx, rowIdx, labels)
  #                                structure(setdiff(rowIdx, idx),
  #                                          .Label=as.factor(
  #                                            labels[setdiff(rowIdx, idx)])),
  #                              rowIdx=seq_len(nrow(cohortMatrix)),
  #                              labels=cohortMatrixGroups)
  #
  #
  # predictions <- bplapply(seq_along(testingDataRowIdxs),
  #                         function(i, testIdxs, data, models)
  #                           SWAP.KTSP.Classify(t(data[testIdxs[[i]], ]),
  #                                              models[[i]]),
  #                         testIdxs=testingDataRowIdxs,
  #                         data=cohortMatrix,
  #                         models=trainedModels
  # )
  #
  #
  # confusionMatrices <- bplapply(seq_along(predictions),
  #                               function(i, predictions, labels)
  #                                 confusionMatrix(predictions[[i]],
  #                                                 levels(labels[[i]]),
  #                                                 mode="prec_recall"),
  #                               predictions=predictions,
  #                               labels=testingDataRowIdxs
  # )
  #
  # modelStats <- bplapply(confusionMatrices,
  #                        function(confMat) confMat$byClass)
  #
  # balancedAcc <- unlist(bplapply(modelStats,
  #                                function(model) model[c('Balanced Accuracy')]))
  #
  # selectedModels <- trainedModels[which(balancedAcc > 0.60)]

  return(RGAmodels)
}

#' Uses non parallelized computation of the random gene assignment model to
#'     replicate the results from the original PCOSP paper. This method is much
#'     slower when multiple cores are available and should only be used to
#'     reproduce the results of the PCOSP paper.
#'
#' @param
#'
#'
.buildRGAmodels <- function(cohortMatrix, cohortMatrixGroups, numModels) {
  randomGeneModels <- lapply(rep(40, numModels),
                             .fitRGAModel,
                             data=cohortMatrix,
                             labels=cohortMatrixGroups
  )
  return(randomGeneModels)
}

.fitRGAModel <- function(n, data, labels) {
  idx <- unlist(mapply(function(grp, labs) sample(which(labs == grp), n, replace=FALSE),
                       grp=sort(unique(labels)),
                       MoreArgs=list(labs=labels),
                       SIMPLIFY=FALSE))
  data <- data[idx, ]
  labels <- labels[idx]
  model <- SWAP.KTSP.Train(t(data), as.factor(labels))
  model$TSPs <- cbind(sample(colnames(data), nrow(model$TSPs)),
                      sample(colnames(data), nrow(model$TSPs)))
  return(model)
}
