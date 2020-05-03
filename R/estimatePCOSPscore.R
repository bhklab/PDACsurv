#' Calculate the probability
#'
##TODO:: HEEWON add function documentation
#'
#' @param valMat A \code{matrix} of expression values for a validation
#'     validation cohort
#' @param selectedModels a \code{list} of models selected using the
#'     `buildPCOSPmodel` function
#'
#' @return A \code{vector} containing the probabilities per patient
#'
#' @importFrom switchBox SWAP.KTSP.Classify
estimatePCOSPprob <- function(valMat, selectedModels) {
  predictions <- bplapply(selectedModels,
                        function(model, valMat)
                          as.numeric.factor(SWAP.KTSP.Classify(t(valMat), model)),
                        valMat=valMat)

  allPredictions <- do.call(rbind, predictions)
  colnames(allPredictions) <- rownames(valMat)

  predProbabilities <- (length(selectedModels) - colSums(allPredictions)) / length(selectedModels)
  return(predProbabilities)
}