#' Calculate the probability
#'
##TODO:: HEEWON add function documentation
#'
#' @param A A \code{matrix} of values
#' @param selectedModels a \code{list} of models selected using the
#'   `buildPCOSPmodel` function
#'
#' @return A \code{list} containing the probabilities
#'
#' @importFrom switchBox SWAP.KTSP.Classify
estimatePCOSPprob <- function(val_mat, selectedModels){
  val_pred <- list()

  for(i in seq_along(selectedModels)) {
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), selectedModels[[i]])
  }

  list_z1=vector()
  freq_early=vector()
  i=1

  for (i in seq_len(nrow(val_mat))) {
    for (k in seq_along(selectedModels)) {
      list_z1 <- append(list_z1, as.numeric(val_pred[[k]][i]))
    }

    freq_early[i] <- length(list_z1[list_z1==1])/length(selectedModels)
    list_z1 <- vector()
  }
  ret_list=list(predicted_probabilities=freq_early)
  return(ret_list)
}