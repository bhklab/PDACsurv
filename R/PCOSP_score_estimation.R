#' Calculate the probability
#' 
##TODO:: HEEWON add function documentation
#' 
#' @param val_mat A \code{matrix} of values
#'
#' @return A \code{list}
#' 
#' @importFrom switchBox SWAP.KTSP.Classify
#' @export
pcosp_prob <- function(val_mat){
  val_pred <- list()
  
  for(i in seq_along(selected_model)) {
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), selected_model[[i]])  
  }
  
  list_z1=vector()
  freq_early=vector()
  i=1
  
  for (i in seq_len(nrow(val_mat))) {
    for (k in seq_along(selected_model)) {
      list_z1 <- append(list_z1, as.numeric(val_pred[[k]][i]))
    }
    
    freq_early[i] <- length(list_z1[list_z1==1])/length(selected_model) 
    list_z1 <- vector()
  }
  ret_list=list(predicted_probabilities=freq_early)
  return(ret_list)
}