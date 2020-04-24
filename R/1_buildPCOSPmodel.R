#' Build PCOSP Model from input data
#'
##TODO:: HEEWON: What does this do? Write a description (two or three senetences max)
#'
#' The script is used for building Pancreatic cancer overall survival predictor
#'   using unique 89 samples profiled using microarray and sequencing platform.
#'
#' @examples
#' # To ensure reproducible results
#' set.seed(1234)
#'
#' # Load the data
#' data(training_cohort)
#'
#' # Return the object
#' model <- buildPCOSPmodel(training_cohort)
#'
#' # Save the object to disk
#' buildPCOSPmodel(training_cohort, saveDir=tempdir())
#'
##TODO:: Determine where this dataset came from? Is it the ouput of another
#   package? From a publication?
#' @param data \code{list} A dataset for which to build the PCOSP model
#' @param saveDir \code{character} A path to a directory to save the model. If you
#'   exclude this the function will return the model object instead.
#'
#' @return \code{?} Either returns the model object or, is \code{saveDir} is
#'   specified it saves to disk instead and return the path
#'
#' @warning This function uses random numbers; remeber to \code{set.seed()}
#'   before running to ensure reproducible results
#'
#' @export
buildPCOSPmodel <- function(data, saveDir) {

  icgc_seq_cohort <- data$icgc_seq_cohort
  icgc_array_cohort <- data$icgc_array_cohort

  # Merged common ICGC seq and array data
  merge_common <- mergeCommonData(icgc_seq_cohort, icgc_array_cohort)

  # Training the model on ICGC seq/array common samples cohort
  merge_common_mat <- convertCohortToMatrix(merge_common)
  merge_common_grp <- ifelse(as(merge_common$OS, 'numeric') >= 365, 1, 0)

  selected_model <- .generateTSPmodels(merge_common_mat, merge_common_grp)

  # Save to disk or return
  if (!missing(saveDir)) {
    save(selected_model, file=saveDir)
    return(paste0('Saved model to ', saveDir))
  } else {
    return(selected_model)
  }
}


##TODO:: See if we can refactor part of this to be reused in reshuffleRandomModels
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Train
.generateTSPmodels <- function(merge_common_mat, merge_common_grp) {
  ## Generating 1000 TSP models
  ##TODO:: Set some of these as function parameters?
  pred <- list()
  count <- 0;
  b_acc <- vector()
  F1 <- vector()
  i=1
  model <- list()
  models_no <- 1000
  count <- 1
  selected_model <- list()

  for(i in seq_len(models_no)){

    # Selecting random 40 samples from Low survival group
    x5 <-sample(which(merge_common_grp == 0), 40, replace=FALSE)
    # Selecting random 40 samples from High survival group
    y5 <-sample(which(merge_common_grp == 1), 40, replace=FALSE)
    x1 <- merge_common_mat[c(x5,y5),]
    y_index <- c(x5, y5) # Selecting the classes of re-sampled samples
    y1 <- merge_common_grp[y_index]


    ## Building k-TSP models
    model[[i]] <- SWAP.KTSP.Train(t(x1), as.factor(y1))

    z <- setdiff(1:164,c(x5,y5)) # Finding test samples excluded in training set
    test <- merge_common_mat[z,]
    test_grp <-merge_common_grp[z]

    ### Testing the model on out of bag samples
    ### Predicting the classes of test set
    pred[[i]] <- SWAP.KTSP.Classify(t(test), model[[i]])

    cc=confusionMatrix(pred[[i]], test_grp,  mode="prec_recall")
    b_acc[i]=as.numeric(cc$byClass)[11]
    F1[i]=as.numeric(cc$byClass)[7]
    print(i)
  }

  selected_model <- model[which(b_acc> 0.60)]
  return(selected_model)
}

