#' Build PCOSP Model from input data
#'
##TODO:: HEEWON: What does this do? Write a description (two or three seneteces max)
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
#' buldPCOSPmodel(training_cohort, saveDir=tempdir())
#' 
#' 
##TODO:: Determine where this dataset came from? Is it the ouput of another
#   package? From a publication?
#' @param \code{list} A dataset for which to build the PCOSP model
#' @param \code{character} A path to a directory to save the model. If you 
#'   exclude this the function will return the model object instead.
#' 
#' @return \code{?} Either returns the model object or, is \code{saveDir} is 
#'   specified it saves to disk instead and return the path
#'
#' @export
#' 
##TODO:: Convert these to import from!
#' @import switchBox switchBox vcdExtra caret forestplot ktspair pROC survcomp 
#' @import survival data.table reportROC verification
#'
buildPCOSPmodel <- function(data, saveDir) {
  
  icgc_seq_cohort = training_cohort$icgc_seq_cohort
  icgc_array_cohort = training_cohort$icgc_array_cohort
  
  rownames(icgc_array_cohort) == rownames(icgc_seq_cohort)
  
  ## Excluding samples censored before 1-yr
  g1 <- which(as.numeric(as.character(icgc_seq_cohort$OS))<=365 & 
             as.numeric(as.character(icgc_seq_cohort$OS_Status))==1)
  g2 <- which(as.numeric(as.character(icgc_seq_cohort$OS))>365)
  g_ind <- sort(c(g1,g2))
  
  icgc_seq_cohort <- icgc_seq_cohort[g_ind,]
  icgc_array_cohort <- icgc_array_cohort[g_ind,]
  
  merge_common <- rbind(icgc_seq_cohort,icgc_array_cohort) # Merged common ICGC seq and array data
  
  
  ## Classes for training
  xx = merge_common
  xmat <- xx[seq_len(nrow(xx)), seq_len(ncol(xx)-2)] # Removing survival data columns
  merge_common_mat <- data.matrix(sapply(xmat, function(xx) as.numeric(as.character(xx))))
  rownames(merge_common_mat) <- rownames(merge_common)
  merge_common_grp <- ifelse(as.numeric(as.character(xx$OS)) >= 365, 1, 0)
  
  ## Generating 1000 TSP models
  pred <- list()
  sel_pred <- list()
  count <- 0;
  b_acc <- vector()
  F1 <- vector()
  i=1
  model <- list()
  models_no = 1000
  count <- 1
  selected_model <- list()
  set.seed(1987) ##FIXME:: CRAN will complain about this; probably Bioconductor as well
                  # the recommendation is to have the user set seed OUTSIDE of the function
  sel_b_acc=list()
  
  for(i in seq_len(1000)){
    
    x5 <-sample(which(merge_common_grp==0), 40, replace=F) # Selecting random 40 samples from Low survival group
    y5 <-sample(which(merge_common_grp==1), 40, replace=F) # Selecting random 40 samples from High survival group
    
    x1 <- merge_common_mat[c(x5,y5),]                
    y_index <- c(x5, y5) # Selecting the classes of re-sampled samples
    y1 <- merge_common_grp[y_index]
    
    
    ## Building k-TSP models
    zzz <- paste0('classifier', i) 
    model[[i]] <- SWAP.KTSP.Train(t(x1), as.factor(y1))
    
    z <- setdiff(1:164,c(x5,y5)) # Finding test samples excluded in training set
    test <- merge_common_mat[z,]   
    test_grp <-merge_common_grp[z]
    
    ### Testing the model on out of bag samples
    pred[[i]] <- SWAP.KTSP.Classify(t(test), model[[i]])   ### Predicting the classes of test set
    
    cc=confusionMatrix(pred[[i]], test_grp,  mode = "prec_recall")
    b_acc[i]=as.numeric(cc$byClass)[11]
    F1[i]=as.numeric(cc$byClass)[7]
    print(i)
  }
  
  selected_model = model[which(b_acc> 0.60)]
    
  ##FIXME:: Do we need this? If so I can make it a message; if not delete it.
  length(selected_model)
  
  if (!missing(saveDir)) {
    save(selected_model, file=saveDir)
    return(paste0('Saved model to ', saveDir))
  } else {
    return(selected_model)
  }
}

