#' Exclude samples censored before year 1
#'
#' @param seqCohort A \code{}
#'
#' @return A sorted vector of indices
whichNotCensoredBeforeYearOne <- function(seqCohort) {
  idxNotCensored <- which(as(seqCohort$OS, 'numeric') <= 365 & 
                            seqCohort$OS_status == 1)
  idxNotYearOne <- which(as(seqCohort$OS, 'numeric') > 365)
  return(sort(c(idxNotCensored, idxNotYearOne)))
}

#' Merge common samples not censored before year 1
#' 
#' @param seqCohort A \code{}
#' @param arrayCohort A \code{}
#' 
#' @return A \code{data.frame} containing the merged cohorts
#' 
#' @export
mergeCommonData <- function(seqCohort, arrayCohort) {
  if (rownames(icgc_array_cohort) != rownames(icgc_seq_cohort)) {
    stop("Rownames do not match between the seqCohort and arrayCohort!")
  }
  
  notCensoredYearOne <- whichNotCensoredBeforeYearOne(seqCohort)
  return(rbind(seqCohort[notCensoredYearOne, ],
               arrayCohort[notCensoredYearOne, ])
         )
}

#' Reshuffle to create 1000 random gene models
#' 
#' Take in a \code{data.frame} of merged cohorts and generate 1000 random
#'  gene models from it
#' 
#' @warning This function uses random numbers. Please use \code{set.seed()} 
#'   before running to ensure reporoducibility of the output.
#' 
#' @param merge_common A \code{data.frame} containing the merged cohorts
#' 
#' @return A \code{list} containing 1000 random gene models
#'
#' @export
#' @importFrom caret confusionMatrix
#' @importFrom switchBox SWAP.KTSP.Classify
#FIXME:: See if we can refactor this further
#FIXME:: Can we remove commented out lines?
reshuffleRandomModels <- function(merge_common) {
  
  ## Classes for training
  merge_common_mat <- convertCohortToMatrix(merge_common)
  merge_common_grp <- ifelse(as(merge_common$OS, 'numeric')) >= 365, 1, 0)
  
  ##TODO:: Make these parameters in the function call
  pred <- list()
  #sel_pred <- list()
  count <- 0;
  b_acc <- vector()
  i <- 1
  model <- list()
  models_no <- 1000
  count <- 1
  selected_model <- list()
  ##FIXME:: Can't set seed in functions for Bioconductor, user has to do it
  # set.seed(1987)
  #sel_b_acc <- list()
  
  for (i in seq_len(models_no)) {
    #set.seed(1)
    # Selecting random 40 samples from short survival group
    x5 <-sample(which(merge_common_grp == 0), 40, replace=FALSE)
    # Selecting random 40 samples from long survival group
    y5 <-sample(which(merge_common_grp == 1), 40, replace=FALSE)
    x1 <- merge_common_mat[c(x5,y5),]                
    
    # Selecting the classes of re-sampled samples
    y_index=c(x5, y5)                  
    y1 <- merge_common_grp[y_index]
    
    ### k-TSP models
    model[[i]] <- SWAP.KTSP.Train(t(x1), as.factor(y1))
    
    selected_model[[count]]=model[[i]]
    selected_model[[count]]$TSPs[, 1] <- 
      sample(colnames(merge_common_mat), 
             length(selected_model[[count]]$TSPs[,1]))
    selected_model[[count]]$TSPs[, 2] <- 
      sample(colnames(merge_common_mat), 
             length(selected_model[[count]]$TSPs[,1]))
    
    count <- count + 1
    print(i)
    
    z <- setdiff(seq_len(164), c(x5, y5)) ### Out of bag samples
    test <- merge_common_mat[z,]   
    test_grp <- merge_common_grp[z]
    
    
    ### Predicting on the test samples
    # Testing on out of bag samples
    pred[[i]] <- SWAP.KTSP.Classify(t(test), selected_model[[i]]) 
    
    cc <- confusionMatrix(pred[[i]], test_grp,  mode="prec_recall")
    b_acc[i] <- as.numeric(cc$byClass)[11]
    
  }
  random_gene_model <- selected_model
  return(random_gene_model)
}

#' Extract all cohorts from a list of cohorts
#'
#' Input a named list of cohorts to be modelled and assign each cohort to an
#'   environment (the global environment is recommended) 
#'
#' @param cohortList A named \code{list} containing the cohorts for an analysis 
#' @param environment The \code{environment} into which the variables will be 
#'   assigned. Use \code{globalenv()} to assign them to the global environment
#'   (recommended)
#'
#' @return Side effects only: assigns a variable to the parent environment with
#'   the name of each list element appended with "_cohort"
#'
#' @export
##TODO:: Should I use list2env function instead of assign
extractAllCohorts <- function(cohortList, environment) {
  if (missing(environment)) stop("Please set a target environment to assign
                                the variables into. We recommend 
                                 environment=globalenv()")
  for (i in seq_along(cohortList)) {
    assign(tolower(paste0(names(cohortList)[i], '_cohort')), 
           cohortList[[i]],
           envir=environment 
           )
  }
}

#'
#'
#'
#'
#'
#'
removeNotCensoredBeforeYearOne <- function(cohort) {
  g <- whichNotCensoredBeforeYearOne(cohort)
}


#' Take in cohort data.frame and return the cohort as a numeric matrix
#'
#' @param cohort A \code{data.frame} containing cohort data
#' 
#' @return A numeric \code{matrix} containing the data from the cohort 
#'   \code{data.frame}
#'
#'
#' @export
convertCohortToMatrix <- function(cohort) {
  cohortMatrix <- 
    vapply(cohort[seq_len(nrow(cohort)), seq_len(ncol(cohort) - 2)],
           function(x) as(x, 'numeric'),
           FUN.VALUE=numeric(nrow(cohort))
           )
  rownames(cohortMatrix) <- rownames(cohort)
  return(cohortMatrix)
}














