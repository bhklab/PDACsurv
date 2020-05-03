
#' Convert the levels of a factor to number
#'
#' A convenience function to converting factor levels into numeric
#'
#' @param factor
#'
##FIXME:: Make sure all the data is either a factor, or not a factor?
##FIXME:: This breaks when written with levels because not all of them are factors
as.numeric.factor <- function(factor) { as.numeric(as.character(factor)) }

    #' Exclude samples censored before year 1
#'
#' @param seqCohort A \code{}
#'
#' @return A sorted vector of indices
whichNotCensoredYearOne <- function(seqCohort) {
    idxNotCensored <- which(as.numeric.factor(seqCohort$OS) <= 365 &
                                as.numeric.factor(seqCohort$OS_Status) == 1)
    idxNotYearOne <- which(as.numeric.factor(seqCohort$OS) > 365)
    return(sort(c(idxNotCensored, idxNotYearOne)))
}

#' Merge common samples not censored before year 1
#'
#' @param seqCohort A \code{matrix} containing expression data for the sequence
#'   cohort
#' @param arrayCohort A \code{matrix} containing expression data for the array
#'   cohort
#'
#' @return A \code{data.frame} containing the merged cohorts
#'
#' @export
mergeCommonData <- function(seqCohort, arrayCohort) {
    if (!all(rownames(arrayCohort) == rownames(seqCohort))) {
        stop("Rownames do not match between the seqCohort and arrayCohort!")
    }

    notCensoredYearOne <- whichNotCensoredYearOne(seqCohort)
    return(rbind(seqCohort[notCensoredYearOne, ],
                 arrayCohort[notCensoredYearOne, ])
    )
}

#' Extract all cohorts from a list of cohorts
#'
#' Input a named list of cohorts to be modelled and assign each cohort to an
#'   environment (the global environment is recommended)
#'
#' @param cohortList A named \code{list} containing the cohorts for an analysis
#' @param environment The \code{environment} into which the variables will be
#'   assigned. Use \code{globalenv()} to assign them to the global environment
#'   (recommended).
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
               function(x) as.numeric.factor(x),
               FUN.VALUE=numeric(nrow(cohort))
        )
    rownames(cohortMatrix) <- rownames(cohort)
    return(cohortMatrix)
}

#'
#'
#'
#' @importFrom switchBox SWAP.KTSP.Classify
#' @importFrom pROC reportROC
.predictKTSP <- function(formattedValCohort, selectedModels, nthread){

    # Temporily change number of cores to parallelize over
    opts <- options()
    options("mc.cores"=nthread)
    on.exit(options(opts))

    grpIdx <- formattedValCohort$grpIndex

    predictions <- bplapply(selectedModels,
                          function(model, valMat)
                              SWAP.KTSP.Classify(t(valMat), model),
                          valMat=formattedValCohort$mat[grpIdx,])

    aucReports <- suppressMessages(lapply(predictions,
                        function(pred, group)
                            reportROC(group,
                                      as.numeric.factor(pred),
                                      plot=FALSE),
                        group=formattedValCohort$grp))

    return(list(
        "predictions"=predictions,
        "AUCs"=vapply(aucReports, function(rep) as.numeric(rep$AUC), FUN.VALUE=numeric(1)),
        "aucSEs"=vapply(aucReports, function(rep) as.numeric(rep$AUC.SE), FUN.VALUE=numeric(1))
    ))
}








