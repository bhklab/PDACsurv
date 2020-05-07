#' Convert the levels of a factor to number
#'
#' A convenience function to converting factor levels into numeric
#'
#' @param factor A
#'
#' @export
##FIXME:: Make sure all the data is either a factor, or not a factor?
##FIXME:: This breaks when written with levels because not all of them are factors
as.numeric.factor <- function(factor) { as.numeric(as.character(factor)) }

#' Exclude samples censored before year 1
#'
#' @param seqCohort A \code{}
#'
#' @return A sorted vector of indices
#' @export
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
#' @importFrom reportROC reportROC
#' @importFrom BiocParallel bplapply
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

##TODO:: Should this go in utilties?
#' Merge n lists element-wise into a list of sublists n items long
#'
#' Take n lists and combine them element-wise into a single list with each sublist
#'   containing the corresponding elements of the original lists. Assigns
#'   sublistNames as the names of the new zipped sublists.
#'
#' @param ... Two or more \code{list}s to zip together element-wise. New list
#'     is named using first list.
#' @param sublistNames A \code{character} vector n items long, specifying the
#'     labels for the zipped items
#'
#' @keywords internal
.zipLists <- function(..., sublistNames) {

    ## FIXME:: Error checking
    # lengths <- vapply(as.list(...), length, FUN.VALUE=numeric(1))
    # if (!all(lengths==length(sublistNames)))
    #     stop("Please ensure sublistNames is the same length as the lists you
    #          are zipping!")

    # Merge lists into sublists element-wise
    zipped <- mapply(list, ..., SIMPLIFY=FALSE)

    # Name each zipped item in the sublist
    zipped <- lapply(zipped, function(stat) structure(stat, .Names=sublistNames))
    return(zipped)
}


#'
#'
#'
#'
#' @export
subsetSharedCohortsAndSamples <- function(metaestimateData) {
        # Intersect the cohort names
        cohortNames <- Reduce(intersect, lapply(metaestimateData, names))
        names(cohortNames) <- cohortNames

        # Subset and reorder all data
        metaestimateData <- lapply(metaestimateData, function(data, cohorts) data[cohorts], cohorts=cohortNames)

    classiferSampleNames <- c(lapply(metaestimateData[1:4], function(data) lapply(data, names)),
                              lapply(metaestimateData[5], function(data) lapply(data, rownames)))

    sharedSampleNames <- lapply(cohortNames,
                                function(cohort, classifierSampleNames)
                                    Reduce(intersect, lapply(classiferSampleNames,
                                                             function(class, cohort) class[[cohort]],
                                                             cohort=cohort)),
                                class=classiferSampleNames)
    # Subset to shared samples
    structure(lapply(names(metaestimateData),
           function(class, data, samples) structure(lapply(names(samples),
                                                function(cohort, data, samples, class) {
                                                    if (class=="survival") {
                                                        data[[class]][[cohort]][samples[[cohort]], ]
                                                    } else {
                                                        data[[class]][[cohort]][samples[[cohort]]]
                                                    }
                                                },
                                                data=data, samples=samples, class=class),
                                                .Names=names(samples)),
           data=metaestimateData, samples=sharedSampleNames),
           .Names=names(metaestimateData))
}






