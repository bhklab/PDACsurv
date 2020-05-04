#' Compare clinical models between patient cohorts
#'
## TODO:: HEEWON - description
#'
#' @param clinicalFeatures
#' @param cohortClasses
#' @param cohorts
#' @param models A \code{character} vector names for cohorts in `clinicalFeatures`
#'     to fit a linear model for. All models will be compared against each
#'     cohort in `cohorts`, or all cohorts if `cohorts` is not specified.
#' @param formula
#'
#' @return A \code{list} with the first level representing the model the
#'   data was fit to, the second level representing the cohort the model was
#'   compared to, and the third level the associated statistics for that
#'   comparison.
#'
#' @importFrom verification
#' @importFrom glmnet glm
#' @importFrom stats anova
#' @export
compareClinicalModels <- function(clinicalFeatures, cohortClasses, cohorts,
                                  models, formula="binary_grp ~ Age + Sex +
                                  T_status + N + M + Grade")
{
    modelFits <- summarizeClinicalModels(clinicalFeatures, cohorts=models)

    fitModels <- lapply(modelFits, function(cohort) cohort$model)

    cFeatures <- clinicalFeatures[grep(paste(cohorts, collapse="|"),
                                       names(clinicalFeatures))]

    cClasses <- cohortClasses[grep(paste(cohorts, collapse="|"),
                                   names(cohortClasses))]

    cohortProbs <- lapply(fitModels,
                          function(model, clinicalCohorts)
                              .predictProbPerCohort(clinicalCohorts, model),
                          clinicalCohorts=cFeatures)

    modelAUCs <- .estimateModelAUCs(cFeatures, cohortProbs)

    return(modelAUCs)
}

#' Fit a generalized linear model for each clinical cohort in a list
#'
#' @param clinicalFeatures A list of \code{data.frame}s containing the clinical
#'     features for each cohort
#' @param cohorts A \code{character} vector of cohorts to subset from
#'      the clinicalFeatures argument. If excluded all cohorts are used.
#' @param formula A \code{character} vector containing the model formula
#'      expressed as a string. If excluded defaults to
#'      "binary_grp ~ Age + Sex + T_status +N + M + Grade".
#' @param ... Fallthrough parameters to `glm` function. If used, each cohort
#'      is passed as the data argument and formula as the formula. All other
#'      parameters must be specified in `...`.
#'
#' @return A named \code{list} containing the fitted model, model summary and anova
#'     for each cohort (named by `cohorts` argument or as
#'     `names(clincalFeatures)` if cohorts is absent).
#'
#' @importFrom glmnet glm
#' @importFrom stats as.formula
#' @export
summarizeClinicalModels <- function(clinicalFeatures, cohorts,
                                    formula="binary_grp ~ Age + Sex +
                                    T_status + N + M + Grade",
                                    ...)
{

    if (!missing(cohorts)) {
        cFeatures <- clinicalFeatures[grep(paste(cohorts, collapse="|"),
                                           names(clinicalFeatures))]
        names(cFeatures) <- cohorts
    } else {
        cFeatures <- clinicalFeatures
    }

    ##TODO:: Add error handling check for formula with data columns

    if (!missing(...)) {
        models <- lapply(cFeatures,
                          function(data, formula)
                              glm(as.formula(formula), data=data, ...),
                          formula=formula)
    } else {
        models <- lapply(cFeatures,
                          function(data, formula)
                              glm(as.formula(formula), data=data,
                                  na.action=na.exclude,
                                  family=binomial(link="logit")),
                          formula=formula)
    }

    summaries <- lapply(models, summary)
    anovas <- lapply(models, summary)

    results <- .zipLists(models, summaries, anovas,
                         sublistNames=c("model", "summary", "anova"))
    return(results)
}

#'
#'
#'
#'
.predictProbPerCohort <- function(clinicalCohorts, model) {
    lapply(clinicalCohorts, .predictClinicalCohortProb, model=model)
}

#'
#'
#'
#'
.predictClinicalCohortProb <- function(clinicalCohort, model) {

    rownames(clinicalCohort) <- clinicalCohort$ID

    prediction1 <- predict(model, clinicalCohort, na.action=na.exclude,
                           type="response")

    prediction2 <- 1 - clinicalCohort$pred_prob[which(clinicalCohort$ID %in%
                                                          names(prediction1))]
    return(list(
        "pred1"=prediction1,
        "pred2"=prediction2
    ))
}

#'
#'
#'
#'
.estimateModelAUCs <- function (cFeatures, cohortProbs) {

    res <- list()
    for (i in seq_along(cohortProbs)) {
        res[[i]] <- structure(lapply(seq_along(cFeatures),
               function(j, cFeatures, cohortProbs)
                   .calcModelAUC(cFeatures[[j]], cohortProbs[[i]][[j]]),
               cFeatures=cFeatures,
               cohortProbs=cohortProbs),
               .Names=names(cFeatures))
    }
    return(structure(res, .Names=names(cohortProbs)))

}

.calcModelAUC <- function(coh, preds) {

    rownames(coh) <- coh$ID

    list(
        "roc1"=list(
            "roc"=reportROC(coh[names(preds$pred1),]$binary_grp, preds$pred1, plot=FALSE),
            "pval"=roc.area(coh[names(preds$pred1),]$binary_grp, preds$pred1)$p.value
        ),
        "roc2"=list(
            "roc"=reportROC(coh[names(preds$pred1),]$binary_grp, preds$pred2, plot=FALSE),
            "pval"=roc.area(coh[names(preds$pred1),]$binary_grp, preds$pred2)$p.value
        )
    )
}

#'
#'
#'
#'
#' @export
barplotModelAUCs <- function(modelAUCs, model, names, colours, ...) {
    data <-data.frame(
        lapply(modelAUCs[[model]],
               function(cohort)
                   as.numeric(c(cohort$roc1$roc$AUC, cohort$roc2$roc$AUC))
        ))
    colnames(data) <- names
    rownames(data) <- c("Clinicopathological mode", "PCOSP")

    if (!missing(...)) {
        barplot(as.matrix(data), ...)
    } else {
        barplot(as.matrix(data), main="", ylim=c(0, 0.8), ylab="AUCs",
                beside=TRUE, col=colours, cex.main=1.4,
                space=rep(c(0.6, 0.08), length(modelAUCs[[model]])))
    }
}