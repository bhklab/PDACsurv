#' Density plot AUC distribution of a PCOSP model
#'
#' Create a density plot of the distrubtion of AUCs between the overall,
#'    sequencing and microarray cohorts in a PCOSP model.
#'
#' @param formattedValCohorts A
#' @param selectedModels A
#' @param seqCohorts A
#' @param vlines A
#' @param nthread A
#' @param title A
#' @param filePath A
#' @param fileName A
#'
#' @return A \code{grob} objects from \code{ggplot2}
#'
#' @export
#' @import ggplot2
densityPlotModel <- function(formattedValCohorts, selectedModels, seqCohorts, title,
                             vlines, nthread, filePath, fileName) {

    dStats <- .densityStats(formattedValCohorts, selectedModels, seqCohorts,
                            nthread)

    densityDF <- data.frame(
        "platforms"=rep(c("1. Array-based", "2. Sequencing", "3. Overall"), each=length(selectedModels)),
        "BAC"=unlist(c(dStats$arrayAUCs, dStats$seqAUCs, dStats$metaAUCs))
    )
    vlineDF <- data.frame(
        "platforms"=unique(densityDF$platforms), v1=vlines
    )

    plt <- ggplot(densityDF, aes(x=BAC, color=platforms)) +
        geom_density(aes(fill=platforms), alpha=0.5) +
        geom_vline(data=vlineDF,
                   aes(xintercept=v1,
                       color=platforms),
                   linetype="dashed",
                   size=1.5) +
        labs(x = "Balanced Accuracy", y="Density")+
        ggtitle(title) +
        facet_wrap(~ platforms, ncol=1, strip.position="right") +
        theme(plot.title=element_text(hjust = 0.5),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "grey"),
              legend.position="none")

    if (!missing(filePath) && !missing(fileName)) {
        ggsave(filename=fileName, path=filePath, plot=plt)
    }
    plt
}


#'
#'
#'
#'
#'
#'
.densityStats <- function(formattedValCohorts, selectedModels, seqCohorts,
                          nthread) {
    KTSPs <- predictKTSPs(formattedValCohorts, selectedModels, nthread)

    isSeq <- names(KTSPs) %in% seqCohorts

    nModels <- length(selectedModels)

    metaAUCs <- .summarizeKTSPs(KTSPs, nModels, nthread)
    seqAUCs <- .summarizeKTSPs(KTSPs[isSeq], nModels, nthread)
    arrayAUCs <- .summarizeKTSPs(KTSPs[!isSeq], nModels, nthread)

    return(list(
        "metaAUCs"=metaAUCs,
        "seqAUCs"=seqAUCs,
        "arrayAUCs"=arrayAUCs
    ))

}

#'
#'
#'
#'
#'
predictKTSPs <- function(formattedValCohorts, selectedModels, nthread) {
    lapply(formattedValCohorts,
           function(cohort, models, nthread) .predictKTSP(cohort, models, nthread),
           models=selectedModels,
           nthread=nthread)
}

#'
#'
#'
#' @importFrom BiocParallel bplapply
.summarizeKTSPs <- function(KTSPs, nModels, nthread) {
    # Temporily change number of cores to parallelize over
    opts <- options()
    options("mc.cores"=nthread)
    on.exit(options(opts))

    AUCs <- do.call(rbind, lapply(KTSPs, function(KTSP) KTSP$AUCs))
    aucSEs <- do.call(rbind, lapply(KTSPs, function(KTSP) KTSP$aucSEs))

    bplapply(seq_len(nModels),
             function(i, AUCs, aucSEs)
                 combine.est(
                     AUCs[, i],
                     aucSEs[, i],
                     na.rm=TRUE
                 )$estimate,
             AUCs=AUCs,
             aucSEs=aucSEs)
}