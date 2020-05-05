#' Draw a forest plot using the statistics calcualted in validation metaestimate
#'
#' @examples
#'
#' @param validationStatsDF A \code{data.frame} containing statistics for Dindex
#'     or concordance index, as available in the \code{list} returned by
#'     `validateMetaEstimateCalculation`.
#' @param stat A \code{character} vector containing the statistic to
#'     use in the forest plot. Options are "dIndex" or "cIndex".
#' @param isSummary A \code{logical} vector indicating which cohorts contain
#'     summary data. Must be one boolean value per row in the data.frame.
#' @param filePath A \code{character} vector containing the file path to write
#'     the plotting results to.
#' @param fileName A \code{character} vector containing the desired file name
#'     with the desired file extension. Supported extensions include .pdf, .png,
#'     .svg and .jpeg, for more information see documentation for the `ggsave()`
#'     function from `ggplot2`.
#' @param ... Additional arguments passed to forestplot.
#'
#' @importFrom scales scientific
#' @importFrom forestplot forestplot
#' @importFrom grid unit grid.grabExpr grid.draw
#' @importFrom ggplot2 ggsave
#' @export
forestPlotMetaEstimate <- function(validationStats, stat, isSummary, filePath,
                                   fileName, ...) {

    # Extract necessary statistics for plotting
    validationStatsDF <- validationStats[[stat]]
    validationStatsDF <- rbind(rep(NA, ncol(validationStatsDF)),
                               validationStatsDF)
    PCOSPscoreList <- validationStats$probabilities
    isSeq <- validationStats$isSequencing
    isSummary <- c(TRUE, isSummary)

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=c("Cohorts", rownames(validationStatsDF)[-1]),
        "pvalue"=c("P value", scientific(validationStatsDF$pval[-1], 2))
    )

    # Extract box sizes
    lengthSeq <- length(unlist(PCOSPscoreList[isSeq]))
    lengthArray <- length(unlist(PCOSPscoreList[!isSeq]))
    boxSizes <- c(NA, vapply(PCOSPscoreList, function(score) length(unlist(score)),
                         FUN.VALUE=numeric(1)),
                  lengthSeq, lengthArray,
                  lengthSeq + lengthArray) / 1000

    # Match correct plot function to call
    if(missing(...)) {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    # Allow user to specify custom parameters
    } else {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes, ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes, ...)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    }
    # Decide whether to plot to device or save to disk
    if (missing(filePath) || missing(fileName)) {
        grid.draw(plot)
    } else {
        grid.draw(plot)
        ggsave(filename=fileName, path=filePath, plot=plot)
    }
}

 #' Draw a forest plot using the statistics calcualted in validation metaestimate
    #'
    #' @examples
    #'
    #' @param validationStatsDF A \code{data.frame} containing statistics for Dindex
    #'     or concordance index, as available in the \code{list} returned by
    #'     `validateMetaEstimateCalculation`.
    #' @param stat A \code{character} vector containing the statistic to
    #'     use in the forest plot. Options are "dIndex" or "cIndex".
    #' @param isSummary A \code{logical} vector indicating which cohorts contain
    #'     summary data. Must be one boolean value per row in the data.frame.
    #' @param filePath A \code{character} vector containing the file path to write
#'     the plotting results to.
#' @param fileName A \code{character} vector containing the desired file name
#'     with the desired file extension. Supported extensions include .pdf, .png,
#'     .svg and .jpeg, for more information see documentation for the `ggsave()`
#'     function from `ggplot2`.
#' @param ... Additional arguments passed to forestplot.
#'
#' @importFrom scales scientific
#' @importFrom forestplot forestplot
#' @importFrom grid unit grid.grabExpr grid.draw
#' @importFrom ggplot2 ggsave
#' @export
forestPlotModelComparision <- function(clinicalModelStats, stat, isSummary, filePath,
                                   fileName, ...) {

    indexClinical <- as.matrix(clinicalModelStats$clinical[[stat]])
    indexPCOSP <- as.matrix(clinicalModelStats$PCOSP[[stat]])

    spacing <- rbind(rep(NA, 3), rep(NA, 3))
    plotData <- matrix(nrow=0, ncol=3)
    for (i in seq_len(nrow(indexClinical))) {
        plotData <-
            rbind(
                plotData,
                spacing,
                indexClinical[i, 1:3],
                indexPCOSP[i, 1:3]
        )
    }

    isSummary <- c(TRUE, isSummary)

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=c("Cohorts", rownames(validationStatsDF)[-1]),
        "pvalue"=c("P value", scientific(validationStatsDF$pval[-1], 2))
    )

    # Extract box sizes
    lengthSeq <- length(unlist(PCOSPscoreList[isSeq]))
    lengthArray <- length(unlist(PCOSPscoreList[!isSeq]))
    boxSizes <- c(NA, vapply(PCOSPscoreList, function(score) length(unlist(score)),
                             FUN.VALUE=numeric(1)),
                  lengthSeq, lengthArray,
                  lengthSeq + lengthArray) / 1000

    # Match correct plot function to call
    if(missing(...)) {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
        # Allow user to specify custom parameters
    } else {
        if (stat == "dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes, ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      isSummary, boxSizes, ...)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    }
    # Decide whether to plot to device or save to disk
    if (missing(filePath) || missing(fileName)) {
        grid.draw(plot)
    } else {
        grid.draw(plot)
        ggsave(filename=fileName, path=filePath, plot=plot)
    }
}


.forestPlotDindex <- function(labelText, validationStatsDF, isSeq, isSummary, boxSizes,
                              ...) {

    if (!missing(...)) {
        plot <- grid.grabExpr(forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        normFun <- local({
            i = 0
            isSeq=isSeq
            b_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            l_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        sumFun <- local({
            i = 0
            s_clrs =c("#FF7F00","#1F78B4","grey57")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <- grid.grabExpr(forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               xlab="Log2 HR",
                               is.summary=isSummary,
                               clip=c(-1, 2.5),
                               txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                              ticks=gpar(cex=0.8),
                                              xlab=gpar(fontfamily="Helvetica",
                                                        cex=1)),
                               col=fpColors(box="black"),
                               title=" ",
                               new_page=FALSE,
                               fn.ci_norm=normFun,
                               fn.ci_sum=sumFun,
                               zero=0,
                               graphwidth=unit(2, "inches"),
                               align=c("l"),
                               pch=16,
                               boxsize=boxSizes + 0.2))
    }
    return(plot)
}

.forestPlotCindex <- function(labelText, validationStatsDF, isSeq, isSummary, boxSizes,
                              ...) {
    if (!missing(...)) {
       plot <- grid.grabExpr(forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               ...))
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        normFun <- local({
            i = 0
            isSeq=isSeq
            b_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            l_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        sumFun <- local({
            i = 0
            s_clrs =c("#FF7F00","#1F78B4","grey57")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(...,col=s_clrs[i])
            }
        })
        # Make the plot
        plot <- grid.grabExpr(forestplot::forestplot(labelText,
                                                     validationStatsDF[, c("mean", "lower", "upper")],
                                                     xlab="C-index",
                                                     is.summary=isSummary,
                                                     clip=c(0.3, 0.8),
                                                     txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                                                    ticks=gpar(cex=0.9),
                                                                    xlab=gpar(fontfamily="Helvetica",
                                                                              cex=1)),
                                                     col=fpColors(box="black"),
                                                     title=" ",
                                                     new_page=FALSE,
                                                     fn.ci_norm=normFun,
                                                     fn.ci_sum=sumFun,
                                                     zero=0.5,
                                                     graphwidth=unit(2, "inches"),
                                                     align=c("l"),
                                                     pch=16,
                                                     boxsize=boxSizes + 0.2))
    }
    return(plot)
}
