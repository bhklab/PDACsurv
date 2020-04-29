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
forestPlotMetaEstimate <- function(validationStats, stat, isSummary,
                                   filePath, fileName, ...) {

    # Extract necessary statistics for plotting
    validationStatsDF <- validationStats[[stat]]
    PCOSPscoreList <- validationStats$PCOSPscores
    isSeq <- validationStats$isSequencing

    # Construct the forest plot table
    labelText <- data.frame(
        "cohort"=rownames(validationStatsDF),
        "pvalue"=scientific(validationStatsDF$pval, 2)
    )

    # Extract box sizes
    lengthSeq <- length(unlist(PCOSPscoreList[isSeq]))
    lengthArray <- length(unlist(PCOSPscoreList[!isSeq]))
    boxSizes <- c(vapply(PCOSPscoreList, function(score) length(unlist(score)),
                         FUN.VALUE=numeric(1)),
                  lengthSeq, lengthArray,
                  lengthSeq + lengthArray) / 1000

    # Match correct plot function to call
    if(missing(...)) {
        if (stat=="dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq,
                                      boxSizes)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq,
                                      boxSizes)
        } else {
            stop(paste0("There is no statistic called: ", stat))
        }
    # Allow user to specify custom parameters
    } else {
        if (stat=="dIndex") {
            plot <- .forestPlotDindex(labelText, validationStatsDF, isSeq, boxSizes,
                              ...)
        } else if (stat=="cIndex") {
            plot <- .forestPlotCindex(labelText, validationStatsDF, isSeq, boxSizes,
                              ...)
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


.forestPlotDindex <- function(labelText, validationStatsDF, isSeq, boxSizes,
                              ...) {

    if (!missing(...)) {
        forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               ...)
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        normFun <- local({
            i = 0
            isSeq=isSeq
            b_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            l_clrs=ifelse(isSeq, "#FF7F00", "#1F78B4")
            #s_clrs =c(rep("red",10),"green","pink","yellow","orange")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
                #fpDrawSummaryCI(...,col=s_clrs[i])
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
                               clip=c(-1, 1.25),
                               txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                              ticks=gpar(cex=0.8),
                                              xlab=gpar(fontfamily="Helvetica",
                                                        cex=1)),
                               col=fpColors(box="black"),
                               title=" ",
                               new_page=FALSE,
                               fn.ci_norm=normFun, fn.ci_sum=sumFun,
                               zero=0,
                               graphwidth=unit(2, "inches"),
                               align=c("l"),
                               pch=16,
                               boxsize=boxSizes + 0.2))
        return(plot)
    }
}

.forestPlotCindex <- function(label, validationStatsDF, isSummary, boxSizes,
                              ...) {
    if (!missing(...)) {
        forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               ...)
    } else {
        # Set plot colouring functions
        # ##TODO:: Determine if there is a more readable way to write this?
        .normFun <-  local({
            i = 0
            boxClrs <- ifelse(isSummary, "orange", "blue")
            lineClrs <- ifelse(isSummary, "orange", "blue")
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line=lineClrs[i], clr.marker=boxClrs[i])
            }
        })
        .sumFun <- local({
            i = 0
            summaryClrs =c("orange", "blue", "grey")
            function(..., col){
                i <<- i + 1
                fpDrawSummaryCI(..., col=summaryClrs[i])
            }
        })
        # Make the plot
        plot <- grid.grabExpr(forestplot::forestplot(labelText,
                               validationStatsDF[, c("mean", "lower", "upper")],
                               xlab="Log2 HR",
                               is.summary=isSummary,
                               clip=c(-1, 1.25),
                               txt_gp=fpTxtGp(label=gpar(fontfamily="Helvetica"),
                                              ticks=gpar(cex=0.8),
                                              xlab=gpar(fontfamily="Helvetica",
                                                        cex=1)),
                               col=fpColors(text="black"),
                               title=" ",
                               new_page=FALSE,
                               fun.ci_norm=.normFun,
                               fun.ci_sum=.sumFun,
                               zero=0,
                               graphwidth=unit(2, "inches"),
                               align=c("l"),
                               pch=16,
                               boxsize=boxSizes + 0.2))
        return(plot)
    }
}
