#' Read in .tsv files of expression data and calculate then return a list with
#'   a signature score for each expression file.
#'
#' @examples
#'
#' @param dataDir A \code{character} vector containing the path to the
#'    data directory. We recommend using `file.path` to make your paths
#'    platform agnostic!
#' @param fileNames A \code{character} vector of file names (with extensions)
#'    to read from and extract the gene signature scores.
#' @param geneCoefFile A \code{character} vector with the file name of the
#'    gene coefficients.
#' @param ... Fallthrough parameter to `sig.scores` function from `genefu`.
#'
#' @return A \code{list} of signature scores per specified file, named
#'    based on the file name.
#'
#' @importFrom genefu sig.score
#' @export
parseGenefuFiles <- function(dataDir, fileNames, geneCoefFile, ...) {

    geneCoefficients <- read.table(file.path(dataDir, geneCoefFile),
                                   sep="\t", header=TRUE)

    expressionData <- lapply(fileNames,
                             function(file, dataDir)
                                 read.table(file.path(dataDir, file),
                                            sep="\t", header=T),
                             dataDir=dataDir)

    if (!missing(...)) {
        sigScores <- lapply(expressionData,
                            function(data, geneCoef, ...)
                                sig.score(x=geneCoef, data=data, ...)$score,
                            geneCoef=geneCoefficients,
                            ...=...
                            )
    } else {
        sigScores <- lapply(expressionData,
                            function(data, geneCoef)
                                sig.score(x=geneCoef, data=data,
                                          do.mapping=FALSE, mapping,
                                          )$score,
                            geneCoef=geneCoefficients)
    }
    names(sigScores) <- gsub("\\_.*txt|.txt", "", fileNames)
    return(sigScores)
}