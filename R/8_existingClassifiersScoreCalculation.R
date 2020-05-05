######### Birnbaum
library(genefu)
gene_coeff=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/bmc_genes.txt", sep="\t", header =T )

#' Read in .tsv files of expression data and calculate then return a list with
#'   a signature score for each expression file.
#'
#' @examples
#'
#' @param dataDir
#' @param fileNames
#' @param geneCoefFile
#' @param ...
#'
#' @importFrom genefu sig.score
#' @export
parseGeneFuFiles <- function(dataDir, fileNames, geneCoefFile, ...) {

    geneCoefficients <- read.table(file.path(dataDir, geneCoefFile),
                                   sep="\t", header=TRUE)

    expressionData <- lapply(fileNames,
                             function(file)
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
    return(sigScores)
}





