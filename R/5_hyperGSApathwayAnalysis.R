library(piano)
reference_genes=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/reference_genes.txt")
pcosp_genes=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/pcosp_genes.txt")

##### HALLMARK PATHWAYS
hallmark= loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/h.all.v6.1.symbols.gmt")
results_hallmark=runGSAhyper(pcosp_genes,gsc= hallmark,  adjMethod="fdr",  universe= reference_genes)
results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/Users/vandanasandhu/Desktop/kegg.txt")

#####  Conanical pathway
genesets=loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/c2.cp.v6.1.symbols.gmt", type="auto")
results=runGSAhyper(pcosp_genes,gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$resTab[,2]<0.05),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/Users/vandanasandhu/Desktop/conanical.txt")

##### GO Molecular function
genesets=loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/c5.mf.v6.1.symbols.gmt", type="auto")
results=runGSAhyper(pcosp_genes,gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$resTab[,2]<0.05),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/Users/vandanasandhu/Desktop/conanical.txt")


##### GO Cellular component
genesets=loadGSC("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Gene-sets/c5.cc.v6.1.symbols.gmt", type="auto")
results=runGSAhyper(pcosp_genes,gsc= genesets,  adjMethod="fdr", universe= reference_genes)
results$resTab[which(results$resTab[,2]<0.05),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/Users/vandanasandhu/Desktop/conanical.txt")

#' Fit GSEA model to Gene Set Collection
#'
##TODO:: HEEWON - documentation
#'
#' @examples
#' data(referenceGenes)
#' data(PCOSPgenes)
#' results <- fitGSEAtoGeneSetCollection()
#'
#' @param GSC A \code{gene set collection} of data to fit the GSEA model to
#' @param referenceGenes A \code{data.frame} containing the refence genes for the
#'   GSEA model
#' @param PCOSPgenes A \code{data.frame} containing the PCSOP genes for the
#'   GSEA model
#' @param filePath A \code{character} vector containing the path to
#'   the cohort data in .gmt format. We reccomend using `file.path()` function
#'   to construct this so that scripts work cross platform. If missing the `data`
#'   is used for the calculation.
#'
#' @return A \code{data.frame} of significant results from the model fit
#'
#' @importFrom piano loadGSC runGSAhyper
#' @export
##TODO:: Add ... to allow setting piano method parameters
fitGSEAtoGeneSetCollection <- function(GSC, referenceGenes, PCOSPgenes, filePath)
  {
  if (missing(data) && !missing(filePath)) {
    geneSetCollection <- loadGSC(file.path(filePath))
  } else if (!missing(data) && missing(filePath)) {
    geneSetCollection <- data
  } else {
    stop("Please either pass in a data object or specific a file path from
         which to load the data!")
  }
  results <- runGSAhyper(PCOSPgenes, geneSetCollection)
  results <- results[which(results$resTab[, 2] > 0.05), ]
}
