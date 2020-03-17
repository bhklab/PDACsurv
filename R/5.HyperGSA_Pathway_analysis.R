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




