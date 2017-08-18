xx= read.table("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Pathway_analysis/PANTHER/panther-biological.txt",sep="\t",header=T)
df=data.frame(Biological_Processes=xx[,1],No_of_genes=xx[,4])
pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure5a.pdf")

BB=barplot(df$No_of_genes, main="Biological Processes", names.arg=df$Biological_Processes,las=2,
        col=c(  "#FDBF6F", "#FF7F00", "#E31A1C","#F1E2CC","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E6F5C9" ,"grey40","#FFF2AE" , "#CCCCCC" , "gold"),
        xlab="No. of genes", border = NA,space=0,cex.lab = 1, cex.axis = 1, cex.names = 0.6, horiz = TRUE)


text(df$No_of_genes,BB, label = round(xx[,8],2),  pos=4,cex = 0.7, col = "black")
legend("bottomright", legend= xx[,2], fill =c(  "#FDBF6F", "#FF7F00", "#E31A1C","#F1E2CC","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E6F5C9" ,"grey40","#FFF2AE" , "#CCCCCC" , "gold"), bty='n' )



########################


xx= read.table("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Pathway_analysis/PANTHER/molecular_function.txt",sep="\t",header=T)
df=data.frame(Molecular_function=xx[,2],No_of_genes=xx[,4])

pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/Figure5b.pdf")

BB=barplot(df$No_of_genes, main="Molecular function", names.arg=df$Molecular_function,las=2,
           col=c( "#8DD3C7" ,"#FFFFB3" ,"cyan" ,"#FB8072" ,"magenta" ,"#FDB462", "#B3DE69" ,"#FCCDE5",  "#8DA0CB" ),
           xlab="No. of genes", border = NA,space=0,cex.lab = 1, cex.axis = 1, cex.names = 0.6, horiz = TRUE)


text(df$No_of_genes,BB, label = round(xx[,8],2),  pos=4,cex = 0.7, col = "black")
legend("bottomright", legend= xx[,1], fill =c( "#8DD3C7" ,"#FFFFB3" ,"cyan" ,"#FB8072" ,"magenta" ,"#FDB462", "#B3DE69" ,"#FCCDE5",  "#8DA0CB" ), bty='n' )



