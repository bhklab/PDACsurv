load("/Users/vandanasandhu/Desktop/os_predictor_project/pdac_classifier.RData")

############# HEATMAPS ##################################
library("foreach")
source("/Users/vandanasandhu/Documents/Scripts/clustermap.R") # Make package functions available
tmp = blocks()         # Contains the 110x40-matrix X and 40-vectors y1,y2,y3


################ Subtyping PDACs based on Collisson classification ################################################################# 
xx<- combine_list$collisson

c<-colnames(xx)
r<-rownames(xx)


mydatamatrix <- data.matrix(xx) 
mydatascale <- t(scale(t(mydatamatrix))) 


X= mydatascale
plot.init(tree=c(1,3), cbar=2)
hcluster(X, clust="row", distance="spearman", linkage="complete",ncomp=5)
hcluster(X, clust="col", distance="spearman", linkage="complete",ncomp=5)
subclust(3, clust="col",B=100, method="part", min.size=5)
subclust(NA, clust="row",B=100, method="part", min.size=5)

plot.hmap(X,colorscale="blue-white-red")
legend("topright",col=c("Blue","Red","green"), c("Exocrine","QM","Classical"),lty=1,cex=0.5)   
plot.text(rownames(X), cex=0.5, side=4)
plot.text(colnames(X), cex=0.2, side=1)
plot.tree(side=3)
plot.hmap.key()
dev.off()
gc <- get.subclust(clust="col", labels=c, order="tree") #clust="col" er default Evt clust="row"
gc=gc[order(gc$label),]
coll=gc

ss=combine_list$survival

mf <- Surv(as.numeric(as.character(ss[,2])),as.numeric(as.character( ss[,3]))== 1)
class(mf)

fit.mf<- survfit(mf ~ gc$cluster)
sd=survdiff(mf ~ gc$cluster)
p.val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p.val
plot(fit.mf,col=c("Blue","green","Red"),xlab="Days after surgery (days)",ylab="Survival")
legend(3500,0.6,col=c("Blue","green","Red"), c("Exocrine","Classical","QM-PDA"),lty=1,cex=0.5)   
legend("topright",paste("P=",round(p.val, digit=2), sep="") )  


################ Subtyping PDACs based on Moffitt classification ################################################################# 

xx<- combine_list$moffitt
c<-colnames(xx)
r<-rownames(xx)


mydatamatrix <- data.matrix(xx) 
mydatascale <- t(scale(t(mydatamatrix))) 


X= mydatascale
plot.init(tree=c(1,3), cbar=2)

hcluster(X, clust="row", distance="spearman", linkage="complete",ncomp=5)
hcluster(X, clust="col", distance="spearman", linkage="complete",ncomp=5)
subclust(2, clust="col",B=100, method="part", min.size=5)
subclust(NA, clust="row",B=100, method="part", min.size=5)

plot.hmap(X,colorscale="darkblue-white-red")
legend("topright",col=c("Red","green"), c("Basal","Classical"),lty=1,cex=0.5)   

plot.text(rownames(X), cex=0.5, side=4)
plot.text(colnames(X), cex=0.2, side=1)
plot.tree(side=3)
plot.hmap.key()

gc <- get.subclust(clust="col", labels=c, order="tree") #clust="col" er default Evt clust="row"
gc=gc[order(gc$label),]

moff=gc
ss=combine_list$survival

mf <- Surv(as.numeric(as.character(ss[,2])),as.numeric(as.character( ss[,3]))== 1)
class(mf)

fit.mf<- survfit(mf ~ gc$cluster)
sd=survdiff(mf ~ gc$cluster)
p.val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p.val
plot(fit.mf,col=c("Red","green"),xlab="Days after surgery (days)",ylab="Survival")

legend(3500,0.6,col=c("Red","green"), c("Basal","Classical"),lty=1,cex=0.5)   
legend("topright",paste("P=",round(p.val, digit=2), sep="") )  




################ Subtyping PDACs based on Bailey classification ################################################################# 

xx<- combine_list$bailey
c<-colnames(xx)
r<-rownames(xx)


mydatamatrix <- data.matrix(xx) 
mydatascale <- t(scale(t(mydatamatrix))) 


X= mydatascale
plot.init(tree=c(1,3), cbar=2)

hcluster(X, clust="row", distance="spearman", linkage="complete",ncomp=5)
hcluster(X, clust="col", distance="spearman", linkage="complete",ncomp=5)
subclust(4, clust="col",B=100, method="part", min.size=5)
subclust(4, clust="row",B=100, method="part", min.size=5)

plot.hmap(X,colorscale="darkblue-white-red")
legend("topright",col=c("Red","green","blue","cyan"), c("ADEX","Pancreatic-Progenitor","Squamous","Immunogenic"),lty=1,cex=0.5)   
plot.text(rownames(X), cex=0.5, side=4)
plot.text(colnames(X), cex=0.2, side=1)
plot.tree(side=3)
plot.tree(side=2, lwd=3)

plot.hmap.key()
gc <- get.subclust(clust="col", labels=c, order="tree") #clust="col" er default Evt clust="row"
gc=gc[order(gc$label),]
bal=gc

ss=combine_list$survival

mf <- Surv(as.numeric(as.character(ss[,2])),as.numeric(as.character( ss[,3]))== 1)
class(mf)

fit.mf<- survfit(mf ~ gc$cluster)
sd=survdiff(mf ~ gc$cluster)

p.val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p.val
plot(fit.mf,col=c("Red","green","Blue","cyan"),xlab="Days after surgery (days)",ylab="Survival")

legend(3500,0.6,col=c("Red","green","Blue","cyan"), c("ADEX","Pancreatic-Progenitor","Squamous","Immunogenic"),lty=1,cex=0.5)   
legend("topright",paste("P=",round(p.val, digit=2), sep="") )  

################ Subtyping PDACs based on PAM50 classification ################################################################# 

xx<- combine_list$pam50
c<-colnames(xx)
r<-rownames(xx)

mydatamatrix <- data.matrix(xx) 
mydatascale <- t(scale(t(mydatamatrix))) 


X= mydatascale
plot.init(tree=c(1,3), cbar=2)

hcluster(X, clust="none", distance="spearman", linkage="complete",ncomp=5)
hcluster(X, clust="none", distance="spearman", linkage="complete",ncomp=5)
subclust(2, clust="col",B=100, method="part", min.size=5)

plot.hmap(X,colorscale="darkblue-white-red")
legend("topright",col=c("Red","green"), c("Basal","Classical"),lty=1,cex=0.5)   

plot.text(rownames(X), cex=0.5, side=4)
plot.text(colnames(X), cex=0.2, side=1)
plot.tree(side=3)

plot.hmap.key()
dev.off()
gc <- get.subclust(clust="col", labels=c, order="tree") #clust="col" er default Evt clust="row"
gc=gc[order(gc$label),]
pam=gc

ss=combine_list$survival

mf <- Surv(as.numeric(as.character(ss[,2])),as.numeric(as.character( ss[,3]))== 1)
class(mf)

fit.mf<- survfit(mf ~ gc$cluster)
sd=survdiff(mf ~ gc$cluster)
p.val <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
p.val

plot(fit.mf,col=c("red","green"),xlab="Days after surgery (days)",ylab="Survival")
legend(3500,0.6,col=c("red","green"), c("Basal","Classical"),lty=1,cex=0.5)   
legend("topright",paste("P=",round(p.val, digit=2), sep="") )  

########## Plotting PCSI subtyping heatmap using different classifiers ########################################################################4
library(ComplexHeatmap)
library(circlize)

xx=read.table("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Subtypes/ALL_SUBTYPES.txt",header=T,sep="\t")

y1=sapply(xx[1 ,2:ncol(xx)], function(xx) as.character(xx))
y2=sapply(xx[2 ,2:ncol(xx)], function(xx) as.character(xx))
y3=sapply(xx[3 ,2:ncol(xx)], function(xx) as.character(xx))
y4=sapply(xx[4 ,2:ncol(xx)], function(xx) as.character(xx))

z1 = set.color(y1, type="discrete", color=palette())
z2 = set.color(y2, type="discrete", color=palette())
z3 = set.color(y3, type="discrete", color=palette())
z4 = set.color(y3, type="discrete", color=palette())

ha = HeatmapAnnotation(df = data.frame(Bailey = y1, Moffitt= y2, PAM50=y3, Collisson=y4),
                       col = list(Bailey  = c("Pancreati_Progenitor" = "blue", "Immunogenic" = "green", "Squamous"="red", "ADEX"="yellow"),
                                  Moffitt= c("Basal" = "red", "Classical" = "blue"),
                                  PAM50= c("Basal" = "red", "Classical" = "blue"),
                                  Collisson= c("QM-PDA" = "red", "Classical" = "blue","Exocrine"="yellow")
                       ))



zero_row_mat = matrix(nrow = 0, ncol = 113)
colnames(zero_row_mat) = rep("",113)
ht = Heatmap(zero_row_mat, top_annotation = ha, column_title = "PDAC classifiers")
draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))
decorate_annotation("Bailey", {grid.text("Bailey", unit(-2, "mm"), just = "right")})
decorate_annotation("Moffitt", {grid.text("Moffitt", unit(-2, "mm"), just = "right")})
decorate_annotation("PAM50", {grid.text("PAM50", unit(-2, "mm"), just = "right")})
decorate_annotation("Collisson", {grid.text("Collisson", unit(-2, "mm"), just = "right")})









