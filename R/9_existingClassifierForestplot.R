source("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Github/Scripts/predict_meta_estimates_function.R")

#################################################################################################################
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/pcosp_Scores.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/birnbaum_Scores.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/haider_Scores.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_Scores_sig.score.RData")
load("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/survival_comparison_published_model.RData")

################################################################################################################
################################################################################################################

pcosp= predict_meta_estimates( pcosp_prob$pcsi, pcosp_prob$tcga, pcosp_prob$kirby, pcosp_prob$icgc_arr, pcosp_prob$unc, pcosp_prob$winter, pcosp_prob$collisson, pcosp_prob$ouh,
                               survival$pcsi_OS, survival$pcsi_OS_Status,
                               survival$tcga_OS, survival$tcga_OS_Status,
                               survival$kirby_OS, survival$kirby_OS_Status,
                               survival$icgc_array_OS, survival$icgc_array_OS_Status,
                               survival$unc_OS, survival$unc_OS_Status,
                               survival$winter_OS, survival$winter_OS_Status,
                               survival$collisson_OS, survival$collisson_OS_Status,
                               survival$ouh_OS, survival$ouh_OS_Status
)

birnbaum= predict_meta_estimates( birnbaum_prob$pcsi, birnbaum_prob$tcga, birnbaum_prob$kirby, birnbaum_prob$icgc_arr, birnbaum_prob$unc, birnbaum_prob$winter, birnbaum_prob$collisson, birnbaum_prob$ouh,
                                  survival$pcsi_OS, survival$pcsi_OS_Status,
                                  survival$tcga_OS, survival$tcga_OS_Status,
                                  survival$kirby_OS, survival$kirby_OS_Status,
                                  survival$icgc_array_OS, survival$icgc_array_OS_Status,
                                  survival$unc_OS, survival$unc_OS_Status,
                                  survival$winter_OS, survival$winter_OS_Status,
                                  survival$collisson_OS, survival$collisson_OS_Status,
                                  survival$ouh_OS, survival$ouh_OS_Status
)

haider=  predict_meta_estimates( haider_prob$pcsi, haider_prob$tcga, haider_prob$kirby, haider_prob$icgc_arr, haider_prob$unc, haider_prob$winter, haider_prob$collisson, haider_prob$ouh,
                                 survival$pcsi_OS, survival$pcsi_OS_Status,
                                 survival$tcga_OS, survival$tcga_OS_Status,
                                 survival$kirby_OS, survival$kirby_OS_Status,
                                 survival$icgc_array_OS, survival$icgc_array_OS_Status,
                                 survival$unc_OS, survival$unc_OS_Status,
                                 survival$winter_OS, survival$winter_OS_Status,
                                 survival$collisson_OS, survival$collisson_OS_Status,
                                 survival$ouh_OS, survival$ouh_OS_Status
)

chen= predict_meta_estimates( chen_prob$pcsi, chen_prob$tcga, chen_prob$kirby, chen_prob$icgc_arr, chen_prob$unc, chen_prob$winter, chen_prob$collisson, chen_prob$ouh,
                              survival$pcsi_OS, survival$pcsi_OS_Status,
                              survival$tcga_OS, survival$tcga_OS_Status,
                              survival$kirby_OS, survival$kirby_OS_Status,
                              survival$icgc_array_OS, survival$icgc_array_OS_Status,
                              survival$unc_OS, survival$unc_OS_Status,
                              survival$winter_OS, survival$winter_OS_Status,
                              survival$collisson_OS, survival$collisson_OS_Status,
                              survival$ouh_OS, survival$ouh_OS_Status
)



##### Plotting Forestplot of  D index ###############
############################################### ##################################################################
############################################### ##################################################################
############################################### ##################################################################
############################################### ##################################################################

r.mean <- c( NA,log2(haider$dindex_micro$estimate), log2(chen$dindex_micro$estimate),  log2(birnbaum$dindex_micro$estimate), log2(pcosp$dindex_micro$estimate)  , NA,
             NA, log2(haider$dindex_seq$estimate), log2(chen$dindex_seq$estimate), log2(birnbaum$dindex_seq$estimate),log2(pcosp$dindex_seq$estimate),    NA,
             NA,log2(haider$dindex_meta$estimate), log2(chen$dindex_meta$estimate), log2(birnbaum$dindex_meta$estimate), log2(pcosp$dindex_meta$estimate))

r.lower <- c( NA,log2(haider$dindex_micro_lower),  log2(chen$dindex_micro_lower),   log2(birnbaum$dindex_micro_lower),log2(pcosp$dindex_micro_lower), NA,
              NA,log2(haider$dindex_seq_lower),log2(chen$dindex_seq_lower),  log2(birnbaum$dindex_seq_lower), log2(pcosp$dindex_seq_lower), NA, 
              NA,log2(haider$dindex_meta_lower),log2(chen$dindex_meta_lower),  log2(birnbaum$dindex_meta_lower),log2(pcosp$dindex_meta_lower) )

r.upper <-c( NA,log2(haider$dindex_micro_upper), log2(chen$dindex_micro_upper),log2(birnbaum$dindex_micro_upper),  log2(pcosp$dindex_micro_upper), NA,
             NA,log2(haider$dindex_seq_upper),log2(chen$dindex_seq_upper),   log2(birnbaum$dindex_seq_upper), log2(pcosp$dindex_seq_upper), NA,
             NA, log2(haider$dindex_meta_upper), log2(chen$dindex_meta_upper), log2(birnbaum$dindex_meta_upper),  log2(pcosp$dindex_meta_upper))


r.pval <- round(c(NA,haider$dindex_micro_pval,chen$dindex_micro_pval,  birnbaum$dindex_micro_pval, pcosp$dindex_micro_pval, NA,
                  NA, haider$dindex_seq_pval, chen$dindex_seq_pval,birnbaum$dindex_seq_pval,pcosp$dindex_seq_pval,  NA, 
                  NA,haider$dindex_meta_pval,  chen$dindex_meta_pval, birnbaum$dindex_meta_pval,pcosp$dindex_meta_pval),2)

r.pval1 <- c(NA, sprintf("%.1E", haider$dindex_micro_pval), sprintf("%.1E", chen$dindex_micro_pval), sprintf("%.1E", birnbaum$dindex_micro_pval),sprintf("%.1E", pcosp$dindex_micro_pval), NA,
             NA, sprintf("%.1E", haider$dindex_seq_pval), sprintf("%.1E", chen$dindex_seq_pval), sprintf("%.1E", birnbaum$dindex_seq_pval),sprintf("%.1E", pcosp$dindex_seq_pval),  NA,
             NA, sprintf("%.1E", haider$dindex_meta_pval), sprintf("%.1E", chen$dindex_meta_pval), sprintf("%.1E", birnbaum$dindex_meta_pval),sprintf("%.1E", pcosp$dindex_meta_pval)
)

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("Microarry Cohorts","Haider, 2014", "Chen, 2015","Birnbaum, 2017","PCOSP",NA, 
                  "Sequencing cohorts","Haider, 2014", "Chen, 2015","Birnbaum, 2017","PCOSP",NA, 
                  "Overall","Haider, 2014", "Chen, 2015","Birnbaum, 2017","PCOSP"
)

data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -18L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Classifiers",rownames(t)),
  c("P values",r.pval1))
#pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/New_Figure_DINDEX.pdf")
fn1 <- local({
  i = 0
  
  s_clrs =c(c("#666666","#666666","#666666","#E7298A"),c("#666666","#666666","#666666","#E7298A"),c("#666666","#666666","#666666","#E7298A") )
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})
forestplot(tabletext2,data2,xlab="Log2 D-index", new_page=FALSE, is.summary=c( rep(TRUE,18)),
           clip=c(-1,4),txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  xlab  = gpar(fontfamily = "Helvetica", cex = 1)),
           col = fpColors(text="black"),title=" ",zero=0,graphwidth=unit(2, "inches"),  align=c("l"),  fn.ci_sum = fn1,boxsize = 0.25)
#dev.off()
################################################################################################################
################################################################################################################
################################################################################################################

##### Plotting Forestplot of Concordance index ###############



r.mean <- c( NA, (haider$con_micro$estimate),  (chen$con_micro$estimate),(birnbaum$con_micro$estimate), (pcosp$con_micro$estimate),NA,
             NA, (haider$con_seq$estimate),  (chen$con_seq$estimate), (birnbaum$con_seq$estimate), (pcosp$con_seq$estimate),  NA,
             NA,(haider$con_meta$estimate) ,  (chen$con_meta$estimate),(birnbaum$con_meta$estimate),  (pcosp$con_meta$estimate)
             
)

r.lower <- c( NA, (haider$con_micro_lower),  (chen$con_micro_lower),(birnbaum$con_micro_lower), (pcosp$con_micro_lower),NA,
              NA, (haider$con_seq_lower),  (chen$con_seq_lower), (birnbaum$con_seq_lower), (pcosp$con_seq_lower),  NA,
              NA,(haider$con_meta_lower) ,  (chen$con_meta_lower),(birnbaum$con_meta_lower),  (pcosp$con_meta_lower)
              
)
r.upper <- c( NA, (haider$con_micro_upper),  (chen$con_micro_upper),(birnbaum$con_micro_upper), (pcosp$con_micro_upper),NA,
              NA, (haider$con_seq_upper),  (chen$con_seq_upper), (birnbaum$con_seq_upper), (pcosp$con_seq_upper),  NA,
              NA,(haider$con_meta_upper) ,  (chen$con_meta_upper),(birnbaum$con_meta_upper),  (pcosp$con_meta_upper)
              
)

r.pval <- round(c(NA,haider$con_micro_pval,chen$con_micro_pval,  birnbaum$con_micro_pval, pcosp$con_micro_pval, NA,
                  NA, haider$con_seq_pval, chen$con_seq_pval,birnbaum$con_seq_pval,pcosp$con_seq_pval,  NA, 
                  NA,haider$con_meta_pval,  chen$con_meta_pval, birnbaum$con_meta_pval,pcosp$con_meta_pval),2)


r.pval1 <- c(NA, sprintf("%.1E", haider$con_micro_pval), sprintf("%.1E", chen$con_micro_pval), sprintf("%.1E", birnbaum$con_micro_pval),sprintf("%.1E", pcosp$con_micro_pval), NA,
             NA, sprintf("%.1E", haider$con_seq_pval), sprintf("%.1E", chen$con_seq_pval), sprintf("%.1E", birnbaum$con_seq_pval),sprintf("%.1E", pcosp$con_seq_pval),  NA,
             NA, sprintf("%.1E", haider$con_meta_pval), sprintf("%.1E", chen$con_meta_pval), sprintf("%.1E", birnbaum$con_meta_pval),sprintf("%.1E", pcosp$con_meta_pval)
)

t <- cbind(r.mean ,r.lower,r.upper,r.pval)
rownames(t) <-  c("Microarry Cohorts","Haider, 2014", "Chen, 2015","Birnbaum, 2017","PCOSP",NA, 
                  "Sequencing cohorts","Haider, 2014", "Chen, 2015","Birnbaum, 2017","PCOSP",NA, 
                  "Overall","Haider, 2014", "Chen, 2015","Birnbaum, 2017","PCOSP"
)

data2 <- 
  structure(list(
    mean  = c(NA,t[,1]),
    lower = c(NA,t[,2]),
    upper = c(NA,t[,3])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -18L), 
    class = "data.frame")


tabletext2<-cbind(
  c("Classifiers",rownames(t)),
  c("P values",r.pval1))
#pdf("/Users/vandanasandhu/Desktop/Project1-Metadatasubtyping/Figures/New_Figure_CONCORDANCE.pdf")

fn1 <- local({
  i = 0
  
  s_clrs =c(c("#666666","#666666","#666666","#E7298A"),c("#666666","#666666","#666666","#E7298A"),c("#666666","#666666","#666666","#E7298A") )
  function(..., col){
    i <<- i + 1
    fpDrawSummaryCI(...,col=s_clrs[i])
  }
})

forestplot(tabletext2,data2,xlab="Concordance Index",new_page=FALSE, is.summary=c( rep(TRUE,18)),
           fn.ci_sum = fn1, clip=c(0.4,0.9), txt_gp = fpTxtGp(label = gpar(fontfamily = "Helvetica"),ticks = gpar(cex=0.8),  
          xlab  = gpar(fontfamily = "Helvetica", cex = 1)), col = fpColors(text="black"),title="",zero=0.5,graphwidth=unit(2, "inches"),align=c("l"), boxsize = 0.25)
                                                                 xlab  = gpar(fontfamily = "Helvetica", cex = 1)), col = fpColors(text="black"),title="",zero=0.5,graphwidth=unit(2, "inches"),align=c("l"), boxsize = 0.25)

############################################### ##################################################################
############################################### ##################################################################
############################################### ##################################################################
############################################### ##################################################################

######################################################################################################################


pcosp_prediction = c(pcosp_prob$pcsi, pcosp_prob$tcga, pcosp_prob$kirby, pcosp_prob$icgc_arr, pcosp_prob$unc, pcosp_prob$winter, pcosp_prob$collisson, pcosp_prob$ouh)
birnbaum_prediction = c(birnbaum_prob$pcsi, birnbaum_prob$tcga, birnbaum_prob$kirby, birnbaum_prob$icgc_arr, birnbaum_prob$unc, birnbaum_prob$winter, birnbaum_prob$collisson, birnbaum_prob$ouh)
chen_prediction = c(chen_prob$pcsi, chen_prob$tcga, chen_prob$kirby, chen_prob$icgc_arr, chen_prob$unc, chen_prob$winter, chen_prob$collisson, chen_prob$ouh)
haider_prediction = c(haider_prob$pcsi, haider_prob$tcga, haider_prob$kirby, haider_prob$icgc_arr, haider_prob$unc, haider_prob$winter, haider_prob$collisson, haider_prob$ouh)

cor.test(pcosp_prediction, birnbaum_prediction, method="spearman")
cor.test(pcosp_prediction, chen_prediction, method="spearman")
cor.test(pcosp_prediction, haider_prediction, method="spearman")

par(mfrow=c(3,1))
plot(pcosp_prediction, birnbaum_prediction, col="green", pch=15)
plot(pcosp_prediction, chen_prediction, col="red", pch=15)
plot(pcosp_prediction, haider_prediction, col="blue", pch=15)

zz=cbind(PCOSP= pcosp_prediction, BIRNBUAM= birnbaum_prediction, CHEN= chen_prediction, HAIDER=haider_prediction)
pairs(zz, upper.panel = panel.cor)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method="spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


#### Meta concordance index comparison
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

pcosp_birnbaum_cindex = cindex.comp.meta(list.cindex1 = list("pcsi"=pcosp$con_pcsi, "tcga"= pcosp$con_tcga, "kirby"= pcosp$con_kirby, "icgc_array" =pcosp$con_icgc_array, "unc"= pcosp$con_unc, "ouh"= pcosp$con_ouh, "winter" = pcosp$con_winter, "collisson"= pcosp$con_collisson),
                 list.cindex2 = list("pcsi"=birnbaum$con_pcsi, "tcga"= birnbaum$con_tcga, "kirby"= birnbaum$con_kirby, "icgc_array" =birnbaum$con_icgc_array, "unc"= birnbaum$con_unc, "ouh"= birnbaum$con_ouh, "winter" = birnbaum$con_winter, "collisson"= birnbaum$con_collisson))
pcosp_birnbaum_cindex_pval=pcosp_birnbaum_cindex$p.value/2


pcosp_haider_cindex =cindex.comp.meta(list.cindex1 = list("pcsi"=pcosp$con_pcsi, "tcga"= pcosp$con_tcga, "kirby"= pcosp$con_kirby, "icgc_array" =pcosp$con_icgc_array, "unc"= pcosp$con_unc, "ouh"= pcosp$con_ouh, "winter" = pcosp$con_winter, "collisson"= pcosp$con_collisson),
                 list.cindex2 = list("pcsi"=haider$con_pcsi, "tcga"= haider$con_tcga, "kirby"= haider$con_kirby, "icgc_array" =haider$con_icgc_array, "unc"= haider$con_unc, "ouh"= haider$con_ouh, "winter" = haider$con_winter, "collisson"= haider$con_collisson))
pcosp_haider_cindex_pval = pcosp_haider_cindex$p.value/2

pcosp_chen_cindex= cindex.comp.meta(list.cindex1 = list("pcsi"=pcosp$con_pcsi, "tcga"= pcosp$con_tcga, "kirby"= pcosp$con_kirby, "icgc_array" =pcosp$con_icgc_array, "unc"= pcosp$con_unc, "ouh"= pcosp$con_ouh, "winter" = pcosp$con_winter, "collisson"= pcosp$con_collisson),
                 list.cindex2 = list("pcsi"=chen$con_pcsi, "tcga"= chen$con_tcga, "kirby"= chen$con_kirby, "icgc_array" =chen$con_icgc_array, "unc"= chen$con_unc, "ouh"= chen$con_ouh, "winter" = chen$con_winter, "collisson"= chen$con_collisson))
pcosp_chen_cindex_pval=pcosp_chen_cindex$p.value/2


#########################################################################################################################
#########################################################################################################################

pcosp_birnbaum_dindex = dindex.comp.meta(list.dindex1 = list("pcsi"=pcosp$dindex_pcsi, "tcga"= pcosp$dindex_tcga, "kirby"= pcosp$dindex_kirby, "icgc_array" =pcosp$dindex_icgc_array, "unc"= pcosp$dindex_unc, "ouh"= pcosp$dindex_ouh, "winter" = pcosp$dindex_winter, "collisson"= pcosp$dindex_collisson),
                 list.dindex2 = list("pcsi"=birnbaum$dindex_pcsi, "tcga"= birnbaum$dindex_tcga, "kirby"= birnbaum$dindex_kirby, "icgc_array" =birnbaum$dindex_icgc_array, "unc"= birnbaum$dindex_unc, "ouh"= birnbaum$dindex_ouh, "winter" = birnbaum$dindex_winter, "collisson"= birnbaum$dindex_collisson))
pcosp_birnbaum_dindex_pval=pcosp_birnbaum_dindex$p.value/2


pcosp_haider_dindex = dindex.comp.meta(list.dindex1 = list("pcsi"=pcosp$dindex_pcsi, "tcga"= pcosp$dindex_tcga, "kirby"= pcosp$dindex_kirby, "icgc_array" =pcosp$dindex_icgc_array, "unc"= pcosp$dindex_unc, "ouh"= pcosp$dindex_ouh, "winter" = pcosp$dindex_winter, "collisson"= pcosp$dindex_collisson),
                 list.dindex2 = list("pcsi"=haider$dindex_pcsi, "tcga"= haider$dindex_tcga, "kirby"= haider$dindex_kirby, "icgc_array" =haider$dindex_icgc_array, "unc"= haider$dindex_unc, "ouh"= haider$dindex_ouh, "winter" = haider$dindex_winter, "collisson"= haider$dindex_collisson))
pcosp_haider_dindex_pval = pcosp_haider_dindex$p.value/2


pcosp_chen_dindex= dindex.comp.meta(list.dindex1 = list("pcsi"=pcosp$dindex_pcsi, "tcga"= pcosp$dindex_tcga, "kirby"= pcosp$dindex_kirby, "icgc_array" =pcosp$dindex_icgc_array, "unc"= pcosp$dindex_unc, "ouh"= pcosp$dindex_ouh, "winter" = pcosp$dindex_winter, "collisson"= pcosp$dindex_collisson),
                 list.dindex2 = list("pcsi"=chen$dindex_pcsi, "tcga"= chen$dindex_tcga, "kirby"= chen$dindex_kirby, "icgc_array" =chen$dindex_icgc_array, "unc"= chen$dindex_unc, "ouh"= chen$dindex_ouh, "winter" = chen$dindex_winter, "collisson"= chen$dindex_collisson))
pcosp_chen_dindex_pval = pcosp_chen_dindex$p.value/2


