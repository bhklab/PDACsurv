load("PDAC_Expression_dataset.RData")
########################## ###################### ###################### ###################### ###################### ###################### 
#### Validation cohorts 
###################### ICGC

icgc_cohort=rs_coh$ICGC_seq
icgc_mat<-data.matrix(sapply(icgc_cohort[1:nrow(icgc_cohort) ,1:(ncol(icgc_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(icgc_mat)=rownames(icgc_cohort)

g1=which(as.numeric(as.character(icgc_cohort$OS))<=365 &  as.numeric(as.character(icgc_cohort$OS_Status))==1);g2=which(as.numeric(as.character(icgc_cohort$OS))>365)
g_icgc_seq=sort(c(g1,g2))
icgc_seq_grp=ifelse(as.numeric(as.character(icgc_cohort$OS))>=365,1,0)[g_icgc_seq]

###################### TCGA
tcga_cohort=rs_coh$TCGA
tcga_mat<-data.matrix(sapply(tcga_cohort[1:nrow(tcga_cohort) ,1:(ncol(tcga_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(tcga_mat)=rownames(tcga_cohort) 
g1=which(as.numeric(as.character(tcga_cohort$OS))<=365 &  as.numeric(as.character(tcga_cohort$OS_Status))==1);g2=which(as.numeric(as.character(tcga_cohort$OS))>365)
g_tcga=sort(c(g1,g2))
tcga_grp=ifelse(as.numeric(as.character(tcga_cohort$OS))>=365,1,0)[g_tcga]

###################### UNC
unc_cohort=rs_coh$UNC
unc_mat<-data.matrix(sapply(unc_cohort[1:nrow(unc_cohort) ,1:(ncol(unc_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(unc_mat)=rownames(unc_cohort)
g1=which(as.numeric(as.character(unc_cohort$OS))<=365 &  as.numeric(as.character(unc_cohort$OS_Status))==1);g2=which(as.numeric(as.character(unc_cohort$OS))>365)
g_unc=sort(c(g1,g2))
unc_grp=ifelse(as.numeric(as.character(unc_cohort$OS))>=365,1,0)[g_unc]

###################### OUH
ouh_cohort=rs_coh$OUH
ouh_mat<-data.matrix(sapply(ouh_cohort[1:nrow(ouh_cohort), 1:(ncol(ouh_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(ouh_mat)=rownames(ouh_cohort)
g1=which(as.numeric(as.character(ouh_cohort$OS))<=365 &  as.numeric(as.character(ouh_cohort$OS_Status))==1);g2=which(as.numeric(as.character(ouh_cohort$OS))>365)
g_ouh=sort(c(g1,g2))
ouh_grp=ifelse(as.numeric(as.character(ouh_cohort$OS))>=365,1,0)[g_ouh]

###################### zhang
zhang_cohort=rs_coh$Zhang
zhang_mat<-data.matrix(sapply(zhang_cohort[1:nrow(zhang_cohort) ,1:(ncol(zhang_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(zhang_mat)=rownames(zhang_cohort)
g1=which(as.numeric(as.character(zhang_cohort$OS))<=365 &  as.numeric(as.character(zhang_cohort$OS_Status))==1);g2=which(as.numeric(as.character(zhang_cohort$OS))>365)
g_zhang=sort(c(g1,g2))
zhang_grp=ifelse(as.numeric(as.character(zhang_cohort$OS))>=365,1,0)[g_zhang]

###################### Winter
winter_cohort=rs_coh$Winter
winter_mat<-data.matrix(sapply(winter_cohort[1:nrow(winter_cohort) ,1:(ncol(winter_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(winter_mat)=rownames(winter_cohort) 
g1=which(as.numeric(as.character(winter_cohort$OS))<=365 &  as.numeric(as.character(winter_cohort$OS_Status))==1);g2=which(as.numeric(as.character(winter_cohort$OS))>365)
g_winter=sort(c(g1,g2))
winter_grp=ifelse(as.numeric(as.character(winter_cohort$OS))>=365,1,0)[g_winter]

###################### ICGC_array
icgc_array_cohort=rs_coh$ICGC_arr
icgc_array_mat<-data.matrix(sapply(icgc_array_cohort[1:nrow(icgc_array_cohort) ,1:(ncol(icgc_array_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(icgc_array_mat)=rownames(icgc_array_cohort)
g1=which(as.numeric(as.character(icgc_array_cohort$OS))<=365 &  as.numeric(as.character(icgc_array_cohort$OS_Status))==1);g2=which(as.numeric(as.character(icgc_array_cohort$OS))>365)
g_icgc_arr=sort(c(g1,g2))
icgc_array_grp=ifelse(as.numeric(as.character(icgc_array_cohort$OS))>=365,1,0)[g_icgc_arr]

###################### Collisson
collisson_cohort=rs_coh$Collisson
collisson_mat<-data.matrix(sapply(collisson_cohort[1:nrow(collisson_cohort) ,1:(ncol(collisson_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(collisson_mat)=rownames(collisson_cohort)
g1=which(as.numeric(as.character(collisson_cohort$OS))<=365 &  as.numeric(as.character(collisson_cohort$OS_Status))==1);g2=which(as.numeric(as.character(collisson_cohort$OS))>365)
g_coll=sort(c(g1,g2))
collisson_grp=ifelse(as.numeric(as.character(collisson_cohort$OS))>=365,1,0)[g_coll]

###################### Chen
chen_cohort=rs_coh$Chen
chen_mat<-data.matrix(sapply(matrix(chen_cohort[1:nrow(chen_cohort) ,1:(ncol(chen_cohort)-2)]), function(xx) as.numeric(as.character(xx))))
rownames(chen_mat)=chen_cohort[,1]
colnames(chen_mat)=colnames(chen_cohort)[1:(ncol(chen_cohort)-2)]
g1=which(as.numeric(as.character(chen_cohort$OS))<=365 &  as.numeric(as.character(chen_cohort$OS_Status))==1); g2=which(as.numeric(as.character(chen_cohort$OS))>365)
g_chen=sort(c(g1,g2))
chen_grp=ifelse(as.numeric(as.character(chen_cohort$OS))>=365,1,0)[g_chen]

###################### Kirby
kirby_cohort=rs_coh$Kirby
kirby_mat<-data.matrix(sapply(matrix(kirby_cohort[1:nrow(kirby_cohort) ,1:(ncol(kirby_cohort)-2)]), function(xx) as.numeric(as.character(xx))))
rownames(kirby_mat)=kirby_cohort[,1]
colnames(kirby_mat)=colnames(kirby_cohort)[1:(ncol(kirby_cohort)-2)]
g1=which(as.numeric(as.character(kirby_cohort$OS))<=365 &  as.numeric(as.character(kirby_cohort$OS_Status))==1); g2=which(as.numeric(as.character(kirby_cohort$OS))>365)
g_kirby=sort(c(g1,g2))
kirby_grp=ifelse(as.numeric(as.character(kirby_cohort$OS))>=365,1,0)[g_kirby]

###################### PCSI
pcsi_cohort=rs_coh$PCSI
pcsi_mat<-data.matrix(sapply(pcsi_cohort[1:nrow(pcsi_cohort) ,1:(ncol(pcsi_cohort)-2)], function(xx) as.numeric(as.character(xx))))
rownames(pcsi_mat)=rownames(pcsi_cohort)
g1=which(as.numeric(as.character(pcsi_cohort$OS))<=365 &  as.numeric(as.character(pcsi_cohort$OS_Status))==1);g2=which(as.numeric(as.character(pcsi_cohort$OS))>365)
g_pcsi=sort(c(g1,g2))
pcsi_grp=ifelse(as.numeric(as.character(pcsi_cohort$OS))>=365,1,0)[g_pcsi]


###################################################################################################################################################################
###################################################################################################################################################################
