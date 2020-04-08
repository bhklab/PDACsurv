load("PDAC_Expression_dataset.RData")
tcga_cohort=rs_coh$TCGA

######### Birnbaum
library(genefu)
gene_coeff=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/bmc_genes.txt", sep="\t", header =T )

### Extracting 36 Birnbaum classifier genes
#### Calculating signature score using sig.score for each cohort
data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/pcsi_dataset.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
pcsi_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/tcga_dataset.txt",sep="\t", header=T)
data=data[which(tcga_cohort$OS > 180), ]
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
tcga_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/kirby_dataset.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
kirby_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/icgc_arr_dataset.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
icgc_arr_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/unc_dataset.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
unc_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/ouh_dataset.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
ouh_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/winter_dataset.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
winter_prob=ss$score

data=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/bmc_med/collisson.txt",sep="\t", header=T)
ss= sig.score(x=gene_coeff,data= data,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = FALSE, verbose = FALSE)
collisson_prob=ss$score

birnbaum_prob= list(pcsi=pcsi_prob, tcga =tcga_prob, icgc_arr =icgc_arr_prob, kirby= kirby_prob, unc= unc_prob, ouh=ouh_prob, winter= winter_prob, collisson= collisson_prob)
save(birnbaum_prob, file="birnbaum_Scores.RData")

########################################################################################################################################################
############## Chen

### Extracting 15 Chen classifier genes
#### Calculating signature score using sig.score for each cohort

gene_coeff1=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/plos_2015/chen_genes.txt", sep="\t", header =T )

pp_pcsi=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/pcsi.txt", header = T, sep="\t")
ss_pcsi= sig.score(x=gene_coeff1,data=pp_pcsi,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score

pp_tcga=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/tcga.txt", header = T, sep="\t")
ss_tcga= sig.score(x=gene_coeff1,data=pp_tcga,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score
ss_tcga=ss_tcga[which( tcga_cohort$OS > 180)]

pp_kirby=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/kirby.txt", header = T, sep="\t")
ss_kirby= sig.score(x=gene_coeff1,data=pp_kirby,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score

pp_unc=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/unc.txt", header = T, sep="\t")
ss_unc= sig.score(x=gene_coeff1,data=pp_unc,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed =TRUE,verbose = TRUE)$score

pp_ouh=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/ouh.txt", header = T, sep="\t")
ss_ouh= sig.score(x=gene_coeff1,data=pp_ouh,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score

pp_winter=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/winter.txt", header = T, sep="\t")
ss_winter= sig.score(x=gene_coeff1,data=pp_winter,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score

pp_collisson=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/collisson.txt", header = T, sep="\t")
ss_collisson= sig.score(x=gene_coeff1,data=pp_collisson,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score

pp_icgc_arr= read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/chen_sig.score/icgc_arr.txt", header = T, sep="\t")
ss_icgc_arr= sig.score(x=gene_coeff1,data=pp_icgc_arr,  do.mapping = FALSE, mapping, size = 0,cutoff = NA, signed = TRUE, verbose = TRUE)$score

chen_prob= list(pcsi=ss_pcsi, tcga =ss_tcga, icgc_arr =ss_icgc_arr, kirby= ss_kirby, unc= ss_unc, ouh=ss_ouh, winter= ss_winter, collisson= ss_collisson)
save(chen_prob, file="chen_Scores_sig.score.RData")

########################################################################################################################################################
#### PCOSP

pcosp_prob= list(pcsi=unlist(pcsi_list), tcga =unlist(tcga_list), icgc_arr = unlist(icgc_array_list), kirby= unlist(kirby_list), unc= unlist(unc_list), ouh=unlist(ouh_list), winter=unlist(winter_list), collisson= unlist(collisson_list))
save(pcosp_prob, file="pcosp_Scores.RData")


########################################################################################################################################################
############# Haider


### Scores are provided by the author 


pcsi=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/PCSI/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header = T)
tcga=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/TCGA/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header=T)
icgc_array=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/IA/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header=T)
kirby=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/Kirby/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header = T)
unc=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/UNC/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header = T)
ouh=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/OUH/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header = T)
winter=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/Winter/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header = T)
collisson=read.table("/Users/vandanasandhu/Desktop/RS_Remodelling_TSP_project/Data_from_Haider/Collisson/mRNA_abundance_geneNames_TRAINING_mRNA_abundance_TN_DE_limmaFDR_0.01_FC_1_uv_cox_P_0.05_patient_posterior_probability.txt", sep=" ", header = T)


pcsi_prob=pcsi$X1
tcga_prob=tcga$X1
tcga_OS = tcga_cohort$OS[which(tcga_cohort$OS > 180) ]
tcga_prob=tcga_prob[which(tcga_cohort$OS > 180) ]


icgc_arr_prob=icgc_array$X1
kirby_prob=kirby$X1
unc_prob=unc$X1
ouh_prob=ouh$X1
winter_prob=winter$X1
collisson_prob=collisson$X1

haider_prob= list(pcsi=pcsi_prob, tcga =tcga_prob, icgc_arr =icgc_arr_prob, kirby= kirby_prob, unc= unc_prob, ouh=ouh_prob, winter= winter_prob, collisson= collisson_prob)
save(haider_prob, file="haider_Scores.RData")





