# PDACsurv
Meta-analysis of prognostic models for PDAC

Introduction

Total of 118 pancreatic cancer samples from Pancreatic cancer sequencing initiative (PCSI) cohort was used for training the TSP model to predict patients with early death (within >1 yr) after surgery. Further, to validate the model we used genomic profiles of 610 samples curated from two sequencing cohorts and four array based cohort from public domain. The simplistic top scoring pair approach was used to build the model, where we looked at relative expression of gene pairs within the patient. 

Scripts

1.Consensus_TSP_model.R: The script builds 1000 TSP model using random bootstrapping of PCSI samples and using subset of samples for building the kTSP model. The excluded samples from traning cohort was used for testing where they were dichotomised into two classes with OS > 1 year or less. The models with balanced accuracy more than 0.6 were selected else rejected. 

2.Meta_estimates_calculations_AUC_Dindex_Concordance.R: This script is used to calculate the meta-estimates of the prognostic model using 7 independent cohorts. The AUC, D-index and concordance index was calculated for each cohort and the meta-estimate was calculated and plotted using this script.
