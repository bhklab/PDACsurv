# PDACsurv
Meta-analysis of prognostic models for PDAC

Introduction

Total of 118 pancreatic cancer samples from Pancreatic cancer sequencing initiative (PCSI) cohort was used for training the TSP model to predict patients with early death (within >1 yr) after surgery. Further, to validate the model we used genomic profiles of 610 samples curated from two sequencing cohorts and four array based cohort from public domain. The simplistic top scoring pair approach was used to build the model, where we looked at relative expression of gene pairs within the patient. 

Scripts

1.Consensus_TSP_model.R: The script builds 1000 TSP model using random bootstrapping of PCSI samples and using subset of samples for building the kTSP model. The excluded samples from traning cohort was used for testing where they were dichotomised into two classes with OS > 1 year or less. The models with balanced accuracy more than 0.6 were selected else rejected. 

2.Meta_estimates_calculations_AUC_Dindex_Concordance.R: This script is used to calculate the meta-estimates of the prognostic model using 7 independent cohorts. The AUC, D-index and concordance index was calculated for each cohort and the meta-estimate was calculated and plotted using this script.

3.Clinical_model_building_comparison.R : This script is used to build the linear regression model using 1. clinical features, 2. clinical and OS-TSP models probabilities and 3.OS-TSP model's probabilities alone. The AUC for all the validation cohorts were calculated and the meta-estimate of AUC was calculated for all three models.

4.Clinical_Models_compariosn_Dindex_Cindex.R : The concordance index and C-index was calaculated for all the three models (1.linear regression model using clinical features, 2. linear regression model using clinical and OS-TSP models probabilities, 3. linear regression model using OS-TSP model's probabilities alone). The three models were compared, meta-estimate was calculated and plotted.

5.Subtyping_PDACs.R : The script is used to subtype PCSI samples based on Moffitt, Bailey, PAM50 and Collisson classification system of PDACs and to plot the heatmaps and survival curves based on each subtyping method. A small snipplet is added to compare the classification of samples based on all four classification methods. 

6.ICGC_array_vs_sequencing_comparison.R : The 82 samples common between ICGC-array and ICGC-sequencing were used for prediction and comparing the prediction probabilities using this script. The AUC curves and scatter plot of predicted probabilities were plotted.
