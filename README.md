# PDACsurv
Meta-analysis of prognostic models for PDAC

Introduction

A unique set of 89 pancreatic cancer samples profiled using both sequencing and microarray platform was used for training the PCOSP (Pancreatic cancer overall survival predictor) to predict patients with early death (within >1 yr) after surgery. Further, to validate the model we used genomic profiles of 823 samples curated from public domain. The simplistic k-top scoring pair approach was used to build the model, where we looked at relative expression of gene pairs within the patient. 

Scripts

1. PCOSP_model_building.R : The script is used for building Pancreatic cancer overall survival predictor using unique 89 samples profiled using microarray and sequencing platform.

2. Meta_estimate_Calculation_Validation_Cohorts.R : The script calculates D-index and Concordance index for independent cohorts and calculate the meta-estimae of C-index and D-index. The forestplot of C-index and D-index was plotted for all the cohorts.

3. PCOSP_score_estimation.R: It is used for calculating the PCOSP score for independent datasets.

4. Validation cohorts formatting.R: Used for formatting independent datasets for input to PCOSP.



