# Reject-to-Remeasure Rep
This repository is associated to the paper "Rejections Based on Predictive Uncertainty Enable Reliable Routine Soil Spectroscopy". It contains all the code and data to reproduce the analysis of the study.

# Structure
0_Descriptive_Analysis -> Contain the R-scripts to plot dataset properties. Contain scripts for Fig. 6b and 6c as well as parts of Fig. 1.

1_Prepare_Analysis -> Contain R and Python-scripts to prepare CV splits and hyperparameter search space.

2_Apply_models -> Contain the Python-scripts to run PLSQR, TabPFN, TabICL and QRF in the outer and inner CV loop. 

3_Analyse -> Contain the R-scripts to analyse the results and reproduce all Figures associated to the analysis.

Data -> Contain folds, hyperparameter search space and the dataset.

Figures -> Contain created Figures by R and Python-scripts.

Results -> Contain the predictions of the different models for the outer folds and the tuning results.

docs -> Contain the html for the R-tutorial.

# Contact
If you have questions or trouble to run the code, please contact me (Jonas Schmidinger) at "jonas.schmidinger@uni-osnabrueck.de". 
