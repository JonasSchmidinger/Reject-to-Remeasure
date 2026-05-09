# Reject-to-Remeasure Rep
This repository is associated to the paper "Reject-to-Remeasure: Uncertainty-Guided Quality Control for Reliable Soil Spectroscopy".

# Structure
0_Descriptive_Analysis -> Contains the R-scripts to plot dataset properties.

1_Prepare_Analysis -> Contains R and Python-scripts to prepare CV splits and hyperparameter search space.

2_Apply_models -> Contains the Python-scripts to run PLSQR, TabPFN, TabICL and QRF in the outer and inner CV loop.

3_Analyse -> Contains the R-scripts to analyse the results.

Data -> Contains folds, hyperparameter search space and in the future the dataset. At the moment "Full_Dataset.csv" contains only the spectral data and shows the dataset structure. However, soil data publication is currently prepapred. 

Figures -> Contains created Figures by R and Python-scripts.

Results -> Contains the predictions of the different models

docs -> Contains the html for the R-tutorial.

# Contact
If you have trouble to run the code, please contact me (Jonas Schmidinger) at "jonas.schmidinger@uni-osnabrueck.de". Data will be accessible in the future but you may contact Marc-Olivier Gasser at "marc-o.gasser@irda.qc.ca" or MAPAQ (Ministère de l’Agriculture, des Pêcheries et de l’Alimentation) for earlier data access under closed agreements.
