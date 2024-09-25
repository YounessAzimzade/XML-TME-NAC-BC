This repository accompanies the paper titled *"Explainable Machine Learning Unveils the Role of the Breast Tumor Microenvironment in Neoadjuvant Chemotherapy Outcome."* It contains all the code and processed data used in the paper. The steps to recreate the results are as follows:

1) **Estimating Cell Fractions:**  
   Given the bulk gene expression profile, the first step is to estimate the cell fractions in your samples. To do this, use the Signature Matrices available in the corresponding folder. CIBERSORTx is recommended as a tool for this purpose. After estimating the cell fractions and retrieving the clinical outcome data for your samples, you will have a dataset similar to `AllSamples.csv`.

2) **Exploring the Role of Cell Types in Clinical Outcomes:**  
   Once you have decided how to divide your data into discovery and validation cohorts, explore the role of cell type frequency in the clinical outcome of your discovery cohort. Use the series of codes in `Discovery_SVM.ipynb`. These scripts take cell fractions, clinical features, and clinical outcomes, then train a model to predict the outcome. They employ SHAP to extract feature importance, returning SHAP values in a dataset similar to `AllSamples.csv`. You will also need to run a similar process on your validation cohort, as outlined in `Validation_SVM.ipynb`.

3) **Calculating pCR Scores for the General Population:**  
   After calculating the SHAP values, run the `Pipeline.R` script. This script takes the SHAP values and cell fractions and calculates the "pCR Scores" for the general population.

4) **pCR Scores for Subtypes:**  
   You can repeat steps 2 and 3 for your subtype of interest to calculate the pCR Scores for that specific subtype.
