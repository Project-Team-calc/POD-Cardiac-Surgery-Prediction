# Prediction of Postoperative Delirium in Adult Cardiac Surgery Patients

This repository contains the dataset and analysis code for the study: **"Prediction of Postoperative Delirium in Adult Cardiac Surgery Patients Without Overt Preoperative Cognitive Impairment: A Prospective, Dual-Center Study and Web-Based Calculator"**

## 📌 Project Overview
The objective of this study was to develop and validate a prediction model for postoperative delirium (POD) specifically in cardiac surgery patients who screen negative for baseline cognitive deficits.

Key methodological features include:
* **Bayesian Logistic Regression:** Applied to analyze predictors with adjustment for clinical centers.
* **Prior Sensitivity Analysis:** Rigorous evaluation of model stability across different prior specifications (Uninformative, Weakly Informative, and Skeptical priors).
* **External Validation:** Performance assessment on an independent clinical center to ensure model generalizability.
* **Landmark Analysis:** Assessment of risk trajectories at the 24-hour postoperative mark.

## 📂 File Descriptions
* **POD_Revision.Rproj**: RStudio Project file for correct environment configuration.
* **POD_Revision.R**: The comprehensive analysis script, including:
  * Data preprocessing and variable rescaling.
  * Bayesian model construction and Prior Sensitivity Analysis.
  * Visualization of model robustness (Forest plots).
  * Internal validation (Bootstrap) and external validation (ROC).
  * Landmark analysis and clinical decision tools (Nomogram and Calibration plots).
* **my_data_cleaned_for_analysis.csv**: The anonymized dataset used for model development and validation.

## 🛠 How to Reproduce
1. Clone this repository or download the ZIP file.
2. Open **`POD_Revision.Rproj`** in RStudio.
3. Ensure the required R packages are installed: `brms`, `pROC`, `rms`, `tidyverse`, `ggplot2`, `ggsci`, `mice`.
4. Run **`POD_Revision.R`** to reproduce the posterior estimates, sensitivity results, and clinical validation metrics.

## 🌐 Interactive Web Calculator
To facilitate clinical translation, the validated model has been operationalized into a web-based risk calculator:  
👉 [https://project-team-calc.shinyapps.io/Risk_Calculator_App/](https://project-team-calc.shinyapps.io/Risk_Calculator_App/)

## 🔒 Data Privacy & Ethics
The dataset (`my_data_cleaned_for_analysis.csv`) has been thoroughly de-identified. All patient-identifiable information (e.g., names, identifiers, specific dates) has been removed in accordance with institutional ethics guidelines to ensure patient privacy.

## ✉️ Contact
For questions regarding the methodology or analysis code, please contact:
**Fang Chen, M.S.**  
Department of Joint Surgery and Sports Medicine, The Second Qilu Hospital of Shandong University  
Email: yeechen80@126.com / medusayang520@gmail.com
