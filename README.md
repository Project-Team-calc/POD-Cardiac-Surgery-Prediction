# POD-Cardiac-Surgery-Prediction
Data and code for the prediction model of POD in cardiac surgery.
# Prediction of Postoperative Delirium in Adult Cardiac Surgery Patients

This repository contains the dataset and analysis code for the study: 
**"Prediction of Postoperative Delirium in Adult Cardiac Surgery Patients Without Overt Preoperative Cognitive Impairment: A Prospective, Dual-Center Study and Web-Based Calculator"**

## 📌 Project Overview
The objective of this study was to develop and validate a prediction model for postoperative delirium (POD) specifically in cardiac surgery patients who screen negative for baseline cognitive deficits. We employed **Bayesian multilevel logistic regression** to account for center-level heterogeneity and performed **bidirectional external validation** to ensure model robustness.

## 📂 File Descriptions
- **`POD_Revision.Rproj`**: RStudio Project file. We recommend opening this first to ensure the working directory is correctly set.
- **`POD_Revision.R`**: The main R script containing:
  - Data preprocessing and variable rescaling.
  - Bayesian multilevel model construction using weakly informative priors.
  - Internal validation (Bootstrap) and Bidirectional external validation.
  - Landmark analysis at the 24-hour postoperative mark.
  - Generation of the Nomogram and Calibration plots.
- **`my_data_cleaned_for_analysis.csv`**: The anonymized dataset used for all analyses.

## 🛠 How to Reproduce
1. Clone this repository or download the ZIP file.
2. Open **`POD_Revision.Rproj`** in RStudio.
3. Open **`POD_Revision.R`** and run the code. 
   - *Note: Ensure you have the necessary R packages installed (e.g., `brms`, `pROC`, `rms`).*

## 🌐 Interactive Web Calculator
To facilitate clinical translation, we have developed a user-friendly, web-based risk calculator accessible at:
👉 **[https://project-team-calc.shinyapps.io/Risk_Calculator_App/](https://project-team-calc.shinyapps.io/Risk_Calculator_App/)**

## 🔒 Data Privacy & Ethics
The provided dataset (`my_data_cleaned_for_analysis.csv`) has been thoroughly de-identified. All patient-identifiable information (e.g., names, IDs, exact dates) has been removed to comply with institutional ethics guidelines and protect patient privacy.

## ✉️ Contact
For questions regarding the methodology or code, please contact:
**Fang Chen, M.S.**  
Corresponding Author  
Department of Joint Surgery and Sports Medicine, The Second Qilu Hospital of Shandong University  
Email: yeechen80@126.com/medusayang520@gmail.com
