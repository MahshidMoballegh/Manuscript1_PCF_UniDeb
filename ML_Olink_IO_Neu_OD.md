[ML_Olink_IO_Neu_OD.md](https://github.com/user-attachments/files/29704430/ML_Olink_IO_Neu_OD.md)
Proteomics Machine Learning Pipeline 
This file contains explaination of R pipeline for feature selection, biomarker discovery, and predictive modeling using high-throughput proteomics NPX data.
The project employs three distinct machine learning algorithms (Random Forest, LASSO, and Support Vector Machines) across two clinical classification workflows.
Pipeline Architecture
The script processes data sequentially through two primary analytical modules:  

                      Final_data_3.xlsx
                              │
                     ┌────────┴────────┐
                     ▼                 ▼
           Module 1: Visits    Module 2: CNS Target (VD Only and VR only) 
                 VD vs VR        CNS_1 vs CNS_2+3
                     │                 │
     ┌───────────────┼──────────────┐  └───────────────┬──────────────┐
     ▼               ▼              ▼                  ▼              ▼              ▼
[Random Forest]   [LASSO]       [SVM-RFE]       [Random Forest]    [LASSO]       [SVM-RFE]
     │               │              │                  │              │              │
     └───────────────┬──────────────┘                  └──────────────┬──────────────┘
                     ▼                                                ▼
              [Venn Diagram]                                   [Venn Diagram]
                     │                                                │
       [Bootstrap Screen & Subsets]                              [Bootstrap Screen & Subsets]


🛠️ Installation & Dependencies
Ensure you have the required R packages installed before executing the script. You can initialize them by running:
Rinstall.packages(c("tidyverse", "readxl", "randomForest", "pROC", "glmnet", "e1071", "caret", "VennDiagram", "writexl"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "circlize"))


📂 Project Structure & OutputsThe script dynamically builds structured output directories separating analytical results:
├── Final_data_3.xlsx                 # Input Dataset (Required)
├── MachineLearning/                  # MODULE 1 OUTPUTS (Visits: VD vs VR)
│   ├── RF/                           # Importance metrics, ROC plots, CV Confusion Matrices
│   ├── LASSO/                        # CV Curves, coefficient weights, boxplots
│   ├── SVM/                          # Tuning parameter grids, RFE lists, heatmaps
│   ├── Phase1_Screening_Results.xlsx # Core + 1 Bootstrap evaluations
│   └── Bootstrap_final_result.xlsx   # Combined optimized subset performance metrics
│
└── MachineLearning_CNS/              # MODULE 2 OUTPUTS (CNS Target in VD samples)
    ├── RF/                           # Gini charts, 5-fold CV diagnostics, Target Boxplots
    ├── LASSO/                        # Active feature coefficients and verification matrices 
    ├── SVM/                          # Support Vector RFE rankings and grouped expressions
    └── Venn_Proteins_Summary.txt     # Complete overlap profile intersection report


🔬Core Workflow Details
Data Cleaning & Transformation Handles clinical metadata factoring for variables like DiseaseType (T-ALL = 1, B-ALL = 2) and Gender (F = 1, M = 2).
Eliminates empty strings or string-formatted "NA" records to ensure computational mathematical parsing.
Scales and centers assay features systematically when managing regression configurations (LASSO/SVM).

Modeling & Feature Selection
Random Forest: Extracts Mean Decrease Gini ranks across standard splits and tracks stratification stability across a 5-fold cross-validation engine.
LASSO: Employs penalized regularization (L_1se) utilizing cv.glmnet under the conservative lambda.1se (1-standard-error) threshold to extract sparse, high-confidence biomarkers.
SVM-RFE: Couples a Radial Basis Function (RBF) kernel with Recursive Feature Elimination (caret::rfe) utilizing cross-validation loops to isolate optimal subset widths.

Advanced Combinatorial Screening (Module 1 Only)
Integrates a systematic 1,000-iteration bootstrap screening utility (run_svm_bootstrap) to measure Out-Of-Bag (OOB) AUC performance.
It screens structural combinations of common intersections (Core Group) extended by additional candidates (Winners) to determine the absolute highest-performing predictive protein signature.


📊 Visualizations
GeneratedThe pipeline yields publication-ready diagnostic figures inside their respective subfolders:
Confusion Matrices: Clean, grid-mapped tile charts utilizing customizable fill layers (plot_cm).
Dot plots: Top feature metrics grouped cleanly inside structured performance brackets (SizeBin).
Boxplots: Faceted NPX distribution charts tracking individual biomarker expression changes across cohorts.
Complex Heatmaps: Deep sample column slices sorted strictly by target groups alongside annotated scaling arrays.
Venn Diagrams: Triple-set overlap charts identifying algorithm-agnostic biomarker signatures.

🚀 Getting Started
Place your data file named Final_data_3.xlsx into your working directory.Source the R file or execute it line-by-line in your environment:
Rsource("ml_proteomics_pipeline.R")
Monitor progress readouts printed to the R console. 
Once finalized, check the newly generated directories (MachineLearning/ and MachineLearning_CNS) for all tables, logs, and graphics.
