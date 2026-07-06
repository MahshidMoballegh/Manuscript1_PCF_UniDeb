[Statistic_Olink_IO_Neu_OD.md](https://github.com/user-attachments/files/29704408/Statistic_Olink_IO_Neu_OD.md)

Statistical Analysis Pipeline

This file contains explaination of  the  statistical workflow supporting our manuscript. The pipeline is designed for high-throughput proteomic data analysis (Olink NPX values) using automated univariate screening, false discovery rate (FDR) control, and advanced visualization engines.

 📂 Codebase Reference
Core Script: `Statistic_Olink_IO_Neu_OD_2.R` implements the entire exploratory and differential abundant pipeline.



 🛠️ Data Pipeline & Architecture

The workflow ingests the proteomic data matrix and branches into multiple systematic sub-group comparisons:

                      [ Final_data_3.xlsx ]
                                │
                                ▼
                  【 Dynamic Patient Screening 】
           Filters out samples with entirely missing assays
                                │
     ┌──────────────────────────┼──────────────────────────┐
     ▼                                                     ▼
【 Cohort Comparisons 】                        【 Longitudinal Trajectories 】
• Overall VR vs VD                               Tracks specific targets
Disease Type, Gender, and  
CNS Levels 1, 2, and 3    




🔬 Statistical Methodology & Parameters

Pre-processing & Quality Control
Missing Data Handling:Samples with entirely missing data across all proteomic assays are programmatically excluded before comparison.
Log2 Fold Change (FC): Since Olink NPX values are delivered on a log_2scale, the Normalized Fold Change is computed directly as the difference between group medians: Delta NPX = median(X) - median(Y).

Hypothesis Testing & Multiple Testing Correction
Univariate Testing: Group comparisons are conducted using the non-parametric Wilcoxon Rank-Sum Test , ensuring robustness against non-normally distributed proteomic expression values.
FDR Management: Multiple testing correction is handled via Storey's `qvalue` package to maintain high statistical power. The script automatically falls back to standard Benjamini-Hochberg (FDR) adjustment if the q-value optimization target fails to converge due to highly flat or restricted p-value distributions.
Significance Thresholds: Features are isolated using strict custom logical limits:
  - `FC_THRESHOLD <- 0.5` (|Log2FC| > 0.5)
  - `Q_THRESHOLD  <- 0.05` (q-value < 0.05)

Automated Visualizations
Volcano Plots: Generated via `ggplot2` and annotated using `ggrepel` to isolate and highlight statistically dynamic Increased and Decreased  proteins in VR compare to VD.
Significance Boxplots:Automatically processes all filtered significant proteins. Advanced boxplots are rendered with built-in significance stars mapped directly from the computed q-values:
  - `***` (q value < 0.001)
  - `**` ( q value < 0.01)
  - `*` (qvalue < 0.05)
  - `ns` (Not Significant)
Expression Heatmaps:Implemented using both `ComplexHeatmap` and `pheatmap` libraries to display row-scaled, group-annotated expression profiles across clinical subsets.



📊 Comprehensive Matrix of Computed Comparisons

The pipeline executes a multi-loop evaluation across these key clinical vectors:

| Comparison Stream | Group 1 vs Group 2 | Analytical Context |
| :--- | :--- | :--- |
| **Visits soecifics** | `VR` vs `VD` | Overall response trajectory |
| **Disease Specific** | `T-ALL` vs `B-ALL` | Run separately inside `VD` and `VR` windows and response trajectory |
| **Sex Specific** | `Female (F)` vs `Male (M)` | Run separately inside `VD` and `VR` windows and response trajectory |
| **CNS Stratification** | `CNS 1` vs `2` \| `CNS 1` vs `3` \| `CNS 2` vs `3` | Run separately inside `VD` and `VR` windows and response trajectory |
| **Combined Targets** | `CNS Group 1` vs `CNS Group 2+3` | Consolidated CNS involvement inside `VD` and `VR` |

---

## 📦 System Dependencies

To initialize the statistical environment, execute the following script in your R console:


# Install CRAN packages
install.packages(c("tidyverse", "readxl", "readr", "dplyr", "tidyr", "ggplot2", "ggrepel", "tibble", "purrr", "pheatmap"))

# Install Bioconductor framework dependencies
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("qvalue", "ComplexHeatmap", "grid"))

🚀 Execution Guide
Place your data file named exactly Final_data_3.xlsx inside your active R working directory.

Source and run Statistic_Olink_IO_Neu_OD_2.R.

Review Outputs: The script automatically generates your results directory structures:

3/Tables/ — Contains detailed .csv files for every single comparison layer.

3/Volcano_Plots/ — Individual volcano plots mapping out specific cohort shifts.

3/Boxplots/ & 3/Heatmaps/ — Abundance of profiles for downstream candidate proteins.

CNS1 vs 2+3/ & CNS1 vs 2+3 (VR)/ — Target sub-directories for consolidated CNS streams.
