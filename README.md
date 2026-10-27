# BCBM RF-SOM Radiomic Phenotyping

R code for the manuscript **"Unsupervised Radiomic Phenotyping of Breast Cancer Brain Metastases via Hybrid Machine Learning"**.

The pipeline ranks radiomic features with a Random Forest, trains a 5 × 5 hexagonal Self-Organizing Map (SOM) on the 15 highest-ranked features, and groups the SOM codebook vectors into four radiomic phenotypes (RP-1 to RP-4) with Ward's D2 clustering. It then tests associations with molecular subtype, receptor status and magnetic field strength, and runs the robustness analyses reported in the paper.

## Data

The data are **not** included in this repository. They are publicly available from The Cancer Imaging Archive (TCIA):

- Collection: <https://www.cancerimagingarchive.net/collection/BCBM-RadioGenomics/>
- Dataset DOI: <https://doi.org/10.7937/RRSE-W278>
- Data descriptor: Taha et al., *Scientific Data* 2025, <https://doi.org/10.1038/s41597-025-06066-6>

Download these two files and place them in the same folder as the script:

- `BCBM-RadioGenomic_Radiomics_Data_v01_20250110.xlsx`
- `BCBM-RadioGenomics-Clinical-data_v01_20250321v2.xlsx`

The data are released under CC BY 4.0 and are subject to the TCIA Data Usage Policy.

## Requirements

- R 4.4.2 (other recent versions should work)
- Packages (installed automatically by the script if missing): `readxl`, `dplyr`, `tidyr`, `stringr`, `tibble`, `ggplot2`, `caret`, `randomForest`, `kohonen`, `corrplot`, `viridis`, `ggdendro`, `gridExtra`, `scales`
- Optional: `cluster` (used in the k-sensitivity analysis)

## How to run

1. Put `BCBM_RF_SOM_analysis.R` and the two `.xlsx` files in one folder.
2. Set that folder as the R working directory.
3. Run:

```r
source("BCBM_RF_SOM_analysis.R")
```

All outputs are written to a new folder, `BCBM_SOM_Outputs/` (figures as 600 dpi PNG files and tables as CSV files).

## What the script does

| Section | Content |
|---|---|
| 0 | Environment, packages, file checks |
| 1 | Data ingestion, cohort flow, quality control |
| 2 | Random Forest feature ranking (1,000 trees, mean decrease in accuracy) |
| 3 | SOM training (5 × 5 hexagonal, 500 epochs) |
| 4 | Figures: correlation heatmap, SOM maps, U-matrix, dendrogram |
| 5 | Chi-square tests with Benjamini–Hochberg correction and Cramér's V |
| 6 | Feature profiles of the four phenotypes |
| 7 | ComBat harmonization of field strength and repeated SOM |
| 8 | Label-free feature-selection controls (Laplacian score, maximum variance) |
| 9 | Patient-level sensitivity analysis |
| 10 | Additional analyses: patient contribution, stratified tests, patient-grouped cross-validation, RF and SOM stability, k sensitivity, alternative pipelines, lesion-size analysis |

Every block in Section 10 runs inside `safe_run()`, so one failing block does not stop the others.

## Reproducibility

- Random seeds: 2025 (main analysis), 2026 (ComBat SOM), 2028 (Laplacian control), 2029 (maximum-variance control). The additional analyses use the seeds listed in the Methods of the paper.
- Main settings: `ntree = 1000`, `mtry = floor(sqrt(p))`, SOM learning rate 0.05 → 0.01, neighbourhood radius 5 → 1, `k = 4`.
- Results can differ slightly across R and package versions.

## Repository contents

```
BCBM_RF_SOM_analysis.R   Main analysis script
results/                 Result tables (CSV) produced by the script
figures/                 Figures (PNG) produced by the script
README.md                This file
LICENSE                  MIT license (code)
```

## Citation

If you use this code, please cite the manuscript (full reference to be added after publication) and the dataset:

Taha, B., Wu, D., Sabal, L., Kollitz, M., Venteicher, A., & Watanabe, Y. (2025). *MRI Dataset of Metastatic Breast Cancer to the Brain with Expert-reviewed Segmentations and Tumor-derived Radiomic Features (BCBM-RadioGenomics)* (Version 1) [Dataset]. The Cancer Imaging Archive. https://doi.org/10.7937/RRSE-W278

## License

The code is released under the MIT License. The data are not part of this repository and follow the TCIA terms (CC BY 4.0).
