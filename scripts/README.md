# Scripts Overview

This folder contains all R scripts used to retrieve data, process climatic and pigment information, run statistical analyses, and generate the figures included in the publication:

> Ochoa-Hueso R., Chatzopoulos P., Casimiro-Soriguer R. **Pigment disorganization in *Salvia rosmarinus* under drought and nitrogen addition in a Mediterranean shrubland.** *Functional Plant Biology*, in revision.

All paths referenced inside the scripts are relative to the repository root.

------------------------------------------------------------------------

## Script descriptions

### **00_pull_climatic_data.R**

Retrieves climatic data directly from **AEMET** using the `climaemet` package.
The script: - Downloads daily precipitation and temperature data for the Ocaña station.
- Formats dates and variables.
- Exports the final dataset to: `data/climatic.data.csv`

**Note:** To run this script, you must obtain an API key from AEMET and set it in your R environment.

------------------------------------------------------------------------

### **01_figure1.R**

Processes the climatic dataset and generates **Figure 1** of the manuscript.
(Data are loaded separately so users without an AEMET API key can still run the analysis.)

The script: - Cleans and aggregates climatic variables.
- Identifies the extended drought period.
- Visualises precipitation and temperature trends.
- Exports Figure 1.

------------------------------------------------------------------------

### **02_analysis.R**

Main analysis workflow for the study.
The script: - Imports processed pigment data.
- Calculates pigment ratios (e.g., de-epoxidation index).
- Computes pigment coupling using the `eco_coupling()` function.
- Fits linear mixed-effects models (LMMs) to test:
  - Drought vs. recovery effects
  - Nitrogen treatment effects
- Exports tables and figures to: `exports/tables/` and `exports/figures/`

------------------------------------------------------------------------

### **03_pigment_networks_figure3.R**

Produces the pigment network visualisation corresponding to **Figure 3**.
The script: 
- Loads pigment correlation matrices.
- Builds network objects (nodes = pigments, edges = correlations).
- Visualises networks under:
  - Summer drought
  - Recovery
- Exports:
  - Network figure: `exports/figures/figure3_networks.pdf`

------------------------------------------------------------------------

## How to run these scripts

1.  Open the R project: `PigmentCoordination_MedShrubland.Rproj`
2.  Run the scripts in the following order:
    -   `00_pull_climatic_data.R` (if you have an AEMET API key)
    -   `01_figure1.R`
    -   `02_analysis.R`
    -   `03_pigment_networks_figure3.R`

All outputs will be saved automatically to their respective folders.

------------------------------------------------------------------------

## Dependencies

The scripts require the following R packages:

### Install CRAN packages

``` r
install.packages(c(
  "tidyverse", "climaemet", "ggrepel", "nlme", "patchwork", "broom", "Hmisc",
  "purrr", "viridis", "svglite", "httr2", "jsonlite", "ggpubr", "igraph",
  "ggraph", "tidygraph", "cowplot"
))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("paschatz/EcoCoupleR")
```

