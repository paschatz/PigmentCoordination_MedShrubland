# Pigment disorganization in *Salvia rosmarinus* under drought and nitrogen addition in a Mediterranean shrubland

This repository contains data and R code accompanying the manuscript:

> Ochoa-Hueso R., Chatzopoulos P., Casimiro-Soriguer R. **Pigment disorganization in *Salvia rosmarinus* under drought and nitrogen addition in a Mediterranean shrubland**. *Functional Plant Biology*, in revision.

We quantify how an extended summer drought and experimental nitrogen (N) addition affect photosynthetic and photoprotective pigments in rosemary (*Salvia rosmarinus*) in a semi-arid Mediterranean shrubland in central Spain. Beyond changes in individual pigment concentrations, we use an ecological coupling framework to evaluate how drought induces pigment disorganization (decoupling) and how recovery after autumn rainfall restores pigment coupling.

## Repository contents
```text
├── data/
│   ├── climatic.data.csv                   # Climatic data pulled from AEMET (scripts/00_pull_climatic_data.R)
│   └── network_data.rdata                  # Data used for network visualisation (scripts/03_pigment_networks_figure3.R)
├── scripts/
│   ├── 00_pull_climatic_data.R             # Pull climatic data directly from AEMET
│   ├── 01_figure1.R                        # Analyse climatic data and export Figure 1
│   ├── 02_analysis.R                       # Main analysis, figures, and tables
│   └── 03_pigment_networks_figure3.R       # Analyse networks and construct Figure 3
├── exports/
│   ├── tables/                             # Tables accompanying the publication
│   └── figures/                            # Figures accompanying the publication
├── LICENSE
├── PigmentCoordination_MedShrubland.Rproj  # R project (run scripts within the R project environment)
└── README.md
```
## Project Overview
This workflow analyzes pigment composition and nitrogen recovery of *Salvia rosmarinus*. The code performs the following:

- Data Retrieval: Automatically fetches raw data from Figshare (API) to ensure immutability and reproducibility.

- Processing: Calculates pigment ratios (e.g., De-epoxidation index) and pigment coupling.

- Statistics: Fits Linear Mixed-Effects Models (nlme) to test for Nitrogen × Treatment interactions and calculates coupling.

- Visualization: Generates all the figures that accompany the publication.

## Usage
Clone or Download this repository.

Open the .Rproj file (e.g., PigmentCoordination_MedShrubland.Rproj) in RStudio. This sets the working directory correctly.

Open the main script: 02_analysis.R 

Run the script from top to bottom.
