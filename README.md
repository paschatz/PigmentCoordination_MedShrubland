# Data and Code for: [Insert Paper Title Here]
This repository contains the R script used to generate the statistical analyses and figures for the manuscript:

Authors: [List Authors] Title: [Insert Paper Title] Journal: [Insert Journal Name] (Year) DOI: [Insert DOI Link]

## Project Overview
This workflow analyzes pigment composition and nitrogen recovery in [Plant Species Name]. The code performs the following:

- Data Retrieval: Automatically fetches raw data from Figshare (API) to ensure immutability and reproducibility.

- Processing: Calculates pigment ratios (e.g., De-epoxidation index).

- Statistics: Fits Linear Mixed-Effects Models (nlme) to test for Nitrogen × Treatment interactions and calculates coupling.

- Visualization: Generates all the figures that accompany the publication.

## Usage
Clone or Download this repository.

Open the .Rproj file (e.g., PigmentCoordination_MedShrubland.Rproj) in RStudio. This sets the working directory correctly.

Open the main script: analysis.R.

Run the script from top to bottom.
