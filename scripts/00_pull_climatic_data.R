################################################################################
#                                                                              #
#   PigmentCoordination_MedShrubland                                           #
#   ---------------------------------------------------                        #
#                                                                              #
#   Author:      [Paschalis Chatzopoulos]                                      #
#   Contact:     [paschatzop@gmail.com]                                        #
#   Date:        November 2025                                                 #
#   Paper:       [Pigment disorganization in Salvia rosmarinus under drought   #
#                 and nitrogen addition in a Mediterranean shrubland]          #
#   Data DOI:    [https://doi.org/10.6084/m9.figshare.30636824.v1]             #
#   Manuscript DOI: [10.XXXX/XXXXX]                                            #
#                                                                              #
#   DESCRIPTION:                                                               #
#   This script downloads historical climatic data from the AEMET Open Data    #
#   service for the meteorological station located closest to the study site.  #
#   After identifying suitable stations near the target coordinates, it        #
#   retrieves daily climate records for Ocaña (station ID: 3099Y) across       #
#   multiple time intervals to ensure full temporal coverage (2009–2024). The  #
#   downloaded datasets are merged and exported as a single .csv file for use  #
#   in the associated manuscript analyses.                                     #
#                                                                              #
################################################################################

# Load required packages
library(tidyverse)
library(climaemet)

# Identify meteorological stations near the study coordinates (40° N, -3.6° W)
# The filter selects stations within ±0.5 degrees of the target location.
aemet_stations() %>%
  filter(
    abs(latitud - 40.000000) < 0.5,
    abs(longitud - -3.600000) < 0.5
  )

# The Aranjuez station ("3100B") contains substantial missing data,
# therefore it is excluded in favour of the Ocaña station.
# aranj <- "3100B"

# Selected station ID (Ocaña), used for all subsequent data downloads
ocana <- "3099Y"

## Obtain the API key required for accessing AEMET data
# The following link redirects to the AEMET API key request page:
# browseURL("https://opendata.aemet.es/centrodedescargas/obtencionAPIKey")

## Register the API key
#aemet_api_key("# SET API KEY HERE", install = TRUE)

# The earliest available daily climate data begin on 2009-01-01. 
# To retrieve the complete series up to December 2024, 
# data are downloaded in three consecutive ranges (function limit).

data_daily_0 <-
  aemet_daily_clim(ocana, start = "2009-01-01", end = "2009-12-31")

data_daily_1 <-
  aemet_daily_clim(ocana, start = "2010-01-01", end = "2016-12-31")

data_daily_2 <-
  aemet_daily_clim(ocana, start = "2017-01-01", end = "2024-12-31")

# Combine all data ranges into a single dataset:
data <- bind_rows(data_daily_0, data_daily_1, data_daily_2)

# Export merged dataset to .csv for subsequent climatic and ecological analyses
write.csv(data, "data/climatic.data.csv", row.names = FALSE)

# Export:
write.csv(data, "data/climatic.data.csv", row.names = FALSE)