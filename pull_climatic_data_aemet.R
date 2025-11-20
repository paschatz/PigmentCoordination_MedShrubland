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
#   In this script I obtain data from AEMET.                                   #
#                                                                              #
#                                                                              #
################################################################################

# Load libraries:
library(tidyverse)
library(climaemet)

# Find nearest station:
aemet_stations() %>%
  # find stations close to 40.000000, -3.600000
  filter(
    abs(latitud - 40.000000) < 0.5,
    abs(longitud - -3.600000) < 0.5
  )

# Aranjuez has many missing values:
#aranj <- "3100B"

# Proceed with Ocaña station:
ocana <- "3099Y"

## Get api key from AEMET
#browseURL("https://opendata.aemet.es/centrodedescargas/obtencionAPIKey")

## Use this function to register your API Key temporarly or permanently
aemet_api_key("# SET API KEY HERE", install = TRUE)

# As of 20-Nov-2025, the earliest available daily data is from 2009-01-01.
# Therefore, we will pull data in three chunks to cover up to 2024-12
data_daily_0 <-
  aemet_daily_clim(ocana, start = "2009-01-01", end = "2009-12-31")

data_daily_1 <-
  aemet_daily_clim(ocana, start = "2010-01-01", end = "2016-12-31")

data_daily_2 <-
  aemet_daily_clim(ocana, start = "2017-01-01", end = "2024-12-31")

# Merge datasets
data <- bind_rows(data_daily_0, data_daily_1, data_daily_2)

# Export:
write.csv(data, "data/climatic.data.csv", row.names = FALSE)