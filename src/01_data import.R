## ________________________________________________________________________

## Title:   Import data
## Purpose: Import and process occurrence data for freshwater fish
## Author:  Dominic Henry & Mo Kajee
## Date:    02/12/2020

## Libraries
library(tidyverse)
library(readxl)
library(sf)
## ________________________________________________________________________

## renv::init() - run this once project code is ready for first commit
## renv::snapshot() - run when new packages are installed (to update renv.lock file)


# Import occurrence data ---------------------------------------------------
# occ_data_all <- map_dfr(dir("data input/species data/", full.names = TRUE),
#                         read_csv)

occ_data_all <- read_csv("data input/Occurrence data_Threatened_14 Sep 2021.csv") %>%
  janitor::clean_names() %>%
  mutate(occurrence_id = str_sub(uuid, 1, 8)) %>%
  select(occurrence_id, taxon, latitude, longitude, sampling_date, conservation_status) %>%
  rename(scientific_name = taxon,
         decimal_latitude = latitude,
         decimal_longitude = longitude)


## Check number of records
occ_data_all %>%
  group_by(scientific_name) %>%
  tally() %>%
  print(n = 25)

# Filter records ----------------------------------------------------------

## Add quality filters (obvs will need to keep more of the original cols to do this)


# Create spatial data -----------------------------------------------------
latlongCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
za <- st_read("data input/RSA_fixed.shp")

## Create sf object (first filter non-valid lat/long values)
occ_data_sf <- occ_data_all %>%
  filter(!(is.na(decimal_latitude) | is.na(decimal_longitude))) %>%
  filter(!(decimal_latitude == "NULL" | decimal_longitude == "NULL")) %>%
  mutate_at(c("decimal_latitude", "decimal_longitude"), as.numeric) %>%
  st_as_sf(coords = c("decimal_longitude","decimal_latitude"), crs = latlongCRS)

## Plot check
ggplot()+
  geom_sf(data = occ_data_sf)+
  geom_sf(data = za, fill = NA, color = "black")

# Write workspaces --------------------------------------------------------
save(list = c("occ_data_all","occ_data_sf"), file = "data output/occ_data_clean.RData")

