## ________________________________________________________________________

## Title:    EST Very High sensitivity
## Purpose:  Process FBIS occurrence data for development of Very High level of EST
## Author:   Dominic Henry
## Date:     24/01/2022

## Libraries
library(raster)
library(tidyverse)
library(sf)
library(glue)
library(readxl)
library(lubridate)
## ________________________________________________________________________

latlongCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
aeaproj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

spp_list <- read_xlsx("data input/clipping_parameters_EST.xlsx") %>%
  filter(iucn_status == "CR") %>%
  pull(scientific_name)

load("data output/occ_data_clean.RData")
occ_data_sf

fepa_rivers <- st_read("data input/River FEPAs/River_FEPAs.shp")

points_to_poly <- function(spp){

  spp_data <- occ_data_sf %>%
    filter(scientific_name %in% spp)

  poly_id <- st_intersection(spp_data, fepa_rivers) %>%
    pull(UNIT_ID) %>%
    unique(.)

  very_high_polys <- fepa_rivers %>%
    filter(UNIT_ID %in% poly_id)

  pp_cast <- sf::st_cast(very_high_polys, "POLYGON")

  pp_cast <- pp_cast %>%
    select(-1) %>%
    mutate(OBJECTID = str_c("a", 1:nrow(pp_cast)),
           SENSFEAT = spp,
           SENSITIVIT = "High",
           THEME = "Animals",
           DEVELOP = NA)
  polyarea <- as.numeric(st_area(pp_cast))   # Remove small segments
  keep <- which(polyarea > 800000)
  pp_cast <- pp_cast[keep,]

  pp_cast <- pp_cast %>%
    select(OBJECTID,SENSITIVIT,THEME,SENSFEAT,DEVELOP) %>%
    st_transform(crs = latlongCRS)

  write_sf(pp_cast, glue("data output/screening tool layers/very high/{spp}.shp"),
           delete_dsn = TRUE)


}

map(.x = spp_list,
    .f= safely(points_to_poly))

# Combine all high shapefiles --------------------------------------------
shpfiles <- dir("data output/screening tool layers/very high", full.names = TRUE, pattern =".shp")
shpfiles

map(shpfiles,
    read_sf) %>%
  do.call(rbind, .) %>%
  write_sf("data output/screening tool layers/FWF_very_high.shp", delete_layer = TRUE,
           overwrite = TRUE)
