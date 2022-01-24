## ________________________________________________________________________

## Title:    EST Medium sensitivity
## Purpose:  Process SDMs for inclusion into Medium level of EST
## Author:   Dominic Henry
## Date:     21/01/2022

## Libraries
library(raster)
library(tidyverse)
library(sf)
library(glue)
library(readxl)
## ________________________________________________________________________

latlongCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
aeaproj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

spp_clip_ref <- read_xlsx("data input/clipping_parameters_EST.xlsx") %>%
  filter(iucn_status != "CR")

# Clip using IUCN range ---------------------------------------------------
iucn_spp <- spp_clip_ref %>%
  filter(clipping_extent == "iucn range") %>%
  pull(scientific_name)

binpaths <- glue("data output/sdm BART results/{iucn_spp}/SDM_bin_mean{iucn_spp}.tif")

bins <- map(binpaths,
            raster)

ras_to_poly <- function(spp, spp_bin){

  iucn_poly <- st_read(glue("data output/sdm data processing/{spp}/IUCN_range_{spp}.shp")) %>%
    st_transform(crs = aeaproj)

  range <- raster::mask(spp_bin, iucn_poly)
  range <- reclassify(range, c(-Inf,0.999999,NA,1,1,1))
  range <- rasterToPolygons(range, dissolve=TRUE)
  pp <- st_as_sf(range)
  pp_cast <- sf::st_cast(pp, "POLYGON")

  pp_cast <- pp_cast %>%
    select(-1) %>%
    mutate(OBJECTID = str_c("a", 1:nrow(pp_cast)),
           SENSFEAT = spp,
           SENSITIVIT = "Medium",
           THEME = "Animals",
           DEVELOP = NA)
  polyarea <- as.numeric(st_area(pp_cast))   # Remove small segments
  keep <- which(polyarea > 800000)
  pp_cast <- pp_cast[keep,]

  pp_cast <- pp_cast %>%
    select(OBJECTID,SENSITIVIT,THEME,SENSFEAT,DEVELOP) %>%
    st_transform(crs = latlongCRS)

  write_sf(pp_cast, glue("data output/screening tool layers/medium/{spp}.shp"),
           delete_dsn = TRUE)

}

## Write individual file
map2(iucn_spp, bins,
      .f = safely(ras_to_poly))


# Clip using Primary Catchment --------------------------------------------

catchments <- st_read("data output/primary_catchment_latlong.shp", crs = latlongCRS) %>%
  st_transform(crs = aeaproj)

prim_spp <- spp_clip_ref %>%
  filter(clipping_extent == "primary catchment") %>%
  pull(scientific_name)

binpaths <- glue("data output/sdm BART results/{prim_spp}/SDM_bin_mean{prim_spp}.tif")

bins <- map(binpaths,
            raster)

ras_to_poly <- function(spp, spp_bin){

  prim_ids <- spp_clip_ref[which(spp_clip_ref$scientific_name == spp), "catchment_id"]

  regions <- str_c("Region ", str_split(prim_ids, pattern = ",")[[1]]) %>%
    str_replace(., "  ", " ")

  prim_poly <- catchments %>%
    filter(name %in% regions)

  range <- raster::mask(spp_bin, prim_poly)
  range <- reclassify(range, c(-Inf,0.999999,NA,1,1,1))
  range <- rasterToPolygons(range, dissolve=TRUE)
  pp <- st_as_sf(range)
  pp_cast <- sf::st_cast(pp, "POLYGON")

  pp_cast <- pp_cast %>%
    select(-1) %>%
    mutate(OBJECTID = str_c("a", 1:nrow(pp_cast)),
           SENSFEAT = spp,
           SENSITIVIT = "Medium",
           THEME = "Animals",
           DEVELOP = NA)
  polyarea <- as.numeric(st_area(pp_cast))   # Remove small segments
  keep <- which(polyarea > 800000)
  pp_cast <- pp_cast[keep,]

  pp_cast <- pp_cast %>%
    select(OBJECTID,SENSITIVIT,THEME,SENSFEAT,DEVELOP) %>%
    st_transform(crs = latlongCRS)

  write_sf(pp_cast, glue("data output/screening tool layers/medium/{spp}.shp"),
           delete_dsn = TRUE)

}

## Write individual file
map2(prim_spp, bins,
     .f = safely(ras_to_poly))


# Clip using Secondary Catchment ------------------------------------------
secondary <- st_read("data input/Catchment Layers/secondary_catchment_area/secondary_catchment_area.shp") %>%
  st_transform(crs = aeaproj) %>%
  st_zm(drop = TRUE)

second_spp <- spp_clip_ref %>%
  filter(clipping_extent == "secondary catchment") %>%
  pull(scientific_name)

binpaths <- glue("data output/sdm BART results/{second_spp}/SDM_bin_mean{second_spp}.tif")

bins <- map(binpaths,
            raster)

# spp <- second_spp[1]
# spp_bin <- bins[[1]]

ras_to_poly <- function(spp, spp_bin){

  second_ids <- spp_clip_ref[which(spp_clip_ref$scientific_name == spp), "catchment_id"]

  regions <- str_split(second_ids, pattern = ",")[[1]] %>%
    str_replace(., "  ", " ") %>%
    str_trim()

  second_poly <- secondary %>%
    filter(name %in% regions)

  range <- raster::mask(spp_bin, second_poly)
  range <- reclassify(range, c(-Inf,0.999999,NA,1,1,1))
  range <- rasterToPolygons(range, dissolve=TRUE)
  pp <- st_as_sf(range)
  pp_cast <- sf::st_cast(pp, "POLYGON")

  pp_cast <- pp_cast %>%
    select(-1) %>%
    mutate(OBJECTID = str_c("a", 1:nrow(pp_cast)),
           SENSFEAT = spp,
           SENSITIVIT = "Medium",
           THEME = "Animals",
           DEVELOP = NA)
  polyarea <- as.numeric(st_area(pp_cast))   # Remove small segments
  keep <- which(polyarea > 800000)
  pp_cast <- pp_cast[keep,]

  pp_cast <- pp_cast %>%
    select(OBJECTID,SENSITIVIT,THEME,SENSFEAT,DEVELOP) %>%
    st_transform(crs = latlongCRS)

  write_sf(pp_cast, glue("data output/screening tool layers/medium/{spp}.shp"),
           delete_dsn = TRUE)

}

## Write individual file
map2(second_spp, bins,
     .f = safely(ras_to_poly))


# Clip using Quartenery ---------------------------------------------------
quartenery <- st_read("data input/Catchment Layers/quaternary_catchment/quaternary_catchment.shp") %>%
  st_transform(crs = aeaproj) %>%
  st_zm(drop = TRUE)

quart_spp <- spp_clip_ref %>%
  filter(clipping_extent == "quartenery") %>%
  pull(scientific_name)

binpaths <- glue("data output/sdm BART results/{quart_spp}/SDM_bin_mean{quart_spp}.tif")

bins <- map(binpaths,
            raster)

ras_to_poly <- function(spp, spp_bin){

  quart_ids <- spp_clip_ref[which(spp_clip_ref$scientific_name == spp), "catchment_id"]

  regions <- str_split(quart_ids, pattern = ",")[[1]] %>%
    str_replace(., "  ", " ") %>%
    str_trim()

  quart_poly <- quartenery %>%
    filter(name %in% regions)

  range <- raster::mask(spp_bin, quart_poly)
  range <- reclassify(range, c(-Inf,0.999999,NA,1,1,1))
  range <- rasterToPolygons(range, dissolve=TRUE)
  pp <- st_as_sf(range)
  pp_cast <- sf::st_cast(pp, "POLYGON")

  pp_cast <- pp_cast %>%
    select(-1) %>%
    mutate(OBJECTID = str_c("a", 1:nrow(pp_cast)),
           SENSFEAT = spp,
           SENSITIVIT = "Medium",
           THEME = "Animals",
           DEVELOP = NA)
  polyarea <- as.numeric(st_area(pp_cast))   # Remove small segments
  keep <- which(polyarea > 800000)
  pp_cast <- pp_cast[keep,]

  pp_cast <- pp_cast %>%
    select(OBJECTID,SENSITIVIT,THEME,SENSFEAT,DEVELOP) %>%
    st_transform(crs = latlongCRS)

  write_sf(pp_cast, glue("data output/screening tool layers/medium/{spp}.shp"),
           delete_dsn = TRUE)

}

## Write individual file
map2(quart_spp, bins,
     .f = safely(ras_to_poly))



# Clip using IUCN + Primary -----------------------------------------------
com_spp <- spp_clip_ref %>%
  filter(clipping_extent == "iucn_and_primary") %>%
  pull(scientific_name)

binpaths <- glue("data output/sdm BART results/{com_spp}/SDM_bin_mean{com_spp}.tif")

bins <- map(binpaths,
            raster)

# spp <- com_spp[1]
# spp_bin <- bins[[1]]

ras_to_poly <- function(spp, spp_bin){

  prim_ids <- spp_clip_ref[which(spp_clip_ref$scientific_name == spp), "catchment_id"]

  regions <- str_c("Region ", str_split(prim_ids, pattern = ",")[[1]]) %>%
    str_replace(., "  ", " ")

  prim_poly <- catchments %>%
    filter(name %in% regions)

  iucn_poly <- st_read(glue("data output/sdm data processing/{spp}/IUCN_range_{spp}.shp")) %>%
    st_transform(crs = aeaproj)

  poly <- st_intersection(prim_poly, iucn_poly)

  range <- raster::mask(spp_bin, poly)
  range <- reclassify(range, c(-Inf,0.999999,NA,1,1,1))
  range <- rasterToPolygons(range, dissolve=TRUE)
  pp <- st_as_sf(range)
  pp_cast <- sf::st_cast(pp, "POLYGON")

  pp_cast <- pp_cast %>%
    select(-1) %>%
    mutate(OBJECTID = str_c("a", 1:nrow(pp_cast)),
           SENSFEAT = spp,
           SENSITIVIT = "Medium",
           THEME = "Animals",
           DEVELOP = NA)
  polyarea <- as.numeric(st_area(pp_cast))   # Remove small segments
  keep <- which(polyarea > 800000)
  pp_cast <- pp_cast[keep,]

  pp_cast <- pp_cast %>%
    select(OBJECTID,SENSITIVIT,THEME,SENSFEAT,DEVELOP) %>%
    st_transform(crs = latlongCRS)

  write_sf(pp_cast, glue("data output/screening tool layers/medium/{spp}.shp"),
           delete_dsn = TRUE)

}

## Write individual file
map2(com_spp, bins,
     .f = safely(ras_to_poly))


# Combine all medium shapefiles -------------------------------------------

shpfiles <- dir("data output/screening tool layers/medium", full.names = TRUE, pattern =".shp")
shpfiles

map(shpfiles,
    read_sf) %>%
  do.call(rbind, .) %>%
  write_sf("data output/screening tool layers/FWF_medium.shp", delete_layer = TRUE,
           overwrite = TRUE)

