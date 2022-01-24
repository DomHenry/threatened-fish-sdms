## ________________________________________________________________________

## Title:   SDM data import
## Purpose: Process raster and occurrence data in prep for fish SDMs
## Author:  Dominic Henry & Mo Kajee
## Date:    08/09/2021

## Libraries
library(lubridate)
library(gridExtra)
library(ggpubr)
library(sf)
library(viridis)
library(readxl)
library(raster)
library(tidyverse)
library(glue)
## ________________________________________________________________________

x <- "data input/Queries/SDM_query_Oreochromis mossambicus.xlsx"

load("data output/occ_data_clean.RData")
start <- Sys.time()

# Import query details ----------------------------------------------------
query <- read_xlsx(x) %>%
  print

sppselect <- query$Value[which(query$Input == "Species")]

occ_data_sf_all <- occ_data_sf %>%
  filter(scientific_name %in% sppselect)

plot(occ_data_sf_all$geometry)

# Temporal filter ---------------------------------------------------------

## Keep all records for now -- No temporal filter ##
occ_data_sf_all <- occ_data_sf_all %>%
  mutate(year = year(sampling_date)) %>%
  filter(year >= as.numeric(query$Value[which(query$Input == "Start year")])) %>%
  filter(year <= as.numeric(query$Value[which(query$Input == "End year")]))

# Create output folder ----------------------------------------------------
if(dir.exists(glue("data output/sdm data processing/{sppselect}"))) {
  print("Folder exists")
} else {
  print("Folder created")
  dir.create(glue("data output/sdm data processing/{sppselect}"))
}

# Spatial filter ----------------------------------------------------------

## Read in IUCN map
tempfile <- tempfile()
unzip(zipfile = glue("data input/Fish Distribution Map_IUCN/{sppselect}.zip"),
      exdir = tempfile)

range_files <- list.files(tempfile, pattern = ".shp$",full.names=TRUE)
range <- st_read(range_files)

## Write range map to file
range %>%
  st_write(glue("data output/sdm data processing/{sppselect}/IUCN_range_{sppselect}.shp"),
           delete_dsn = TRUE)

## Check plots
plot(occ_data_sf_all$geometry)
plot(range$geometry, add = TRUE)
plot(st_buffer(range,0.1)$geometry, add = TRUE)

## KEEP ALL POINTS FOR CERTAIN SPECIES (LIKE TILAPIA)
keep_all = TRUE

if(keep_all){

  occ_data_sf <-  occ_data_sf_all
} else{

  ## Intersect range and keep points within range and buffer
  occ_data_sf <- st_intersection(occ_data_sf_all, st_buffer(range,0.1))

  ## Write outlying points to file
  occ_data_sf_all %>%
    filter(!(occurrence_id %in% occ_data_sf$occurrence_id)) %>%
    st_write(glue("data output/sdm data processing/{sppselect}/Outlier_occ_{sppselect}.shp"),
             delete_dsn = TRUE)
  }


# Spatial projections  ----------------------------------------------------
geo_proj <- query$Value[which(query$Input == "Geographic projection")]

aeaproj <- "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

latlongCRS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

if (geo_proj == "Yes") {
  occ_points <- occ_data_sf %>%
    st_transform(aeaproj)
} else {
  occ_points <- occ_data_sf %>%
    st_transform(latlongCRS)
}

## Remove duplicate coordinate values
occ_points <- occ_points %>%
  distinct(geometry, .keep_all = TRUE)

# Create background points ------------------------------------------------

## Read in catchment file
catchments <- st_read("data output/primary_catchment_latlong.shp", crs = latlongCRS) %>%
  st_transform(crs = aeaproj)

## Select catchments of interest
catch_vals <- query$Value[which(query$Input == "Catchments")]
regions <- str_c("Region ", str_split(catch_vals, pattern = ",")[[1]]) %>%
  str_replace(., "  ", " ")

catchment <- catchments %>%
  filter(name %in% regions) %>%
  st_union()

catchment %>%
  st_transform(crs = latlongCRS) %>%
  st_write(glue("data output/sdm data processing/{sppselect}/catchments_{sppselect}.shp"), delete_layer = TRUE)

# plot(catchment)

## Read in value of background points from query
n_points <- as.numeric(query$Value[which(query$Input == "Background points")])

## Randomly sample points within buffer (returns sfc object - use st_sf() to convert to sf)
bck_points <- st_sample(catchment, n_points) %>%
  st_sf()

## Write temporary files used in Python processing
bck_points %>%
  st_transform(crs = latlongCRS) %>%
  st_write("data output/temp_bck_dir/bck_points.shp", delete_layer = TRUE)

catchment %>%
  st_transform(crs = latlongCRS) %>%
  st_write("data output/temp_bck_dir/occ_buffer.shp", delete_layer = TRUE)

occ_points %>%
  st_transform(crs = latlongCRS) %>%
  st_write("data output/temp_bck_dir/occ_points.shp", delete_layer = TRUE)

## Run python background point processing
system('"C:/Program Files/ArcGIS/Pro/bin/Python/envs/arcgispro-py3/python" src/bck_points_processing.py',
       intern = TRUE)

## Read in adjusted background points
bck_points <- st_read("data output/temp_bck_dir/bck_points_updated.shp") %>%
  st_transform(aeaproj)

## Read in RSA grid
if (geo_proj == "Yes") {
  zagrid_sf <- readRDS("data input/zagrid_aea_sf.rds")
} else {
  ## Add non-projected sf ZA grid at some point
}

## Create a background point unique ID
bck_points <- mutate(bck_points, ID = str_c("B", 1:nrow(bck_points)))

## Assign a polygonID to each background point
bck_grid <- st_join(bck_points, zagrid_sf, join = st_intersects, left = FALSE)

## Group by polygons and keep first point (i.e. remove duplicate points)
bck_points <- bck_grid %>% #
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  ungroup

## Import RSA shapefile
if (geo_proj == "Yes") {
  za <- st_read("data input/RSA_fixed.shp", crs = latlongCRS) %>%
    st_transform(aeaproj)
} else {
  za <- st_read("data input/RSA_fixed.shp", crs = latlongCRS)
}

# Occurrence and BP plot checks --------------------------------------------

## Create a bounding box for map inset
occ_bbox <- st_make_grid(catchment, n = 1)

p1 <- ggplot() +
  geom_sf(data = za, fill = alpha("grey",0.5))+
  geom_sf(data = occ_points, size = 1.5)+
  geom_sf(data = occ_bbox, fill = NA, col = "dodgerblue", size = 1.3) +
  theme_bw()+
  theme(legend.text = element_text(size = 14))+
  ggtitle(glue("{query$Value[which(query$Input == 'Species')]} - {nrow(occ_points)} occurrence points"))

p2 <- ggplot() +
  geom_sf(data = za, fill = alpha("grey",0.5), size = 1.0)+
  geom_sf(data = zagrid_sf, fill = alpha("grey",0.3))+
  geom_sf(data = bck_points, size = 0.2, col = "red")+
  geom_sf(data = occ_points, size = 0.8, col = "blue")+
  coord_sf(xlim = c(st_bbox(occ_bbox)$xmin,st_bbox(occ_bbox)$xmax),
           ylim = c(st_bbox(occ_bbox)$ymin,st_bbox(occ_bbox)$ymax))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "dodgerblue", fill=NA, size=1.6))+
  ggtitle(query$Value[which(query$Input == "Species")])

p3 <- ggplot() +
  geom_sf(data = za, fill = alpha("grey",0.5), size = 1.0)+
  geom_sf(data = occ_points, size = 1.5)+
  geom_sf(data = catchment, fill = alpha("red",0.1))+
  coord_sf(xlim = c(st_bbox(occ_bbox)$xmin,st_bbox(occ_bbox)$xmax),
           ylim = c(st_bbox(occ_bbox)$ymin,st_bbox(occ_bbox)$ymax))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "dodgerblue", fill=NA, size=1.6))+
  theme(legend.text = element_text(size = 14))+
  ggtitle(glue("{query$Value[which(query$Input == 'Species')]}"))

p4 <- ggplot() +
  geom_sf(data = za, fill = alpha("grey",0.5), size = 1.0)+
  geom_sf(data = bck_points, size = 0.3)+
  coord_sf(xlim = c(st_bbox(occ_bbox)$xmin,st_bbox(occ_bbox)$xmax),
           ylim = c(st_bbox(occ_bbox)$ymin,st_bbox(occ_bbox)$ymax))+
  theme_bw()+
  theme(panel.border = element_rect(colour = "dodgerblue", fill=NA, size=1.6))+
  ggtitle(glue("{query$Value[which(query$Input == 'Species')]} - {nrow(bck_points)} background points"))

pdf(glue("data output/sdm data processing/{sppselect}/fig1.pdf"), width = 16, height = 9)
grid.arrange(grobs = list(p1,p2,p3,p4), ncol = 2) # Also see ggarrange
dev.off()

# p <- arrangeGrob(grobs = list(p1,p2,p3,p4), ncol = 2)
# ggsave(glue("data output/sdm data processing/{sppselect}/fig1.png"),
#        p,
#        width = 12, height = 6)

# Environmental data ------------------------------------------------------
env_layer_list <- read_xlsx(x, sheet = "env_var_list") %>%
  filter(.data[[glue("select - {sppselect}")]] == 1) %>% # .data[[]] is tidy eval for a string input
  dplyr::select(folder,layer) %>%
  mutate(path = glue("data input/Environmental data 30s reduced/{folder}/{layer}")) %>%
  dplyr::select(path) %>%
  pull

env_layer_list

hydro_layer_list <- read_xlsx(x, sheet = "hydro_var_list") %>%
  filter(.data[[glue("select - {sppselect}")]] == 1) %>% # .data[[]] is tidy eval for a string input
  dplyr::select(folder,layer) %>%
  mutate(path = glue("data input/Environmental data 30s reduced/{folder}/{layer}")) %>%
  dplyr::select(path) %>%
  pull

hydro_layer_list

envstack <- readAll(stack(env_layer_list))
hydrostack <- readAll(stack(hydro_layer_list))

envstack <- stack(envstack,hydrostack)
names(envstack)

# Project rasters ---------------------------------------------------------

## Project the rasters to AEA and crop (20% buffer added) to occurrence buffer
projenv_aea <- function(x){
  projection(x) <- latlongCRS
  projectRaster(x, crs=aeaproj)
}

projenv_latlong <- function(x){
  projection(x) <- latlongCRS
  projectRaster(x, crs=latlongCRS)
}

if (geo_proj == "Yes") {
  envstack <- stack(envstack)
  envstack <- stack(map(envstack@layers,projenv_aea))
  envstack <- readAll(stack(crop(envstack,extent(bck_points)*1.10))) # Add a 10% buffer
} else {
  envstack <- stack(envstack)
  envstack <- stack(map(envstack@layers,projenv_latlong))
  envstack <- readAll(stack(crop(envstack,extent(bck_points)*1.10)))
}

# Center and scale envstack -----------------------------------------------

## Create an option in the query file as to whether variables should be scaled?
# envstack <- scale(envstack)

# Collinearity of environemental predictors -------------------------------
cormatrix <- raster::layerStats(envstack, "pearson", na.rm = TRUE)

cormatrix <- as_tibble(cormatrix$`pearson correlation coefficient`) %>%
  mutate(var = colnames(cormatrix$`pearson correlation coefficient`)) %>%
  select(var, everything()) %>%
  write_csv(glue("data output/sdm data processing/{sppselect}/raster_corr_matrix.csv"))

## This method seems fairly robust and easy to understand
var_reduced <- usdm::vifcor(envstack, th=0.6)
var_reduced@results$Variables

# Alternative
# vif_reduced <- usdm::vifstep(envstack, th=10)

# Exclude the collinear variables that were identified in the previous step
envstack <- usdm::exclude(envstack,var_reduced)

usdm::vif(envstack) %>%
  as_tibble %>%
  write_csv(glue("data output/sdm data processing/{sppselect}/vif_check.csv")) %>%
  print(n = 20)

# Raster predictor plots --------------------------------------------------
enlen <- 1:nlayers(envstack)
plists <- split(enlen, ceiling(seq_along(enlen)/4))
plotras <- function(x){plot(envstack[[x]])}

pdf(glue("data output/sdm data processing/{sppselect}/fig3.pdf"), width = 16, height = 9)
map(plists, plotras)
dev.off()

# Write query file --------------------------------------------------------
copyfrom <- x
copyto <- glue("data output/sdm data processing/{sppselect}/")
file.copy(copyfrom, copyto, recursive=TRUE)

# Write occurrence points to shapefile -------------------------------------
occ_points %>%
  st_write(glue("data output/sdm data processing/{sppselect}/occ_points_{sppselect}.shp"),
           delete_dsn=TRUE)

# Write workspace ---------------------------------------------------------
envstack <- readAll(stack(envstack))

save(list = c("occ_points","bck_points","envstack"),
     file = glue("data output/sdm data processing/{sppselect}/sdm_input_data.RData"))

end <- Sys.time()

print(glue("SUCCESSFULLY COMPLETED - RUNTIME: {round(end-start,2)} {attr(end-start, 'units')}"))
system('CMD /C "ECHO GREAT SUCCESS. I LIKE. THE MODEL HAS FINISHED RUNNING! && PAUSE"',
       invisible=FALSE, wait=FALSE)

