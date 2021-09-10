README
================

# threatened-fish-sdms

# Species distribution models for South Africa’s Threatened Freshwater Fishes

## Files required in input folder

``` r
list.files("data input")
```

    ##  [1] "Environmental data 30s reduced"                 
    ##  [2] "Fish Distribution Map_IUCN"                     
    ##  [3] "Occurrence data_threatened only_08 Dec 2020.csv"
    ##  [4] "Queries"                                        
    ##  [5] "RSA_fixed.cpg"                                  
    ##  [6] "RSA_fixed.dbf"                                  
    ##  [7] "RSA_fixed.prj"                                  
    ##  [8] "RSA_fixed.sbn"                                  
    ##  [9] "RSA_fixed.sbx"                                  
    ## [10] "RSA_fixed.shp"                                  
    ## [11] "RSA_fixed.shp.xml"                              
    ## [12] "RSA_fixed.shx"                                  
    ## [13] "zagrid_aea_sf.rds"

``` r
list.files("data input/Queries")
```

    ## [1] "SDM_query_Labeo rubromaculatus.xlsx"

``` r
list.files("data input/Environmental data 30s reduced")
```

    ##  [1] "Bioclim layers"       "Envirem layers"       "Fragmentation layers"
    ##  [4] "Hydrological rasters" "Landcover layers"     "NDVI layers"         
    ##  [7] "Soils layers"         "Topography layers"    "Vegmap layers"       
    ## [10] "Wetland layers"

## Files required in output folder

``` r
list.files("data output")
```

    ## [1] "occ_data_clean.RData"          "primary_catchment_latlong.dbf"
    ## [3] "primary_catchment_latlong.prj" "primary_catchment_latlong.shp"
    ## [5] "primary_catchment_latlong.shx" "sdm BART results"             
    ## [7] "sdm data processing"           "temp_bck_dir"
