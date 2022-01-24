threatened-fish-sdms
================

# Species distribution models for South Africa’s Threatened Freshwater Fishes

## Files required in input folder

``` r
list.files("data input")
```

    ##  [1] "~$clipping_parameters_EST.xlsx"                           
    ##  [2] "Catchment Layers"                                         
    ##  [3] "clipping_parameters_EST.xlsx"                             
    ##  [4] "Environmental data 30s reduced"                           
    ##  [5] "Fish Distribution Map_IUCN"                               
    ##  [6] "List of species included_Clipping Extent for species.xlsx"
    ##  [7] "Occurrence data_threatened only_08 Dec 2020.csv"          
    ##  [8] "Occurrence data_Threatened.csv"                           
    ##  [9] "Occurrence data_Threatened_14 Sep 2021.csv"               
    ## [10] "Queries"                                                  
    ## [11] "River FEPAs"                                              
    ## [12] "Rivers"                                                   
    ## [13] "RSA_fixed.cpg"                                            
    ## [14] "RSA_fixed.dbf"                                            
    ## [15] "RSA_fixed.prj"                                            
    ## [16] "RSA_fixed.sbn"                                            
    ## [17] "RSA_fixed.sbx"                                            
    ## [18] "RSA_fixed.shp"                                            
    ## [19] "RSA_fixed.shp.CTDH-LT.14876.47212.sr.lock"                
    ## [20] "RSA_fixed.shp.CTDH-LT.16600.47212.sr.lock"                
    ## [21] "RSA_fixed.shp.CTDH-LT.2904.47212.sr.lock"                 
    ## [22] "RSA_fixed.shp.CTDH-LT.51596.47212.sr.lock"                
    ## [23] "RSA_fixed.shp.xml"                                        
    ## [24] "RSA_fixed.shx"                                            
    ## [25] "Threatened Species Catchments_09 Sep 2021.xlsx"           
    ## [26] "zagrid_aea_sf.rds"

``` r
list.files("data input/Catchment Layers")
```

    ## [1] "primary_catchment"        "primary_catchment_area"  
    ## [3] "quaternary_catchment"     "quinary_catchment"       
    ## [5] "secondary_catchment_area" "tertiary_catchment_area"

``` r
list.files("data input/River FEPAs")
```

    ##  [1] "NFEPA_metadata_riverFEPAs.pdf" "NFEPA_riverFEPAs.zip"         
    ##  [3] "River Fepas.lyr"               "River_FEPAs.dbf"              
    ##  [5] "River_FEPAs.prj"               "River_FEPAs.sbn"              
    ##  [7] "River_FEPAs.sbx"               "River_FEPAs.shp"              
    ##  [9] "River_FEPAs.shp.xml"           "River_FEPAs.shx"              
    ## [11] "River_FEPAs_dd.docx"

``` r
list.files("data input/Queries")
```

    ##  [1] "SDM_query_Austroglanis barnardi.xlsx"     
    ##  [2] "SDM_query_Chetia brevis.xlsx"             
    ##  [3] "SDM_query_Chiloglanis bifurcus.xlsx"      
    ##  [4] "SDM_query_Chiloglanis emarginatus.xlsx"   
    ##  [5] "SDM_query_Ctenopoma multispine.xlsx"      
    ##  [6] "SDM_query_Enteromius gurneyi.xlsx"        
    ##  [7] "SDM_query_Enteromius treurensis.xlsx"     
    ##  [8] "SDM_query_Labeo rubromaculatus.xlsx"      
    ##  [9] "SDM_query_Labeo seeberi.xlsx"             
    ## [10] "SDM_query_Marcusenius caudisquamatus.xlsx"
    ## [11] "SDM_query_Oreochromis mossambicus.xlsx"   
    ## [12] "SDM_query_Pseudobarbus afer.xlsx"         
    ## [13] "SDM_query_Pseudobarbus asper.xlsx"        
    ## [14] "SDM_query_Pseudobarbus burgi.xlsx"        
    ## [15] "SDM_query_Pseudobarbus capensis.xlsx"     
    ## [16] "SDM_query_Pseudobarbus erubescens.xlsx"   
    ## [17] "SDM_query_Pseudobarbus phlegethon.xlsx"   
    ## [18] "SDM_query_Pseudobarbus quathlambae.xlsx"  
    ## [19] "SDM_query_Pseudobarbus senticeps.xlsx"    
    ## [20] "SDM_query_Pseudobarbus skeltoni.xlsx"     
    ## [21] "SDM_query_Pseudobarbus swartzi.xlsx"      
    ## [22] "SDM_query_Pseudobarbus verloreni.xlsx"    
    ## [23] "SDM_query_Sandelia bainsii.xlsx"          
    ## [24] "SDM_query_Serranochromis meridianus.xlsx" 
    ## [25] "SDM_query_Silhouettea sibayi.xlsx"

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
    ## [5] "primary_catchment_latlong.shx" "screening tool layers"        
    ## [7] "sdm BART results"              "sdm data processing"          
    ## [9] "temp_bck_dir"

``` r
list.files("data output/screening tool layers")
```

    ##  [1] "FWF_high.dbf"      "FWF_high.prj"      "FWF_high.shp"     
    ##  [4] "FWF_high.shx"      "FWF_medium.dbf"    "FWF_medium.prj"   
    ##  [7] "FWF_medium.shp"    "FWF_medium.shx"    "FWF_very_high.dbf"
    ## [10] "FWF_very_high.prj" "FWF_very_high.shp" "FWF_very_high.shx"
    ## [13] "high"              "medium"            "very high"
