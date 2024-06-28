

# ================ EXTRACT COVARIATES FOR EVERY USA ADMIN2 ===================

# timesaver when modelling later

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# counties shapefile
shp = sf::st_read("./data/spillovers/Lyme/us_cdc/cb_2018_us_county_500k.shp") %>%
  dplyr::mutate(CTYCODE = as.integer(COUNTYFP),
                STCODE = as.integer(STATEFP)) 
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
)
shp$ADMcode = paste(shp$STCODE, shp$CTYCODE, sep="_")

# only mainland contiguous us excluding alaska
study_area <- st_read("./data/shapefiles/cb_2018_us_state_5m.shp") %>% 
  dplyr::filter(NAME %in% c("Hawaii", "Guam", "American Samoa", "Puerto Rico", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands", "Alaska"))
fip = as.numeric(unique(study_area$STATEFP))
shp = shp %>%
  dplyr::filter(!STCODE %in% fip) %>%
  st_transform(crs = 4326) 

# extract covariates
print("Extracting covariate data")
covs_df = read.csv("./scripts/03_modelling/covariate_lookup.csv")
shpx = shp %>% st_drop_geometry() %>% dplyr::select(ADMcode, NAME, Longitude, Latitude) %>% dplyr::rename("Admin_name"=NAME)
for(cc in covs_df$cov[ 20:21 ]){
  print(cc)
  extr_x = extractCovariateVals(shp, covariate_name=cc) # function stored in '00_covar_extraction_funcs.R'
  shpx = cbind(shpx, extr_x)
}

filename = "./output/model_df/usa_admin2_covars.csv"
write.csv(shpx, filename, row.names=FALSE)




