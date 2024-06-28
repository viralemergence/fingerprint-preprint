

# ================ EXTRACT COVARIATES FOR EVERY BRAZIL ADMIN2 ===================

# timesaver when modelling later

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# municipios shapefile
shp = sf::st_read("./data/spillovers/Chagas/brazil_moh/Brazil_shp_harm_2022.shp") %>%
  dplyr::mutate(ADMcode = IBGE6)

#
disease_name = "brazil_all"

# extract covariates
print("Extracting covariate data")
covs_df = read.csv("./scripts/03_modelling/covariate_lookup.csv")
shpx = shp %>% st_drop_geometry() %>% dplyr::select(code_mn, IBGE6, ADMcode, name_mn)
for(cc in covs_df$cov[ 7:nrow(covs_df)  ]){
  print(cc)
  extr_x = extractCovariateVals(shp, covariate_name=cc) # function stored in '00_covar_extraction_funcs.R'
  shpx = cbind(shpx, extr_x)
}

filename = "./output/model_df/brazil_admin2_covars.csv"
write.csv(shpx, filename, row.names=FALSE)




