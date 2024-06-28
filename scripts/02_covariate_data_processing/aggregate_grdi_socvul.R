


# =============== aggregate global gridded relative deprivation index ================

# (highly collinear with urbanisation at the grid cell level)
# aggregate to 20km grid cells to capture coarse pattern while reducing this issue

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)

# soc vul
sv = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/socialvulnerability_global/povmap-grdi-v1_BETA_September2022/povmap-grdi-v1_SEP2022-BETA/povmap-grdi-v1-Sep2022.tif")

# aggregate
sv_m = raster::aggregate(sv, fact=10, fun=mean, na.rm=TRUE)
writeRaster(sv_m, file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/socialvulnerability_global/socialvul_mean_fact10.tif", format="GTiff")

# aggregate to 20km grids
sv_m2 = raster::aggregate(sv_m, fact=2, fun=mean, na.rm=TRUE)
writeRaster(sv_m2, file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/socialvulnerability_global/socialvul_mean_fact20.tif", format="GTiff")

# aggregate to sd
sv_sd = raster::aggregate(sv, fact=20, fun=sd, na.rm=TRUE)
writeRaster(sv_sd, file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/socialvulnerability_global/socialvul_sd_fact20.tif", format="GTiff")
