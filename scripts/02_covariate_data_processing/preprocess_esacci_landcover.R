


# ======= derives forest, urbanisation and agricultural expansion layers from ESA-CCI land cover ======

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")

# ESA CCI land cover 2000 and 2018
# es = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7.tif")
ee = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/C3S-LC-L4-LCCS-Map-300m-P1Y-2018-v2.1.1.nc")

# create mask raster for NAs on ocean
# mask = es == 210
# writeRaster(mask, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_watermask.tif", format="GTiff")
mask = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_watermask.tif")



# ========== 1. agriculture ==================

# create agriculture layers per epoch
# ag1 = es %in% c(10, 11, 12, 20, 30, 40)
# writeRaster(ag1, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agri_2000.tif", format="GTiff")
# ag2 = ee %in% c(10, 11, 12, 20, 30, 40)
# writeRaster(ag2, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agri_2018.tif", format="GTiff")

# read in layers per epoch
ag1 = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agri_2000.tif")
ag2 = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agri_2018.tif")

# expansion raster
agx = ag1 == 0 & ag2 == 1
writeRaster(agx, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agriexpansion_20002018.tif", format="GTiff")

# mask rasters to account for ocean
agx = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agriexpansion_20002018.tif")
agx = raster::mask(agx, mask, maskvalue=1)
writeRaster(agx, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agriexpansion_masked_20002018.tif", format="GTiff")

# mask rasters to account for ocean
ag1 = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agri_2018.tif")
ag1 = raster::mask(ag1, mask, maskvalue=1)
writeRaster(ag1, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agrimasked_2018.tif", format="GTiff")


# ========== 2. urban ==================

# create urban layers per epoch
# ur1 = es %in% c(190)
# writeRaster(ur1, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urban_2000.tif", format="GTiff")
# ur2 = ee %in% c(190)
# writeRaster(ur2, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urban_2018.tif", format="GTiff")

# read in layers per epoch
ur1 = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urban_2000.tif")
ur2 = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urban_2018.tif")

xx = extent(3, 4.25, 6, 8)
ur1 = crop(ur1, xx)
# ur2 = crop(ur2, xx)
# urx = ur1 == 0 & ur2 == 1
# 
# u1 = ur1 %>% as.data.frame(xy = TRUE)
# u2 = ur2 %>% as.data.frame(xy = TRUE) 
# u1$esacci_urban_2000[ u1$esacci_urban_2000 == 0 ] = NA
# u2$esacci_urban_2018[ u2$esacci_urban_2018 == 0 ] = NA
# ggplot() + 
#   geom_raster(data=u1, aes(x=x, y=y), fill="white") + 
#   geom_raster(data=u2[ !is.na(u2$esacci_urban_2018), ], aes(x=x, y=y), fill="red", alpha=0.7) + 
#   geom_raster(data=u1[ !is.na(u1$esacci_urban_2000), ], aes(x=x, y=y), fill="blue", alpha=0.6) + 
#   maptheme +
#   coord_fixed()

# expansion raster
urx = ur1 == 0 & ur2 == 1
writeRaster(urx, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanexpansion_20002018.tif", format="GTiff")

# mask rasters to account for ocean
urx = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanexpansion_20002018.tif")
urx = raster::mask(urx, mask, maskvalue=1)
writeRaster(urx, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanexpansion_masked_20002018.tif", format="GTiff")

# mask rasters to account for ocean
ur1 = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urban_2018.tif")
ur1 = raster::mask(ur1, mask, maskvalue=1)
writeRaster(ur1, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanmasked_2018.tif", format="GTiff")


# ============ 3. forest cover =============

fo1 = ee %in% c(50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 160, 170)
writeRaster(fo1, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_forest_2018.tif", format="GTiff")
fo1 = raster::mask(fo1, mask, maskvalue=1)
writeRaster(fo1, "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_forestmasked_2018.tif", format="GTiff")


