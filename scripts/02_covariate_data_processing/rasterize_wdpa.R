
# ============= create protected areas raster from World Database of Protected Areas ==================

# reads shapefile polygon PAs
# converts to a binary raster of "protected/not-protected"

library(sf)
library(fasterize)
library(raster)

# 1km template
template_ras = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/mining/global_miningarea_v2_30arcsecond.tif")

# shapefiles
shp1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_WDOECM_Sep2022_Public_all_shp/WDPA_WDOECM_Sep2022_Public_all_shp_0/WDPA_WDOECM_Sep2022_Public_all_shp-polygons.shp") %>%
  dplyr::mutate(id = 1)
ff1 = fasterize::fasterize(shp1, template_ras, field="id")

shp2 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_WDOECM_Sep2022_Public_all_shp/WDPA_WDOECM_Sep2022_Public_all_shp_1/WDPA_WDOECM_Sep2022_Public_all_shp-polygons.shp") %>%
  dplyr::mutate(id = 1)
ff2 = fasterize::fasterize(shp2, template_ras, field="id")

shp3 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_WDOECM_Sep2022_Public_all_shp/WDPA_WDOECM_Sep2022_Public_all_shp_2/WDPA_WDOECM_Sep2022_Public_all_shp-polygons.shp") %>%
  dplyr::mutate(id = 1)
ff3 = fasterize::fasterize(shp3, template_ras, field="id")

# create raster 
ff = ff1 | ff2 | ff3
writeRaster(ff, file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_Sep2022_1km_ras.tif", format="GTiff")

# set NA values to zero
ff = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_Sep2022_1km_ras.tif")
ff = raster::reclassify(x = ff, rcl = cbind(NA, 0))
writeRaster(ff, file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_Sep2022_1km_ras_zeroes.tif", format="GTiff")

# create ocean mask from template raster
mask = !is.na(template_ras)
ff_mk = raster::mask(ff, mask, maskvalue=0)
writeRaster(ff_mk, file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/protected_areas/WDPA_Sep2022_1km_ras_final.tif", format="GTiff")

