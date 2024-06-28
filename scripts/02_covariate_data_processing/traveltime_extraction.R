


# ===================== Extract population within specified health travel time thresholds ====================

setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")

# read in rasters
# hca_m = terra::rast("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/healthcare_travel/2020_motorized_travel_time_to_healthcare/2020_motorized_travel_time_to_healthcare.geotiff")
# hca_w = terra::rast("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/healthcare_travel/2020_walking_only_travel_time_to_healthcare/2020_walking_only_travel_time_to_healthcare.geotiff")
# wp = terra::rast("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_1km_Aggregated.tif")
# 
# # crop resample wp to match hca
# wp = terra::resample(wp, hca_m, method="near")
# raster::writeRaster(raster::raster(wp), file="C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/pop2020_resampled_htt.tif", format="GTiff")

# worldpop
wp = terra::rast("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/pop2020_resampled_htt.tif")

# create rasters of thresholded tthc
# for(thresh in c(30, 60, 90, 120)){
#   
#   tm = hca_m < thresh
#   tw = hca_w < thresh
#   tm = tm * wp
#   tw = tw * wp
#   names(tm) = paste("tt_motorized_", thresh, sep="")
#   names(tw) = paste("tt_walking_", thresh, sep="")
#   
#   writeRaster(raster::raster(tm), file=paste(names(tm), ".tif", sep=""), format="GTiff")
#   writeRaster(raster::raster(tw), file=paste(names(tw), ".tif", sep=""), format="GTiff")
#   
# }


# read in rasters
rr = terra::rast(list.files(pattern=".tif", full.names=TRUE))

# for each disease read study area
ll = rev(list.files("./output/model_outputs/disease_models/aug23/", full.names=TRUE, pattern=".R"))

for(l in ll){
  
  load(l)
  
  dz_l = model_obj$data$Disease[1]
  print(dz_l)
  
  # create mask
  print("cropping...")
  sa = model_obj$study_area %>% dplyr::mutate(dummy=1)
  rr_l = terra::crop(rr, sa)
  wp_l = terra::crop(wp, sa)
  names(wp_l) = "total_pop"
  rr_l = terra::rast(list(rr_l, wp_l))
  layer_names = names(rr_l)
  
  # as dengue is within latitudinal bounds just crop; otherwise mask
  if(!dz_l %in% c("Dengue")){
    print("rasterising...")
    mask = fasterize::fasterize(sa, raster::raster(rr_l[[1]]), field = "dummy")
    rr_l = rr_l * terra::rast(mask)
  }
  names(rr_l) = layer_names

  print("summarising...")
  result = data.frame(
    variable = names(rr_l),
    value = terra::global(rr_l, fun="sum", na.rm=TRUE)
  ) %>%
    dplyr::rename("value" = sum) %>%
    dplyr::mutate(total_pop = value[ variable == "total_pop" ],
                  prop_pop = value / total_pop,
                  type = ifelse(grepl("motorized", variable), "motorized", "walking"),
                  mins_tt = unlist( lapply( strsplit( variable, "_"), "[", 3)),
                  Disease = dz_l) %>%
    dplyr::filter(variable != "total_pop")
  row.names(result) = c()
  
  write.csv(result, paste("./output/model_outputs/tthc_per_disease/", dz_l, ".csv", sep=""), row.names=FALSE)
  
  # remove
  rm(rr_l)
  rm(mask)
  rm(wp_l)
  
}

# check outputs
res = do.call(
  rbind.data.frame, 
  lapply(list.files("./output/model_outputs/tthc_per_disease/", pattern=".csv", full.names=TRUE), read.csv)
)

# plot
res %>% 
  dplyr::filter(variable == "tt_motorized_120") %>%
  ggplot() + 
  geom_point(aes(Disease, (1-prop_pop)*100), size=3) + 
  coord_flip()

res %>% 
  dplyr::filter(variable == "tt_walking_120") %>%
  ggplot() + 
  geom_point(aes(Disease, (1-prop_pop)*100), size=3) + 
  coord_flip()

res %>%
  dplyr::filter() %>%
  dplyr::filter(variable == "tt_motorized_120") %>%
  ggplot() + 
  geom_histogram(aes((total_pop - value)/10^6))
