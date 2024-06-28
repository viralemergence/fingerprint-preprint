
# =============== Generating climate change covariates from ERA5 temperature data ====================

library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
rr = raster::stack("C:/Users/roryj/Dropbox/Research/fingerprint/data/era5_land/era5land_temperature_monthly.grib")

# layer names (1950-2020, months 1-12)
names = data.frame(
  year = rep(1950:2020, each=12),
  month = rep(1:12, length(1950:2020))
)




# =================== temperature change: change in mean between baseline period (1950-70) and study period (2001-2020) =======================

# annual mean temp  
ff = list.files("E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", pattern=".tif", full.names=TRUE)
ff = ff[ -grep(".xml", ff)]
ff = ff[ grep("meantemperature", ff)]

# 1. baseline
# ie. mean and variability in annual temperatures across years in baseline period
ff1 = ff[ 1:21]
t_baseline = stack(ff1)
t_longtermaverage = mean(t_baseline)

# 2. tmean during focal period
ff2 = ff[ 52:71]
tmean_i = stack(ff2)
t_focal = mean(tmean_i)

# calculate change
t_diff = t_focal - t_longtermaverage
writeRaster(t_diff, file="E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temperature_meanchange_baselinetopresent.tif", format="GTiff", overwrite=TRUE)

# save mean temperature in study period
ff2 = ff[ 51:71]
tmean_i = stack(ff2)
t_focal = mean(tmean_i)
writeRaster(t_focal, file="E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temperature_annualmean_20002020.tif", format="GTiff", overwrite=TRUE)


# ne = rnaturalearth::ne_coastline(returnclass = "sf")
# anom_plot = t_diff %>%
#   as.data.frame(xy=TRUE) %>%
#   dplyr::rename("mean_anomaly"=3) %>%
#   ggplot() + 
#   maptheme + 
#   geom_raster(aes(x, y, fill=mean_anomaly)) + 
#   scale_fill_gradient2(high = "darkred", low="darkblue", mid = "white", midpoint=0, na.value="white", name="Anomaly") + 
#   geom_sf(data=ne, fill=NA, color="grey50", size=0.2) +
#   ggtitle("Standardised mean temperature anomaly 1991-2020 compared to baseline 1950-1990") + 
#   theme(plot.title=element_text(size=13, hjust=0.5))
# ggsave(anom_plot, file="./era5land_tempanomaly_19812020.png", device="png", units="in", dpi=600, width=8, height=4)


# =================== temperature anomaly: mean anomaly during 1991-2020 compared to 1951-1990 =======================

# annual anomaly per pixel i and year y calculated as 
# (tmean_i,y - longtermaverage_i) / var(tmean_i,y)

# can then calculate mean/cumulative anomaly
# baseline period for longtermaverage == 1951-1990
# period for calculating anomaly == 1991-2020

# annual mean temp  
ff = list.files("E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", pattern=".tif", full.names=TRUE)
ff = ff[ -grep(".xml", ff)]
ff = ff[ grep("meantemperature", ff)]

# 1. longtermaverage: stack and calculate long-term mean across 1950-1990 and long term sd
# ie. mean and variability in annual temperatures across years in baseline period
ff1 = ff[ 1:41]
t_baseline = stack(ff1)
t_longtermaverage = mean(t_baseline)
t_longtermsd = raster::calc(t_baseline, fun = sd)

# 2. annual tmean during focal period (1991-2020)
ff2 = ff[ 42:71]
tmean_i = stack(ff2)

# calculate anomaly
anom_i = tmean_i - t_longtermaverage
anom_mean = mean(anom_i)

# calculate standardised anomaly (mean anomaly during focal period divided by the variance during the baseline period)
mean_anomaly = anom_mean / t_longtermsd
writeRaster(mean_anomaly, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temperature_meananomaly_19912020.tif", format="GTiff", overwrite=TRUE)


ne = rnaturalearth::ne_coastline(returnclass = "sf")
anom_plot = mean_anomaly %>%
  as.data.frame(xy=TRUE) %>%
  dplyr::rename("mean_anomaly"=3) %>%
  ggplot() + 
  maptheme + 
  geom_raster(aes(x, y, fill=mean_anomaly)) + 
  scale_fill_gradient2(high = "darkred", low="darkblue", mid = "white", midpoint=0, na.value="white", name="Anomaly") + 
  geom_sf(data=ne, fill=NA, color="grey50", size=0.2) +
  ggtitle("Standardised mean temperature anomaly 1991-2020 compared to baseline 1950-1990") + 
  theme(plot.title=element_text(size=13, hjust=0.5))
ggsave(anom_plot, file="./era5land_tempanomaly_19812020.png", device="png", units="in", dpi=600, width=8, height=4)









# ================ OLD CODE ==================


# # ============ change in mean temp ===============
# 
# # mean temp
# for(i in 1950:2020){
#   print(i)
#   yx = rr[[ which(names$year == i) ]]
#   mt_x = mean(yx) - 273.15 # convert to celcius
#   names(mt_x) = paste("era5land_meantemperature_", i, sep="")
#   writeRaster(mt_x, file=paste("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", names(mt_x), sep=""), format="GTiff", overwrite=TRUE)
# }
# rm(rr)
# 
# # read in and use linear models to estimate slope across years 1980 - 2020
# ff = list.files("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", pattern=".tif", full.names=TRUE)
# ff = ff[ -grep(".xml", ff)]
# ff = ff[ 32:71]
# tras = stack(ff)
# 
# # extract
# temp_dat = as.data.frame(tras, xy=TRUE)
# temp_dat$cellid = 1:nrow(temp_dat)
# temp_dat = temp_dat[ !is.na(temp_dat$era5land_meantemperature_1981), ]
# 
# t2 = temp_dat %>%
#   dplyr::select(-x, -y) 
# names(t2)[ 1:40 ] = paste("t", 1981:2020, sep="")
# t2 = t2 %>%
#   reshape2::melt(id.vars=c("cellid")) %>%
#   dplyr::mutate(variable = substr(variable, 2, 7))
# rm(temp_dat)
# 
# # calculate slope
# slopes = t2 %>%
#   dplyr::arrange(cellid, variable) %>%
#   dplyr::group_by(cellid) %>%
#   dplyr::summarise(
#     slope = as.vector(lm(value ~ as.numeric(variable))$coefficients[2])
#   )
# write.csv(slopes, "./slopes_eratemp.csv", row.names=FALSE)
# 
# # create raster of slopes
# tslopes = raster(tras)
# tvals = as.data.frame(tslopes, xy=TRUE)
# tvals$cellid = 1:nrow(tvals)
# tvals = left_join(tvals, slopes)
# tvals = tvals %>% dplyr::arrange(cellid)
# values(tslopes) = tvals$slope
# 
# as.data.frame(tslopes, xy=TRUE) %>%
#   ggplot() + 
#   geom_raster(aes(x, y, fill=layer)) + 
#   scale_fill_gradient2() + 
#   coord_fixed()
# 
# writeRaster(tslopes, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_slope_19812020.tif", format="GTiff")
# 
# 
# 
# # ================ temperature of the warmest month =================
# 
# # temperature of the warmest month
# tras = stack()
# for(i in 1950:2020){
#   print(i)
#   yx = rr[[ which(names$year == i) ]]
#   mt_x = max(yx)
#   mt_x = mt_x - 273.15 # convert to celcius
#   names(mt_x) = paste("era5land_tempwarmest_", i, sep="")
#   tras = stack(tras, mt_x)
#   #writeRaster(mt_x, file=paste("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", names(mt_x), sep=""), format="GTiff", overwrite=TRUE)
# }
# rm(rr)
# 
# writeRaster(tras, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/twarmest_annual", format="raster", overwrite=TRUE)
# 
# 
# # read in and use linear models to estimate slope across years 1980 - 2020
# ff = list.files("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", pattern=".tif", full.names=TRUE)
# ff = ff[ -grep(".xml", ff)]
# ff = ff[ grep("tempwarmest", ff)]
# ff = ff[ 32:71]
# tras = stack(ff)
# 
# tt = tras
# tras = tras[[32:71]]
# 
# # extract
# temp_dat = as.data.frame(tras, xy=TRUE)
# temp_dat$cellid = 1:nrow(temp_dat)
# temp_dat = temp_dat[ !is.na(temp_dat$era5land_tempwarmest_1982), ]
# 
# t2 = temp_dat %>%
#   dplyr::select(-x, -y) 
# names(t2)[ 1:40 ] = paste("t", 1981:2020, sep="")
# t2 = t2 %>%
#   reshape2::melt(id.vars=c("cellid")) %>%
#   dplyr::mutate(variable = substr(variable, 2, 7))
# rm(temp_dat)
# 
# # calculate slope
# slopes = t2 %>%
#   dplyr::arrange(cellid, variable) %>%
#   dplyr::group_by(cellid) %>%
#   dplyr::summarise(
#     slope = as.vector(lm(value ~ as.numeric(variable))$coefficients[2])
#   )
# write.csv(slopes, "./slopes_eratemp_twarmestmonth.csv", row.names=FALSE)
# 
# # create raster of slopes
# tslopes = raster(tras)
# tvals = as.data.frame(tslopes, xy=TRUE)
# tvals$cellid = 1:nrow(tvals)
# tvals = left_join(tvals, slopes)
# tvals = tvals %>% dplyr::arrange(cellid)
# values(tslopes) = tvals$slope
# 
# as.data.frame(tslopes, xy=TRUE) %>%
#   ggplot() + 
#   geom_raster(aes(x, y, fill=layer)) + 
#   scale_fill_gradient2() + 
#   coord_fixed()
# 
# writeRaster(tslopes, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/twarmest_slope_19812020.tif", format="GTiff")
# plot(tslopes)
# 
# 
# # ================ temperature of the coolest month =================
# 
# # temperature of the warmest month
# tras = stack()
# for(i in 1950:2020){
#   print(i)
#   yx = rr[[ which(names$year == i) ]]
#   mt_x = min(yx)
#   mt_x = mt_x - 273.15 # convert to celcius
#   names(mt_x) = paste("era5land_tempcoolest_", i, sep="")
#   tras = stack(tras, mt_x)
#   #writeRaster(mt_x, file=paste("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/temp_annual/", names(mt_x), sep=""), format="GTiff", overwrite=TRUE)
# }
# rm(rr)
# 
# tt = tras
# tras = tras[[32:71]]
# 
# # extract
# temp_dat = as.data.frame(tras, xy=TRUE)
# temp_dat$cellid = 1:nrow(temp_dat)
# temp_dat = temp_dat[ !is.na(temp_dat$era5land_tempcoolest_1981), ]
# 
# t2 = temp_dat %>%
#   dplyr::select(-x, -y) 
# names(t2)[ 1:40 ] = paste("t", 1981:2020, sep="")
# t2 = t2 %>%
#   reshape2::melt(id.vars=c("cellid")) %>%
#   dplyr::mutate(variable = substr(variable, 2, 7))
# rm(temp_dat)
# 
# # calculate slope
# slopes = t2 %>%
#   dplyr::arrange(cellid, variable) %>%
#   dplyr::group_by(cellid) %>%
#   dplyr::summarise(
#     slope = as.vector(lm(value ~ as.numeric(variable))$coefficients[2])
#   )
# write.csv(slopes, "./slopes_eratemp_tcoolestmonth.csv", row.names=FALSE)
# 
# # create raster of slopes
# tslopes = raster(tras)
# tvals = as.data.frame(tslopes, xy=TRUE)
# tvals$cellid = 1:nrow(tvals)
# tvals = left_join(tvals, slopes)
# tvals = tvals %>% dplyr::arrange(cellid)
# values(tslopes) = tvals$slope
# 
# # as.data.frame(tslopes, xy=TRUE) %>%
# #   ggplot() + 
# #   geom_raster(aes(x, y, fill=layer)) + 
# #   scale_fill_gradient2() + 
# #   coord_fixed()
# 
# writeRaster(tslopes, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/tcoolest_slope_19812020.tif", format="GTiff", overwrite=TRUE)
# plot(tslopes)
