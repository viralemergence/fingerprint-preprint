

# =============== Generating climate change covariates from ERA5 precip data ====================

library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
rr = raster::stack("C:/Users/roryj/Dropbox/Research/fingerprint/data/era5_land/era5land_precip_monthly.grib")

# layer names (1950-2020, months 1-12)
# precip comes in mean daily precip for the month so need to multiply by days in month for total
names = data.frame(
  year = rep(1950:2020, each=12),
  month = rep(1:12, length(1950:2020))
)
names$days_in_month = lubridate::days_in_month(as.Date(paste(names$year, names$month, "01", sep="-")))

# annual total precip
for(i in 1958:2020){
  print(i)
  yx = rr[[ which(names$year == i) ]]
  yx = yx * 1000 * names$days_in_month[ names$year == i ]  # convert to total precip mm 
  yx = sum(yx)
  names(yx) = paste("era5land_totalprecipmm_", i, sep="")
  writeRaster(yx, file=paste("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_annual/", names(yx), sep=""), format="GTiff", overwrite=TRUE)
}
rm(rr)

plot(crop(yx, extent(-82, -72, -10, 10)))




# -------------------- precip change: difference between historical baseline and focal period ----------------

# precips
ff = list.files("E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_annual/", pattern=".tif", full.names=TRUE)
ff = ff[ -grep(".xml", ff)]

# define baseline period (1950-70)
ff1 = ff[ 1:21]
p_baseline = stack(ff1)
p_longtermaverage = mean(p_baseline)

# replace extreme high values with NAs (mainly ecuador/colombia) 
# should explore how to do this better, but not a major issue as do not overlap with our spillover data
vv = values(p_longtermaverage) 
vv[ vv> 10000] = NA
values(p_longtermaverage) = vv

# define study period
ff1 = ff[ 52:71]
p_focal = stack(ff1)
p_focal = mean(p_focal)

# difference
p_diff = p_focal - p_longtermaverage
writeRaster(p_diff, file="E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_meanchange_baselinetopresent.tif", format="GTiff", overwrite=TRUE)


# save annual mean precip 2000-2020
ff1 = ff[ 51:71]
p_focal = stack(ff1)
p_focal = mean(p_focal)
writeRaster(p_focal, file="E:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_meanannual_20002020.tif", format="GTiff", overwrite=TRUE)




# --------------- precipitation anomalies per-year: wetness, dryness and average standardised -------------------

ff = list.files("D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_annual/", pattern=".tif", full.names=TRUE)
ff = ff[ -grep(".xml", ff)]

# 1. longtermaverage: stack and calculate long-term mean across 1950-1990 and long term sd
# ie. mean and variability in annual precipitation across years in baseline period
ff1 = ff[ 1:41]
p_baseline = stack(ff1)
p_longtermaverage = mean(p_baseline)
p_longtermsd = raster::calc(p_baseline, fun = sd)

# replace erroneous areas with NAs (mainly ecuador/colombia) - need to explore how to do this better
vv = values(p_longtermaverage) 
vv[ vv> 12000] = NA
values(p_longtermaverage) = vv

# 2. annual precipitation during focal period (1981-2020)
ff2 = ff[ 42:71]
precip_i = stack(ff2)

# calculate annual anomaly and standardise by sd
anom_i = precip_i - p_longtermaverage
anom_i_standardised = anom_i / p_longtermsd

# Metric 1: proportion of years with standardised positive anomaly > 1 sd ("wetness")
wetness = sum( anom_i_standardised > 1 ) / nlayers(anom_i_standardised)
writeRaster(wetness, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_wetnessanomaly_19912020.tif", format="GTiff", overwrite=TRUE)

# Metric 2: number of years with standardised negative anomaly > 1 sd ("dryness")
dryness = sum( anom_i_standardised < -1 ) / nlayers(anom_i_standardised)
writeRaster(dryness, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_drynessanomaly_19912020.tif", format="GTiff", overwrite=TRUE)

# Metric 3: mean standardised anomaly
mean_anomaly = mean(anom_i_standardised)
writeRaster(mean_anomaly, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precip_meananomaly_19912020.tif", format="GTiff", overwrite=TRUE)


ne = rnaturalearth::ne_coastline(returnclass = "sf")
mean_anomaly %>%
  as.data.frame(xy=TRUE) %>%
  dplyr::rename("mean_anomaly"=3) %>%
  ggplot() +
  maptheme +
  geom_raster(aes(x, y, fill=mean_anomaly)) +
  scale_fill_gradient2(high = "darkred", low="darkblue", mid = "white", midpoint=0, na.value="white", name="Anomaly") +
  geom_sf(data=ne, fill=NA, color="grey50", size=0.2) +
  ggtitle("Standardised mean precip anomaly 1986-2020 compared to baseline 1950-1985") +
  theme(plot.title=element_text(size=13, hjust=0.5))
ggsave(anom_plot, file="./era5land_precipanomaly_19812020.png", device="png", units="in", dpi=600, width=8, height=4)



# --------- slope of standardised anomaly (i.e. +ve/-ve trend in anomaly relative to baseline) --------

# extract
p_dat = as.data.frame(anom_i_standardised, xy=TRUE)
p_dat$cellid = 1:nrow(p_dat)
p_dat = p_dat[ !is.na(p_dat$layer.1), ]

p2 = p_dat %>%
  dplyr::select(-x, -y) 
names(p2)[ 1:40 ] = paste("p", 1981:2020, sep="")
p2 = p2 %>%
  reshape2::melt(id.vars=c("cellid")) %>%
  dplyr::mutate(variable = substr(variable, 2, 7))
rm(p_dat)

# p2 %>% 
#   dplyr::filter(cellid == sample(p2$cellid, 1)) %>%
#   ggplot() + 
#   geom_point(aes(as.numeric(variable), value)) + 
#   geom_smooth(aes(as.numeric(variable), value), method="lm")

# calculate slope
slopes = p2 %>%
  dplyr::arrange(cellid, variable) %>%
  dplyr::group_by(cellid) %>%
  dplyr::summarise(
    slope = as.vector(lm(value ~ as.numeric(variable))$coefficients[2])
  )
write.csv(slopes, "./slopes_eraprecanom.csv", row.names=FALSE)

# create raster of slopes
tslopes = raster(p_longtermaverage)
tvals = as.data.frame(tslopes, xy=TRUE)
tvals$cellid = 1:nrow(tvals)
tvals = left_join(tvals, slopes)
tvals = tvals %>% dplyr::arrange(cellid)
values(tslopes) = tvals$slope

writeRaster(tslopes, file="D:/ResearchProjects/202011_fingerprint/environmental_layers/era5_land/precipanom_slope_19812020.tif", format="GTiff")

as.data.frame(tslopes, xy=TRUE) %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=layer)) + 
  scale_fill_gradient2() + 
  coord_fixed()



