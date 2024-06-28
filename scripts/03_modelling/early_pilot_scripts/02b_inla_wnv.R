



# ================ West Nile fever ===================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_setup_functions.R")
source("./scripts/03_modelling/00_covar_extraction_funcs.R")

# set seed for random points
set.seed(20578)


# ------------- read in data and create study area boundary box -------------

# read disease data (n.b. remove )
disease_name = "wnv"
dz = read.csv("C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet.csv") %>%
  dplyr::filter(Disease == "West Nile virus disease" & Presentation == "All cases") %>%
  dplyr::mutate(record_id = 1:length(Disease), 
                ADMcode = as.character(as.integer(ADMcode))) %>%
  dplyr::filter(Year >= 2004)
load("C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet_shp.R")
dz_shp = shp %>%
  dplyr::mutate(ADMcode = as.character(as.integer(ADMcode))) %>%
  dplyr::filter(ADMcode %in% dz$ADMcode) %>%
  st_transform(crs = 4326)

# USA shapefile
study_area <- st_read("./data/shapefiles/cb_2018_us_state_5m.shp") %>% 
  dplyr::filter(!NAME %in% c("Hawaii", "Guam", "American Samoa", "Puerto Rico", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands", "Alaska")) %>%
  st_transform(crs = 4326)

# crop by waf
# bbox <- extent(c(-18, 16, 4, 14))
study_area <- st_union(study_area)
study_area <- st_sf(id = 1, geometry = study_area)

dz %>% 
  ggplot() + 
  geom_sf(data=study_area, fill=NA) + 
  maptheme + 
  geom_point(aes(Longitude, Latitude), size=0.25, color="coral") +
  facet_wrap(~Year)


# ------------ create standardised polygons for entire area --------------

# polygons
polys = dz_shp %>%
  dplyr::select(ADMcode) %>%
  dplyr::full_join(dz[ dz$DataType == "polygon",])
dz = polys
dz$poly_area = as.vector(st_area(dz))/10^6
dz = sf::st_crop(dz, study_area)

# jitter lat-lon for occurrences
locs_jit = data.frame()
for(i in 1:nrow(dz)){
  cat(i)
  locs_jit = rbind(locs_jit, 
                   st_coordinates(sf::st_sample(dz[i, ], size=1)))
}
dz = cbind(dz, locs_jit)
dz$Longitude = dz$X
dz$Latitude = dz$Y


# ---------- create pseudoabsences with same median buffer size -----------

# buffer size for pseudoabsences
med_area_size = median(dz$poly_area)
buf_size = sqrt(med_area_size / pi) * 1000
buf_size 

# num pseudos
npts = 8000
pas <- sp::spsample(as_Spatial(study_area), n=npts, type="random") %>%
  as.data.frame() %>%
  dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
  dplyr::mutate(record_id = paste("pa", 1: npts, sep="_"))
pseudo = st_as_sf(pas, coords = c("Longitude", "Latitude"), crs = 4326)
pbuf = st_buffer(pseudo, buf_size) %>%
  dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude))



# ---------- combine into dataset for modeling ---------

dd = rbind(
  dz %>% dplyr::select(record_id, Longitude, Latitude) %>% dplyr::mutate(presence = 1),
  pbuf %>% dplyr::select(record_id, Longitude, Latitude) %>% dplyr::mutate(presence = 0)
)

plot(study_area$geometry); points(dd$Longitude, dd$Latitude, pch=16, cex=0.25, col="grey50")
points(dd$Longitude[ dd$presence == 1], dd$Latitude[dd$presence==1], pch=16, cex=0.3, col="red")


pl = ggplot() + 
  geom_sf(data=study_area, fill="grey97") + 
  maptheme +
  geom_point(data=as.data.frame(pbuf), aes(Longitude, Latitude), color="coral", size=0.4) +
  geom_point(data=as.data.frame(dz), aes(Longitude, Latitude), color="skyblue4", size=0.5, alpha=0.7)
ggsave(pl, file="./wnv_pts.png", device="png", units="in", width=10, height=7, dpi=600)


 

# ============== covariate processing and extraction ==================

# # biodiversity intactness
#cov1 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/biodiversity_intactness/lbii.asc")
#crs(cov1) = crs(dd)
# cov1 = crop(cov1, study_area)
# names(cov1) = "biodiversity"
#cc = exactextractr::exact_extract(cov1, dd[1, ], fun='mean')

# # poor livestock keepers
# cov2 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/fao_livestockkeepers/rplk.grd")
# names(cov2) = "livestock_keepers"
# cc = exactextractr::exact_extract(cov2, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov2)

# pop dens
cov3 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/popdens_gpw/gpw_v4_population_density_rev11_2010_30_sec.tif")
names(cov3) = "popdens"
cc = exactextractr::exact_extract(cov3, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov3)

# travel time to healthcare (motorised)
cov4 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/healthcare_travel/2020_motorized_travel_time_to_healthcare/2020_motorized_travel_time_to_healthcare.geotiff")
names(cov4) = "healthtraveltime"
cc = exactextractr::exact_extract(cov4, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov4)

# tmean slope
cov5 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/tmean_slope_19812020.tif")
names(cov5) = "tmean_slope"
cc = exactextractr::exact_extract(cov5, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov5)

# tmean anomaly
cov6 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/temperature_meananomaly_19812020.tif")
names(cov6) = "tmean_anomaly"
cc = exactextractr::exact_extract(cov6, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov6)

# agri expansion
cov7 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agriexpansion_masked_20002018.tif")
names(cov7) = "agri_expansion"
cc = exactextractr::exact_extract(cov7, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov7)

# urban expansion
cov8 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanexpansion_masked_20002018.tif")
names(cov8) = "urban_expansion"
cc = exactextractr::exact_extract(cov8, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov8)

# agriculture cover with land mask
cov9 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agri_2018.tif")
mask = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_watermask.tif")
names(cov9) = "agri_cover"; names(mask) = "mask"
cov9 = stack(cov9, mask)

extractcover = function(x){
  x = x[ x$mask == 0, ]
  return( sum(x[, 1] * x$coverage_fraction) / sum(x$coverage_fraction) )
}
cc = exactextractr::exact_extract(terra::rast(cov9), dd)
cc = sapply(cc, extractcover)
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov9)

# urban cover with land mask
cov10 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urban_2018.tif")
mask = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_watermask.tif")
names(cov10) = "urban_cover"; names(mask) = "mask"
cov10 = stack(cov10, mask)
cc = exactextractr::exact_extract(terra::rast(cov10), dd)
cc = sapply(cc, extractcover)
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov10)

# precipitation anomaly
cov11 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_meananomaly_19812020.tif")
names(cov11) = "precip_anomaly"
cc = exactextractr::exact_extract(cov11, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov11)

# precipitation wetness
cov12 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_wetnessanomaly_19812020.tif")
names(cov12) = "precip_wetness"
cc = exactextractr::exact_extract(cov12, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov12)

# log transform extremely overdispersed variables
dd$popdens_log = log(dd$popdens+1)
dd$healthtraveltime_log = log(dd$healthtraveltime+1)
dd = dd %>% 
  dplyr::select(-healthtraveltime, -popdens) %>% 
  st_drop_geometry()

# # view boxplot of covars
plot_covars = dd %>%
  reshape2::melt(id.vars = 1:4) %>%
  dplyr::filter(!variable %in% c("healthtraveltime", "popdens", "y")) %>%
  ggplot() +
  ggforce::geom_sina(aes(x = factor(presence), y=value, group=factor(presence)), width=1, size=0.1, col="grey60") +
  geom_boxplot(aes(x = factor(presence), y=value, group=factor(presence), fill=variable), width=0.25, alpha=0.8) +
  facet_wrap(~variable, scales="free_y") +
  theme_classic() +
  xlab("Presence") +
  ylab("Covariate value") +
  theme(legend.position="none",
        strip.background = element_blank(), strip.text = element_text(size=13))
ggsave(plot_covars, file="./wnv_covars.png", device="png", units="in", dpi=600, width=10, height=10)

# check for multicollinearity - urban cover and urban expansion are reasonably correlated
usdm::vifstep(dd[ , 5:ncol(dd) ])

# scale covariates that are not logged
dd = dd[ complete.cases(dd), ] 
dd[ , c(5:10) ] <- apply(dd[ , c(5:10) ] , 2, scale)

# save 
write.csv(dd, "./output/model_df/wnv_df.csv", row.names=FALSE)
dd = read.csv("./output/model_df/wnv_df.csv")


# function to extract specified covariate
extractCovariateVals = function(sf_data, covariate_name){
  
  # location where all covariates are stored
  predictors_location = "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/"
  
  
  # --- 1. environmental:
  
  # climate: tmean slope
  
  # climate: tmean anomaly
  
  # climate: precip slope
  
  # climate: precip anomaly
  
  # climate: precip wetness
  
  # climate: precip dryness
  
  # deforestation: Hansen via GEE (calls "extractHansen" wrapper function for GEE call)
  
  # biodiversity intactness
  
  # wildfires?
  
  
  # --- 2. human-animal interface
  
  # protected area coverage
  
  # wildlife hunting (defaunation)
  
  # agricultural cover
  
  # agricultural expansion rate
  
  # habitat fragmentation?
  
  # poor livestock keepers
  if(covariate_name == "livestock_keepers"){
    ras = raster::raster(paste(predictors_location, "fao_livestockkeepers/rplk.grd", sep=""))
    names(ras) = "livestock_keepers"
    extracted = exactextractr::exact_extract(ras, sf_data, fun='mean')
  }
  

  # --- 3. social/socio-ecological
  
  # urban cover
  
  # urbanisation rate (either ESA or Liu)
  
  # human population density (log)
  
  # livestock density: cattle
  
  # livestock density: sheep
  
  # livestock density: goats
  
  # livestock density: camels
  
  # travel time to healthcare (log)
  
  # social vulnerability index 
  
}



# ------------------- INLA binomial model (SDM with SPDE ) ---------------------


# --------------- set up INLA objects for spatial field -------------

# specify spatial locations of points for mesh
spatial_locs = dd[ , c("Longitude", "Latitude")]

###  https://haakonbakka.bitbucket.io/btopic104.html
### define max edge, outer bound as relative to extent of our study area
#max.edge = diff(range(spatial_locs$longitude))/20 # 1/15
max.edge = 1.5
bound.outer = diff(range(spatial_locs$Longitude))/7 # outer bound: 1/3 of range diff initially

# create mesh for SPDE with ref to spatial extent of df1
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
plot(mesh1)#; points(spatial_locs, col="red", pch=16)

### define pc matern spde model
# prior median range at approximately 0.5*diff(range(spatial_locs$Longitude)) - roughly extent of study area
# prior prob of median 6 or more is 0.5
# prior probability of marginal sigma of 3 or more is 0.01
spde = inla.spde2.pcmatern(mesh1, 
                           prior.range = c(50, 0.1), 
                           prior.sigma = c(3, 0.01))

# create projector matrix for point locations based on mesh
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

# response variable: cases, with total population as offset
dd$y = dd$presence

# error family
family1 = "binomial"

# priors: 
# iid hyperprior (penalised complexity)
hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))

# fixed effects
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=1, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for fixed effects

# create stack following haakon: https://haakonbakka.bitbucket.io/btopic107.html#3_stationary_models
# define offsets (human rodent contact)
stack1 = inla.stack(data=list(y=dd$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, dd[ , -which(names(dd) %in% c("y"))]),
                                 iidx = 1:nrow(dd)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')




# =============== run inference =================

# function for convenience
# fitINLAmodel = function(formx, verbose=FALSE){
#   return(
#     inla(formx,
#          verbose = verbose,
#          data = inla.stack.data(stack1, spde=spde),
#          family=family1,
#          control.fixed = control.fixed1, 
#          control.predictor=list(A=inla.stack.A(stack1), compute=TRUE, link=1),
#          control.compute=list(cpo=TRUE, waic=TRUE, config=TRUE))
#   )
# }
aa =  fitINLAModel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_anomaly + healthtraveltime_log + precip_anomaly + agri_expansion + urban_expansion"), 
             verbose = TRUE,
             family=family1,
             inla.mode = "experimental")

# fit intercept + spde model ("baseline")
m0 = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde)"), verbose = FALSE)

# fit univariate models with sdpe
m1a = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_slope"))
m1b = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_anomaly"))
m1c = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + popdens_log"))
m1d = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + healthtraveltime_log"))
m1e = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + precip_anomaly"))
m1f = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + precip_wetness"))
m1g = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + agri_expansion"))
m1h = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + urban_expansion"))

# infocrit
exWAIC = function(m) return(m$waic$waic)
ic = data.frame(model=c("baseline",  "tmean_slope", "tmean_anomaly", "popdens_log", "healthtraveltime_log", "precip_anomaly", "precip_wetness", "agri_expansion", "urban_expansion"), 
                waic = sapply(list(m0, m1a, m1b, m1c, m1d, m1e, m1f, m1g, m1h), exWAIC))

# fit multivariate model 
m2 = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_slope + precip_anomaly + popdens_log + healthtraveltime_log + agri_expansion + urban_expansion"))
ic = rbind(ic, data.frame(model="multivariate", waic = exWAIC(m2)))

# fit multivariate nonspatial 
m3 = fitINLAmodel(formx = formula("y ~ -1 + Intercept + tmean_slope + precip_anomaly + popdens_log + healthtraveltime_log + agri_expansion + urban_expansion"))
#ic = rbind(ic, data.frame(model="multivariate_nonspatial", waic = exWAIC(m2)))

waic_plot = ic %>%
  #dplyr::mutate(model = factor(model, levels=c("multivariate", "cattle", "healthtravel", "roads", "forest", "agriculture", "baseline"), ordered=TRUE)) %>%
  ggplot() +
  geom_point(aes(model, waic), size=8, pch=21, fill="skyblue4", alpha=0.8) + 
  theme_bw() + 
  geom_hline(yintercept=ic$waic[ ic$model == "baseline"], lty=2) + 
  coord_flip() +
  ylab("WAIC") + xlab("Model") +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13))

# multivar model without health travel
m4 = fitINLAmodel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_slope + precip_anomaly + popdens_log + agri_expansion + urban_expansion"))
ic = rbind(ic, data.frame(model="multivariate", waic = exWAIC(m2)))



# =============== extract SPDE from multivariate model ============
 
# project SPDE
proj = inla.mesh.projector(mesh1, dims=c(500, 500))
full.proj = expand.grid(x = proj$x, y=proj$y)
full.proj$field = c(inla.mesh.project(proj, m2$summary.random$s$mean))
pras = rasterFromXYZ(full.proj)
mask = fasterize::fasterize(study_area, pras)
pras = mask(pras, mask)
pras = crop(pras, study_area)
pras = as.data.frame(pras, xy=TRUE)

spde_plot = ggplot() + 
  geom_raster(data = pras, aes(x, y, fill=field)) +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black") + 
  geom_point(data=dd[ dd$presence == 1, ], aes(Longitude, Latitude), col="red", size=0.1, alpha=0.3) + 
  theme(legend.position = c(0.05, 0.25), plot.title = element_text(hjust=0.5, size=16), legend.title=element_blank()) +
  ggtitle("West Nile virus disease")
ggsave(spde_plot, file="./wnv_spde.png", device="png", units="in", width=10, height=7, dpi=600)

# ============ extract fixed effects =================

extractFixedINLA = function(model, model_name="mod", transform=FALSE){
  ff = model$summary.fixed
  ff$param = row.names(ff)
  ff$param[ ff$param == "(Intercept)" ] = "Intercept"
  names(ff)[3:5] = c("lower", "median", "upper")
  if(transform == TRUE){
    ff[ 1:5 ] = exp(ff[ 1:5 ])
  }
  ff
}

# from univars (conv to OR)
fu = do.call(
  rbind.data.frame, 
  list(extractFixedINLA(m1a, transform=F), 
       #extractFixedINLA(m1b, transform=F),
       extractFixedINLA(m1c, transform=F), 
       extractFixedINLA(m1d, transform=F),
       extractFixedINLA(m1e, transform=F),
       #extractFixedINLA(m1f, transform=F),
       extractFixedINLA(m1g, transform=F),
       extractFixedINLA(m1h, transform=F)
  )
) %>% 
  dplyr::filter(param != "Intercept") %>%
  dplyr::mutate(type = "Univariate (spatial)")

fm = extractFixedINLA(m4, transform=F) %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::mutate(type = "Multivariate (spatial)")

f_nonspatial =
  extractFixedINLA(m3, transform=F) %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::mutate(type = "Multivariate (non-spatial)")


# viz
params_plot1 = rbind(fu, fm) %>%
  #rbind(f_nonspatial) %>%
  dplyr::mutate(param = replace(param, param == "healthtraveltime_log", "Travel time to\nhealth facility (log)"),
                param = replace(param, param == "livestock_keepers", "Poor livestock\nkeepers"),
                param = replace(param, param == "tmean_slope", "Tmean\nslope"),
                param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
                param = replace(param, param == "agri_expansion", "Agriculture\nexpansion"),
                param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
                param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
                param = replace(param, param == "urban_cover", "Urban\ncover"),
                param = replace(param, param == "agri_cover", "Agriculture\ncover"),
                param = replace(param, param == "popdens_log", "Population\ndensity (log)")) %>%
  #dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", ordered=TRUE))) %>%
  dplyr::smutate(type = factor(type, levels=c("Univariate (spatial)", "Multivariate (spatial)"), ordered=TRUE)) %>%
  ggplot() + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.75) +
  theme_classic() + 
  ylab("Posterior marginal slope (log odds)") + 
  xlab("") + 
  theme(legend.position = "top",
        axis.text = element_text(size=12), 
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13)) +
  scale_color_viridis_d(begin=0, end=0.45)

params_plot2 = rbind(fm, f_nonspatial) %>%
  dplyr::mutate(param = replace(param, param == "healthtraveltime_log", "Travel time to\nhealth facility (log)"),
                param = replace(param, param == "livestock_keepers", "Poor livestock\nkeepers"),
                param = replace(param, param == "tmean_slope", "Tmean\nslope"),
                param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
                param = replace(param, param == "agri_expansion", "Agriculture\nexpansion"),
                param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
                param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
                param = replace(param, param == "urban_cover", "Urban\ncover"),
                param = replace(param, param == "agri_cover", "Agriculture\ncover"),
                param = replace(param, param == "popdens_log", "Population\ndensity (log)")) %>%
  #dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", ordered=TRUE))) %>%
  dplyr::mutate(type = factor(type, levels=c("Multivariate (spatial)", "Multivariate (non-spatial)"), ordered=TRUE)) %>%
  ggplot() + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.75) +
  theme_classic() + 
  ylab("Posterior marginal slope (log odds)") + 
  xlab("") + 
  theme(legend.position = "top",
        axis.text = element_text(size=12), 
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13)) +
  scale_color_viridis_d(begin=0.45, end=0.85)

pp = gridExtra::grid.arrange(spde_plot, params_plot1, params_plot2, ncol=1, heights=c(1, 0.8, 0.8))

ggsave(pp, file="./output/plots/wnv_dashboardtest_jun22.png", device="png", units="in", width=9.5, height=12, dpi=600, scale=0.85)




# ============ visualisation plot ===============

pc1 = gridExtra::grid.arrange(las_spde, plot_covars, nrow=2, heights=c(1, 1))
pc2 = gridExtra::grid.arrange(params_plot, waic_plot, nrow=2, heights=c(1.7, 1))
pc = gridExtra::grid.arrange(pc1, pc2, ncol=2)
ggsave(pc, file="./output/plots/Lassafever_dashboard_2022.png", device="png", units="in", dpi=600, width=18, height=9.5)
