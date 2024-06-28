



# ================ West Nile virus data degradation exercise ===================

# Test data degradation approach for WNV from most to least data rich
# 1. model annual case counts in areal model (BYM2, nbinomial)
#


# Test approach for WNV: examine 
# 1. differences between "outbreak" and num cases definitions visually
# 2. model annual "outbreaks" with polygon absences in IHME way (binomial)
# 3. model annual case counts with polygon absences in IHME way (neg binom/poisson)
# 4. model "outbreaks" with point pseudoabsences (binomial)

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_setup_functions.R")

# set seed for random points
set.seed(7000)


# ------------- read in data and create study area boundary box -------------

# read disease data
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

# USA shapefile for study area
study_area <- st_read("./data/shapefiles/cb_2018_us_state_5m.shp") %>% 
  dplyr::filter(!NAME %in% c("Hawaii", "Guam", "American Samoa", "Puerto Rico", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands", "Alaska")) %>%
  st_transform(crs = 4326)
study_area <- st_union(study_area)
study_area <- st_sf(id = 1, geometry = study_area)

# offshore_areas
offshore_areas = st_read("./data/shapefiles/cb_2018_us_state_5m.shp") %>% 
  dplyr::filter(NAME %in% c("Hawaii", "Guam", "American Samoa", "Puerto Rico", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands", "Alaska"))



# --------- 1. generate dataset (occurrence polygons using US shapefile, and pseudoabsences) -----------

# a. polygons
polys = dz_shp %>%
  dplyr::select(ADMcode) %>%
  dplyr::full_join(dz[ dz$DataType == "polygon",])
dz = polys
dz$poly_area = as.vector(st_area(dz))/10^6
dz = sf::st_crop(dz, study_area)

# visualisation: cases and annual outbreaks
cases_summarised = dz %>%
  st_drop_geometry() %>%
  dplyr::group_by(ADMcode) %>%
  dplyr::summarise(NumOutbreaks = length(NumCases),
                   TotalCases = sum(NumCases),
                   MeanCases = mean(NumCases))
cases_summarised = dz %>%
  dplyr::filter(!duplicated(ADMcode)) %>%
  dplyr::select(-NumCases) %>%
  dplyr::left_join(cases_summarised)

p1 = cases_summarised %>%
  ggplot() +
  geom_sf(aes(fill=NumOutbreaks), col=NA) +
  geom_sf(data=study_area, fill=NA, col="black") +
  maptheme +
  scale_fill_viridis_c(option="magma", name="Num years\nwith outbreaks", direction=-1) +
  ggtitle("Number of outbreak years") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5))

p2 = cases_summarised %>%
  ggplot() +
  geom_sf(aes(fill=log(MeanCases)), col=NA) +
  geom_sf(data=study_area, fill=NA, col="black") +
  maptheme +
  scale_fill_viridis_c(option="magma", direction=-1) +
  ggtitle("Mean annual cases (log)") +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        plot.title = element_text(size=20, hjust=0.5))

pc = gridExtra::grid.arrange(p2, p1, ncol=1)
ggsave(pc, file="./output/wnv_modeltest/00_westnile_maps.png", device="png", units="in", width=8, height=8, dpi=600)




# --------- 2. extract environ data for all USA polygons and for point buffers ----------

# --- a.polygons 

# dd = sf::st_read("C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/arbonet_countyarboviral/cb_2018_us_county_500k.shp") %>%
#   dplyr::mutate(CTYCODE = as.integer(COUNTYFP),
#                 STCODE = as.integer(STATEFP),
#                 ADMcode = paste(STATEFP, COUNTYFP, sep="")) %>% 
#   st_transform(crs = 4326) 
# 
# # latlon centroids
# centroids = st_coordinates(st_centroid(dd)) %>%
#   as.data.frame() %>%
#   dplyr::rename("Longitude"=1, "Latitude"=2) %>%
#   dplyr::mutate(ADMcode = dd$ADMcode)
# 
# # # biodiversity intactness
# # cov1 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/biodiversity_intactness/lbii.asc")
# #crs(cov1) = crs(dd)
# # names(cov1) = "biodiversity"
# #cc = exactextractr::exact_extract(cov1, dd[1, ], fun='mean')
# 
# # # poor livestock keepers
# # cov2 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/fao_livestockkeepers/rplk.grd")
# # names(cov2) = "livestock_keepers"
# # cc = exactextractr::exact_extract(cov2, dd, fun='mean')
# # dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov2)
# 
# # pop dens
# cov3 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/popdens_gpw/gpw_v4_population_density_rev11_2010_30_sec.tif")
# names(cov3) = "popdens"
# cc = exactextractr::exact_extract(cov3, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov3)
# 
# # travel time to healthcare (motorised)
# cov4 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/healthcare_travel/2020_motorized_travel_time_to_healthcare/2020_motorized_travel_time_to_healthcare.geotiff")
# names(cov4) = "healthtraveltime"
# cc = exactextractr::exact_extract(cov4, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov4)
# 
# # tmean slope
# cov5 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/tmean_slope_19812020.tif")
# names(cov5) = "tmean_slope"
# cc = exactextractr::exact_extract(cov5, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov5)
# 
# # tmean anomaly
# cov6 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/temperature_meananomaly_19812020.tif")
# names(cov6) = "tmean_anomaly"
# cc = exactextractr::exact_extract(cov6, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov6)
# 
# # agri expansion
# cov7 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agriexpansion_masked_20002018.tif")
# names(cov7) = "agri_expansion"
# cc = exactextractr::exact_extract(cov7, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov7)
# 
# # urban expansion
# cov8 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanexpansion_masked_20002018.tif")
# names(cov8) = "urban_expansion"
# cc = exactextractr::exact_extract(cov8, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov8)
# 
# # agri cover
# cov9 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agrimasked_2018.tif")
# names(cov9) = "agri_cover"
# cc = exactextractr::exact_extract(cov9, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov9)
# 
# # urban cover
# cov10 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanmasked_2018.tif")
# names(cov10) = "urban_cover"
# cc = exactextractr::exact_extract(cov10, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov10)
# 
# # precipitation anomaly
# cov11 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_meananomaly_19812020.tif")
# names(cov11) = "precip_anomaly"
# cc = exactextractr::exact_extract(cov11, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov11)
# 
# # precipitation wetness
# cov12 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_wetnessanomaly_19812020.tif")
# names(cov12) = "precip_wetness"
# cc = exactextractr::exact_extract(cov12, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov12)
# 
# # population count
# cov12 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/popdens_gpw/gpw_v4_population_count_rev11_2010_30_sec.tif")
# names(cov12) = "population"
# cc = exactextractr::exact_extract(cov12, dd, fun='sum')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov12)
# 
# # log transform extremely overdispersed variables
# dd$popdens_log = log(dd$popdens+1)
# dd$healthtraveltime_log = log(dd$healthtraveltime+1)
# dd = dd %>% 
#   dplyr::select(-healthtraveltime, -popdens) %>% 
#   st_drop_geometry()
# 
# # log pop
# dd$logpop = log(dd$population+1)
# 
# # centroids
# dd = left_join(dd, centroids)
# usdm::vifstep(dd[ , 13:22])
# 
# # save
# write.csv(dd, "./output/wnv_modeltest/usa_covars.csv", row.names=FALSE)




# ================= Four modelling approaches ====================

# Gradually eroding from the richest to the most patchy data
# 1. Annual case counts using all years and polygons (random intercept for year)
# 2. Annual presence/absence of outbreaks using all years and polygons (random intercept for year)
# 3. Annual presence/absence of outbreaks using "pseudoabsence" random sample of polygons (spatial only, no temporal metadata)
# 4. Annual presence/absence of outbreaks using point pseudoabsences with median buffer size (spatial only)

# INLA funcs
# priors
hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=1, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for fixed effects

### fitINLAModel: fit and return INLA model with specified formula and family

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param family likelihood
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, family, verbose=FALSE, config=FALSE, return.marginals=FALSE, inla.mode="classic"){
  return(
    inla(formx,
         verbose = verbose,
         data = inla.stack.data(stack1, spde=spde),
         family=family,
         control.fixed = control.fixed1, 
         control.predictor=list(A=inla.stack.A(stack1), 
                                compute=TRUE, 
                                link=1),
         control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, 
                              config=config, 
                              return.marginals=return.marginals),
         control.inla = list(strategy='adaptive', # adaptive gaussian
                             cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
         inla.mode = inla.mode)
  )
}


# --------- Build dataframe for model -------------

# read in polygon data and exclude offshores
dd = read.csv("./output/wnv_modeltest/usa_covars.csv") %>%
  dplyr::filter(! STATEFP %in% as.numeric(offshore_areas$STATEFP))

# all years
db = expand.grid(unique(dd$ADMcode), 2004:2020) %>%
  as.data.frame() %>%
  dplyr::rename("ADMcode"=1, "Year"=2) %>%
  dplyr::mutate(ADMcode = as.integer(as.character(ADMcode))) %>%
  dplyr::left_join(
    dz %>% dplyr::select(ADMcode, Year, NumCases) %>% dplyr::mutate(ADMcode = as.integer(as.character(ADMcode))) %>% st_drop_geometry()
  ) %>%
  dplyr::mutate(NumCases = replace(NumCases, is.na(NumCases), 0),
                Presence = ifelse(NumCases > 0, 1, 0))
if(sum(db$NumCases) == sum(dz$NumCases)) print("Test passed: case counts match")

# combine with environmental data
db = left_join(db, dd[ , 12:ncol(dd)] %>% dplyr::mutate(ADMcode = as.integer(as.character(ADMcode))))

# scale predictors that are not logged
db[ , c(5:12)] = apply(db[ , c(5:12) ] , 2, scale)



# --------- 1. Model 1: case counts ('full' model) --------------

# specify spatial locations of points for mesh
spatial_locs = db[ , c("Longitude", "Latitude")]

# create mesh for SPDE with ref to spatial extent of df1
max.edge = 1.2
bound.outer = diff(range(spatial_locs$Longitude))/7 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
plot(mesh1)#; points(spatial_locs, col="red", pch=16)

### define pc matern spde model
# prior median range at approximately 0.5*diff(range(spatial_locs$Longitude)) - roughly extent of study area
# prior probability of marginal sigma of 3 or more is 0.01
# create projector matrix for point locations based on mesh
# alpha 3/2 https://groups.google.com/g/r-inla-discussion-group/c/ZhZVu8YPI8I
spde = inla.spde2.pcmatern(mesh1, 
                           alpha=3/2,
                           prior.range = c(50, 0.1), 
                           prior.sigma = c(3, 0.01))
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

# response variable: cases, with total population as offset
db$y = db$NumCases
family1 = "nbinomial"

# create stack following haakon: https://haakonbakka.bitbucket.io/btopic107.html#3_stationary_models
stack1 = inla.stack(data=list(y=db$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, db[ , -which(names(db) %in% c("y"))]),
                                 iidx = 1:nrow(db)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')

# fit intercept + spde model with year random effect
m1 = fitINLAModel(formx = formula("y ~ -1 + Intercept + offset(logpop) + f(s, model=spde) + f(Year, model='iid') + tmean_anomaly + healthtraveltime_log + precip_anomaly + agri_cover + agri_expansion + urban_cover"),
                  family=family1, 
                  verbose = TRUE,
                  inla.mode="experimental")
f1 = extractFixedINLA(m1, transform=TRUE) %>% dplyr::mutate(type = "Annual case counts (all data with year random effect)")

# yearly effect
# extractRandomINLA(m1$summary.random$Year, effect_name = "Year", transform=FALSE) %>%
#   ggplot() + 
#   geom_point(aes(value, median)) + geom_linerange(aes(value, ymin=lower, ymax=upper), alpha=0.4) + 
#   geom_hline(yintercept=0)

# compare to "classic" INLA mode - inference is the same
# m1b = fitINLAModel(formx = formula("y ~ -1 + Intercept + offset(logpop) + f(s, model=spde) + f(Year, model='iid') + tmean_anomaly + healthtraveltime_log + precip_anomaly + agri_cover + agri_expansion + urban_cover"),
#                    family=family1, 
#                    verbose = TRUE, 
#                    inla.mode="classic")
# f1b = extractFixedINLA(m1b, transform=TRUE) %>% dplyr::mutate(type = "Annual case counts (all data, classic)") 
# 
# f1 %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::mutate(param = replace(param, param == "healthtraveltime_log", "Travel time to\nhealth facility (log)"),
#                 param = replace(param, param == "livestock_keepers", "Poor livestock\nkeepers"),
#                 param = replace(param, param == "tmean_slope", "Tmean\nslope"),
#                 param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
#                 param = replace(param, param == "agri_expansion", "Agriculture\nexpansion"),
#                 param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
#                 param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
#                 param = replace(param, param == "urban_cover", "Urban\ncover"),
#                 param = replace(param, param == "agri_cover", "Agriculture\ncover"),
#                 param = replace(param, param == "popdens_log", "Population\ndensity (log)")) %>%
#   ggplot() +
#   geom_hline(yintercept=1, lty=2) +
#   geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) +
#   geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.75) +
#   theme_classic() +
#   ylab("Posterior marginal slope") +
#   xlab("") +
#   theme(legend.position = "bottom",
#         axis.text = element_text(size=12),
#         axis.title = element_text(size=13),
#         legend.title = element_blank(),
#         legend.text = element_text(size=13)) +
#   scale_color_viridis_d(begin=0, end=0.85) +
#   guides(color = guide_legend(ncol=1))




# ----------- 2. Model 2: presence-absence annual "outbreaks" -----------------

# response variable: cases, with total population as offset
db$y = db$Presence
family1 = "binomial"

# create stack following haakon
stack1 = inla.stack(data=list(y=db$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, db[ , -which(names(db) %in% c("y"))]),
                                 iidx = 1:nrow(db)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')

# model and extract fixed effects
m2 = fitINLAModel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + f(Year, model='iid') + tmean_anomaly + healthtraveltime_log + precip_anomaly + agri_cover + agri_expansion + urban_cover"), 
                  verbose = TRUE,
                  family=family1,
                  inla.mode = "experimental")
f2 = extractFixedINLA(m2, transform=TRUE) %>%
  dplyr::mutate(type = "Annual case occurrence (all data with year random effect)") 



# ----------- 3. Model 3: spatial outbreaks with polygon pseudoabsences (no temporal) ------------------

# create presence and pseudoabsence data (10,000)
npts = 10000
pres = db %>% dplyr::filter(Presence == 1)
abs = db[ sample(1:nrow(db), size=npts, replace=TRUE), ] %>% dplyr::mutate(Presence = 0)
db = rbind(pres, abs)

# response variable: cases, with total population as offset
db$y = db$Presence
family1 = "binomial"

#  spatial objects for INLA model
spatial_locs = db[ , c("Longitude", "Latitude")]
max.edge = 1.2
bound.outer = diff(range(spatial_locs$Longitude))/7 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
spde = inla.spde2.pcmatern(mesh1, 
                           alpha=3/2,
                           prior.range = c(50, 0.1), 
                           prior.sigma = c(3, 0.01))
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

#create stack
stack1 = inla.stack(data=list(y=db$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, db[ , -which(names(db) %in% c("y"))]),
                                 iidx = 1:nrow(db)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')

# model and extract fixed effects
m3 = fitINLAModel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_anomaly + healthtraveltime_log + precip_anomaly + agri_cover + agri_expansion + urban_cover"), 
                  verbose = TRUE,
                  family=family1,
                  inla.mode = "experimental")
f3 = extractFixedINLA(m3, transform=TRUE) %>%
  dplyr::mutate(type = "Spatial 'outbreak' occurrence (polygon presence + polygon pseudoabsences)") 






# =============== model 4: annual outbreak points with point pseudoabsences ================

# read in and extract covariates via point pipeline

# ------------ disease data --------------

# disease_name = "wnv"
# dz = read.csv("C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet.csv") %>%
#   dplyr::filter(Disease == "West Nile virus disease" & Presentation == "All cases") %>%
#   dplyr::mutate(record_id = 1:length(Disease), 
#                 ADMcode = as.character(as.integer(ADMcode))) %>%
#   dplyr::filter(Year >= 2004)
# load("C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet_shp.R")
# dz_shp = shp %>%
#   dplyr::mutate(ADMcode = as.character(as.integer(ADMcode))) %>%
#   dplyr::filter(ADMcode %in% dz$ADMcode) %>%
#   st_transform(crs = 4326)
# 
# 
# # ------------ create standardised polygons for entire area --------------
# 
# # polygons
# polys = dz_shp %>%
#   dplyr::select(ADMcode) %>%
#   dplyr::full_join(dz[ dz$DataType == "polygon",])
# dz = polys
# dz$poly_area = as.vector(st_area(dz))/10^6
# dz = sf::st_crop(dz, study_area)
# 
# # jitter lat-lon for occurrences
# locs_jit = data.frame()
# for(i in 1:nrow(dz)){
#   cat(i)
#   locs_jit = rbind(locs_jit, st_coordinates(sf::st_sample(dz[i, ], size=1)))
# }
# dz = cbind(dz, locs_jit)
# dz$Longitude_jit = dz$X
# dz$Latitude_jit = dz$Y
# 
# 
# # ---------- create pseudoabsences with same median buffer size -----------
# 
# # buffer size for pseudoabsences
# med_area_size = median(dz$poly_area)
# buf_size = sqrt(med_area_size / pi) * 1000
# buf_size 
# 
# # num pseudos
# npts = 10000
# pas <- sp::spsample(as_Spatial(study_area), n=npts, type="random") %>%
#   as.data.frame() %>%
#   dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
#   dplyr::mutate(record_id = paste("pa", 1: npts, sep="_"))
# pseudo = st_as_sf(pas, coords = c("Longitude", "Latitude"), crs = 4326)
# pbuf = st_buffer(pseudo, buf_size) %>%
#   dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
#   dplyr::mutate(Longitude_jit = Longitude, Latitude_jit = Latitude)
# 
# 
# 
# # ---------- combine into dataset for modeling ---------
# 
# dd = rbind(
#   dz %>% dplyr::select(record_id, Longitude, Latitude, Longitude_jit, Latitude_jit) %>% dplyr::mutate(Presence = 1),
#   pbuf %>% dplyr::select(record_id, Longitude, Latitude, Longitude_jit, Latitude_jit) %>% dplyr::mutate(Presence = 0)
# )



# -------- extract covars -------------

# # pop dens
# cov3 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/popdens_gpw/gpw_v4_population_density_rev11_2010_30_sec.tif")
# names(cov3) = "popdens"
# cc = exactextractr::exact_extract(cov3, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov3)
# 
# # travel time to healthcare (motorised)
# cov4 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/healthcare_travel/2020_motorized_travel_time_to_healthcare/2020_motorized_travel_time_to_healthcare.geotiff")
# names(cov4) = "healthtraveltime"
# cc = exactextractr::exact_extract(cov4, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov4)
# 
# # tmean slope
# cov5 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/tmean_slope_19812020.tif")
# names(cov5) = "tmean_slope"
# cc = exactextractr::exact_extract(cov5, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov5)
# 
# # tmean anomaly
# cov6 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/temperature_meananomaly_19812020.tif")
# names(cov6) = "tmean_anomaly"
# cc = exactextractr::exact_extract(cov6, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov6)
# 
# # agri expansion
# cov7 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agriexpansion_masked_20002018.tif")
# names(cov7) = "agri_expansion"
# cc = exactextractr::exact_extract(cov7, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov7)
# 
# # urban expansion
# cov8 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanexpansion_masked_20002018.tif")
# names(cov8) = "urban_expansion"
# cc = exactextractr::exact_extract(cov8, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov8)
# 
# # agri cover
# cov9 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agrimasked_2018.tif")
# names(cov9) = "agri_cover"
# cc = exactextractr::exact_extract(cov9, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov9)
# 
# # urban cover
# cov10 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanmasked_2018.tif")
# names(cov10) = "urban_cover"
# cc = exactextractr::exact_extract(cov10, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov10)
# 
# # precipitation anomaly
# cov11 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_meananomaly_19812020.tif")
# names(cov11) = "precip_anomaly"
# cc = exactextractr::exact_extract(cov11, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov11)
# 
# # precipitation wetness
# cov12 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_wetnessanomaly_19812020.tif")
# names(cov12) = "precip_wetness"
# cc = exactextractr::exact_extract(cov12, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov12)
# 
# # population count
# cov12 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/popdens_gpw/gpw_v4_population_count_rev11_2010_30_sec.tif")
# names(cov12) = "population"
# cc = exactextractr::exact_extract(cov12, dd, fun='sum')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov12)
# 
# # log transform extremely overdispersed variables
# dd$popdens_log = log(dd$popdens+1)
# dd$healthtraveltime_log = log(dd$healthtraveltime+1)
# dd = dd %>%
#   dplyr::select(-healthtraveltime, -popdens) 
# 
# # log pop
# dd$logpop = log(dd$population+1)

# save
#save(dd, file="C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/output/wnv_modeltest/wnv_data_pseudopoints.R")
load(file="C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/output/wnv_modeltest/wnv_data_pseudopoints.R")

# drop geometry for modelling
dd = dd %>% st_drop_geometry()

# scale predictors that are not logged
dd = dd[ complete.cases(dd), ] 
dd[ , c(7:14) ] <- apply(dd[ , c(7:14) ] , 2, scale)



# -------- fit model ------

# response variable: cases, with total population as offset
db = dd
db$y = db$Presence
family1 = "binomial"

#  spatial objects for INLA model
spatial_locs = db[ , c("Longitude_jit", "Latitude_jit")]
max.edge = 1.2
bound.outer = diff(range(spatial_locs$Longitude))/7 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
spde = inla.spde2.pcmatern(mesh1, 
                           alpha=3/2,
                           prior.range = c(50, 0.1), 
                           prior.sigma = c(3, 0.01))
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

#create stack
stack1 = inla.stack(data=list(y=db$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, db[ , -which(names(db) %in% c("y"))]),
                                 iidx = 1:nrow(db)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')

# model and extract fixed effects (n.b. no year)
m4 = fitINLAModel(formx = formula("y ~ -1 + Intercept + f(s, model=spde) + tmean_anomaly + healthtraveltime_log + precip_anomaly + agri_cover + agri_expansion + urban_cover"), 
                  verbose = TRUE,
                  family=family1,
                  inla.mode = "experimental")
f4 = extractFixedINLA(m4, transform=TRUE) %>%
  dplyr::mutate(type = "Spatial 'outbreak' occurrence (polygon presence + point-buffer pseudoabsences)") 





# ----------- combine all and view -------------

estimates = rbind(f1, f2) %>% rbind(f3) %>% rbind(f4)
row.names(estimates) = c()
write.csv(estimates, "./output/wnv_modeltest/polygonmodels_parameters.csv", row.names=FALSE)


plotx = estimates %>%
  dplyr::mutate(type = factor(type, levels=unique(estimates$type), ordered=TRUE)) %>%
  dplyr::filter(param != "Intercept") %>%
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
  ggplot() +
  geom_hline(yintercept=1, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=3) +
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.75) +
  theme_classic() +
  ylab("Posterior marginal slope") +
  xlab("") +
  theme(legend.position = "bottom",
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13)) +
  scale_color_viridis_d(begin=0, end=0.85) +
  guides(color = guide_legend(ncol=1))
ggsave(plotx, file="./output/wnv_modeltest/00_fixedeffects_modelcomparison.png", device="png", units="in", width=9, height=6, dpi=600)






# --------- extract and view fitted SPDE for spatial outbreak model -----------

# project SPDE
proj = inla.mesh.projector(mesh1, dims=c(500, 500))
full.proj = expand.grid(x = proj$x, y=proj$y)
full.proj$field = c(inla.mesh.project(proj, m4$summary.random$s$mean))
pras = rasterFromXYZ(full.proj)
mask = fasterize::fasterize(study_area, pras)
pras = mask(pras, mask)
pras = crop(pras, study_area)
pras = as.data.frame(pras, xy=TRUE)
spde_plot = ggplot() +
  geom_raster(data = pras, aes(x, y, fill=field)) +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white", name="Random\nintercept\n(spde)") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black") +
  geom_point(data=dd[ dd$Presence == 1, ], aes(Longitude_jit, Latitude_jit), col="red", size=0.1, alpha=0.2) +
  theme(legend.position = c(0.95, 0.25), plot.title = element_text(hjust=0.5, size=16)) +
  ggtitle("West Nile virus disease")
ggsave(spde_plot, file="./output/wnv_modeltest/00_spde_outbreakmodel.png", device="png", units="in", width=9, height=6, dpi=600)

