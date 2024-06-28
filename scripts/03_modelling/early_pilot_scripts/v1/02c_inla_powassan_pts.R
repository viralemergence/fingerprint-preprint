



# ================ Powassan virus encephalitis ===================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_setup_functions.R")
source("./scripts/03_modelling/00_covar_extraction_funcs.R")

# set seed for random points
set.seed(209)


# ------------- read in data and create study area boundary box -------------

# read disease data
disease_name = "powassan"
dz = read.csv("C:/Users/roryj/Documents/PhD/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet.csv") %>%
  dplyr::filter(Disease == "Powassan virus disease" & Presentation == "All cases") %>%
  dplyr::mutate(record_id = 1:length(Disease), 
                ADMcode = as.character(as.integer(ADMcode)))
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

# create buffered outline around points to clip study area
buf_sa = st_as_sf(dz, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_union() %>%
  st_convex_hull %>%
  st_buffer(200000) %>%
  rmapshaper::ms_simplify(keep = 0.03) %>%
  smoothr::smooth(method = "chaikin")

study_area = st_intersection(study_area, buf_sa)
study_area <- st_sf(id = 1, geometry = study_area)



# ------------ create standardised occurrence polygons for entire area --------------

# polygons
polys = dz_shp %>%
  dplyr::select(ADMcode) %>%
  dplyr::full_join(dz[ dz$DataType == "polygon",])
dz = polys
dz$poly_area = as.vector(st_area(dz))/10^6
dz = sf::st_crop(dz, study_area)




# ---------- create pseudoabsences/"control" background points with same median buffer size -----------

# buffer size for pseudoabsences
med_area_size = median(dz$poly_area)
buf_size = sqrt(med_area_size / pi) * 1000
buf_size 

# using worldpop 
pd = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_1km_Aggregated.tif") %>%
  crop(study_area)
pd = raster::mask(pd, fasterize::fasterize(study_area, raster(pd), field = "id"))
pd = log(pd + 1)
#plot(pd, col=viridis::cividis(200))
vv = values(pd)
vv[ vv == 0 ] = median(vv[ vv != 0 & !is.na(vv)])
values(pd) = vv

npts = 1000
pas = enmSdm::sampleRast(pd, n=npts, replace=TRUE, prob=TRUE)
if(is.null(pas)){ pas = enmSdm::sampleRast(pd, n=npts, replace=TRUE, prob=TRUE) }
pas = as.data.frame(pas) %>%
  dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
  dplyr::mutate(record_id = paste("pa", 1:npts, sep="_"))

# create buffers
pseudo = st_as_sf(pas, coords = c("Longitude", "Latitude"), crs = 4326)
pbuf = st_buffer(pseudo, buf_size) %>%
  dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
  dplyr::mutate(ADMcode = paste("pa", 1:npts, sep="_"))
#plot(pbuf$geometry)



# ---------- combine into dataset for modeling ---------

dd = rbind(
  dz %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 1),
  pbuf %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 0)
) %>%
  dplyr::filter(!st_is_empty(.))

# viz
plot(study_area$geometry); points(dd$Longitude, dd$Latitude, pch=16, cex=0.25, col="grey50")
points(dd$Longitude[ dd$presence == 1], dd$Latitude[dd$presence==1], pch=16, cex=0.9, col="red")

# # view pop distribution of points
# pdx = pd %>%
#   as.data.frame(xy=TRUE) %>%
#   ggplot() +
#   geom_raster(aes(x, y, fill=layer)) +
#   geom_sf(data=study_area, fill=NA, col="black") +
#   maptheme +
#   scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white", name="Population\n(log)") +
#   geom_point(data = pbuf, aes(Longitude, Latitude), col="darkred", size=0.1, alpha=0.6) +
#   geom_point(data = dz, aes(Longitude, Latitude), col="red", size=0.5, alpha=0.8)
# ggsave(pdx, file="./powassan_weighted_background_wp.png", device="png", units="in", dpi=900, width=9, height=5.5)





# ============== covariate processing and extraction ==================

# # biodiversity intactness
# cov1 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/biodiversity_intactness/lbii.asc")
#crs(cov1) = crs(dd)
# names(cov1) = "biodiversity"
#cc = exactextractr::exact_extract(cov1, dd[1, ], fun='mean')

# # poor livestock keepers
# cov2 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/fao_livestockkeepers/rplk.grd")
# names(cov2) = "livestock_keepers"
# cc = exactextractr::exact_extract(cov2, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov2)

# forest loss (hansen); function calls EarthEngine API, function defined in "covar_extraction_funcs"
c0 = extractHansenGEE(shp=dd)
dd = left_join(dd, c0)

# livestock density (sheep, cattle, goat, buffalo)
c1 = raster::stack(list.files("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/livestock (glw3)/unmodelled_aw/", pattern=".tif", full.names=TRUE))
c1 = c1[[ grep("Ct|Gt|Sh|Bf", names(c1)) ]]
cov1 = crop(sum(c1), study_area)
cov1 = cov1 / raster::area(cov1)
names(cov1) = "livestock_dens"
cc = exactextractr::exact_extract(cov1, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov1)

# # pop dens
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
# cov5 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/tmean_slope_19812020.tif")
# names(cov5) = "tmean_slope"
# cc = exactextractr::exact_extract(cov5, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov5)

# tmean anomaly
cov6 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/temperature_metrics/temperature_meananomaly_19912020.tif")
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

# agri cover
cov9 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_agrimasked_2018.tif")
names(cov9) = "agri_cover"
cc = exactextractr::exact_extract(cov9, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov9)

# urban cover
cov10 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_urbanmasked_2018.tif")
names(cov10) = "urban_cover"
cc = exactextractr::exact_extract(cov10, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov10)

# precipitation anomaly
cov11 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precip_meananomaly_19912020.tif")
names(cov11) = "precip_anomaly"
cc = exactextractr::exact_extract(cov11, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov11)

# # precipitation slope
# cov12 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/era5_land/precip_metrics/precipanom_slope_19812020.tif")
# names(cov12) = "precip_slope"
# cc = exactextractr::exact_extract(cov12, dd, fun='mean')
# dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov12)

# forest cover
cov13 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/esa_cci/esacci_forest_2018.tif")
names(cov13) = "forest_cover"
cc = exactextractr::exact_extract(cov13, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov13)

# social vulnerability (currently US only)
cov14 = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/socialvulnerability_usa/usgrid-us-social-vulnerability-index-tract-2018-wgs84-geotiff/usgrid-us-social-vulnerability-index-tract-2018-wgs84-geotiff/svi_2018_tract_overall_wgs84.tif")
names(cov14) = "social_vulnerability"
cc = exactextractr::exact_extract(cov14, dd, fun='mean')
dd = cbind(dd, cc); names(dd)[ ncol(dd)-1 ] = names(cov14)

# log transform extremely overdispersed variables
dd$popdens_log = log(dd$popdens+1)
dd$healthtraveltime_log = log(dd$healthtraveltime+1)
dd$livestock_log = log(dd$livestock_dens+1)
dd = dd %>% 
  dplyr::select(-healthtraveltime, -popdens, -livestock_dens) %>% 
  st_drop_geometry()

# # view boxplot of covars
plot_covars = dd %>%
  reshape2::melt(id.vars = 1:5) %>%
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
ggsave(plot_covars, file="./powassan_covars.png", device="png", units="in", dpi=600, width=10, height=10)

# save mpx df
write.csv(dd, "./output/model_df/powassan_df_wp.csv", row.names=FALSE)






# ------------------- INLA binomial model (SDM with SPDE ) ---------------------

# inla functions

### fitINLAModel: fit and return INLA model with specified formula and family

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param family likelihood
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, family, verbose=FALSE, config=FALSE, return.marginals=FALSE, inla.mode="experimental"){
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

### rasteriseSPDE: generates a raster of the fitted spatial field masked by actual study area

#' @param fitted_mesh_x = model$summary_random$spde$mean
#' @param mesh_x = mesh object used for model
#' @param shp_x = shapefile of study boundaries used for defining mesh (to mask to correct extent)

rasteriseSPDE = function(fitted_mesh_x, mesh_x, shp_x){
  
  proj = inla.mesh.projector(mesh_x, dims=c(500, 500))
  full.proj = expand.grid(x = proj$x, y=proj$y)
  full.proj$field = c(inla.mesh.project(proj, fitted_mesh_x))
  pras = rasterFromXYZ(full.proj)
  mask = fasterize::fasterize(shp_x, pras)
  pras = mask(pras, mask)
  pras = crop(pras, shp_x)
  return(pras)
}

## extractFixedINLA: extracts fixed effects from specified model

#' @param model fitted INLA object
#' @param model_name modelname to provide in dataframe
#' @param transform exponentiate coefficients to RR/OR scale YN?

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


# ----------- set up for model -------------

# data and keep complete cases (all)
dd = read.csv("./output/model_df/powassan_df_wp.csv") %>%
  dplyr::select(-popdens_log)
sum(complete.cases(dd)) / nrow(dd)
dd = dd[ complete.cases(dd), ]

# scale non log predictors
cols_for_scaling = c("forest_loss", "forest_cover", "tmean_anomaly", "agri_expansion", "agri_cover", "urban_expansion", "urban_cover", "precip_anomaly", "social_vulnerability")
dd[ , which(names(dd) %in% cols_for_scaling) ] <- apply(dd[ , which(names(dd) %in% cols_for_scaling) ], 2, scale)

# # check for multicollinearity - urban cover and urban expansion are reasonably correlated
usdm::vifstep(dd[ , 6:ncol(dd) ])
corrplot::corrplot(cor(dd[ , 6:ncol(dd) ]))

# model priors
hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=1, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for fixed effects

# specify spatial locations of points for mesh
spatial_locs = dd[ , c("Longitude", "Latitude")]

###  https://haakonbakka.bitbucket.io/btopic104.html
### define max edge, outer bound as relative to extent of our study area
#max.edge = diff(range(spatial_locs$longitude))/20 # 1/15
max.edge = 1
bound.outer = diff(range(spatial_locs$Longitude))/7 # outer bound: 1/3 of range diff initially

# create mesh for SPDE with ref to spatial extent of df1
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
plot(mesh1); points(spatial_locs, col="red", pch=16, cex=0.3)

### define pc matern spde model
# prior median range at approximately 0.5*diff(range(spatial_locs$Longitude)) - roughly extent of study area
# prior prob of median 6 or more is 0.5
# prior probability of marginal sigma of 3 or more is 0.01
spde = inla.spde2.pcmatern(mesh1, 
                           prior.range = c(20, 0.1), 
                           prior.sigma = c(5, 0.1))

# create projector matrix for point locations based on mesh
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

# response variable: cases, with total population as offset
dd$y = dd$presence

# error family
family1 = "binomial"

# create stack https://haakonbakka.bitbucket.io/btopic107.html#3_stationary_models
# define offsets (human rodent contact)
stack1 = inla.stack(data=list(y=dd$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, dd[ , -which(names(dd) %in% c("y"))]),
                                 iidx = 1:nrow(dd)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')



# --------------- run inference: baseline model, univariate per predictor, all predictors, causally-informed subset of predictors -------------

# baseline formula and predictors
form_base = "y ~ -1 + Intercept + f(s, model=spde)"
predictors = c("tmean_anomaly", "precip_anomaly", 
               "forest_cover", "forest_loss", 
               "agri_cover", "agri_expansion", 
               "urban_cover", #"urban_expansion",
               "healthtraveltime_log", "livestock_log",
               "social_vulnerability")

# baseline SPDE model
m0 = fitINLAModel(formx = formula(form_base), family=family1)

# univariate models (store in list)
m1 = vector("list", length=length(predictors))

for(i in 1:length(predictors)){
  
  pred = predictors[i]
  print(paste("Fitting", pred, sep=" "))
  form_i = formula(paste(form_base, pred, sep = " + "))
  m1_i = fitINLAModel(formx = form_i, family=family1)
  m1[[i]] = m1_i
  
}

# multivariate model (all preds)
form = formula(paste(c(form_base, predictors[ -which(predictors %in% c("healthtraveltime_log", "urban_expansion", "agri_cover")) ]), collapse=" + "))
m2 = fitINLAModel(formx=form, family=family1)

# multivariate model (causal subset)
preds_c = c("tmean_anomaly", "precip_anomaly", "agri_expansion",
            "forest_cover", "forest_loss", "urban_cover", "social_vulnerability")
form = formula(paste(c(form_base, preds_c), collapse=" + "))
m3 = fitINLAModel(formx=form, family=family1)




# ------------ view model results --------

# extract and view spatial field SPDE from full model
spde_plot = rasteriseSPDE(m2$summary.random$s$mean, mesh1, study_area) %>%
  as.data.frame(xy=TRUE) %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=field)) +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black") + 
  geom_point(data=dd[ dd$presence == 1, ], aes(Longitude, Latitude), col="red", size=0.8) + 
  theme(legend.position = c(0.98, 0.15), plot.title = element_text(hjust=0.5, size=16)) +
  ggtitle("Powassan virus encephalitis")


# extract fitted parameters
# univariate models
f1 = do.call(
  rbind.data.frame, 
  lapply(m1, function(x) extractFixedINLA(x))
) %>%
  dplyr::mutate(type = "Univariate")

f2 = extractFixedINLA(m2) %>%
  dplyr::mutate(type = "Multivariate")

f3 = extractFixedINLA(m3) %>%
  dplyr::mutate(type = "Causal")

# combine and save
fx = do.call(rbind.data.frame, list(f1, f2, f3)) %>% 
  dplyr::mutate(Disease = "Powassan encephalitis",
                Num_observations = sum(dd$presence))
row.names(fx) = c()
write.csv(fx, "./output/model_outputs/disease_params/powassan_params.csv", row.names=FALSE)


# plot parameters from 3 models
fx %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::mutate(param = replace(param, param == "healthtraveltime_log", "Travel time to\nhealth facility (log)"),
                param = replace(param, param == "livestock_log", "Livestock\ndensity (log)"),
                param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
                param = replace(param, param == "forest_cover", "Forest\ncover"),
                param = replace(param, param == "forest_loss", "Forest\nloss"),
                param = replace(param, param == "social_vulnerability", "Social\nvulnerability"),
                param = replace(param, param == "agri_expansion", "Agriculture\nexpansion"),
                param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
                param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
                param = replace(param, param == "urban_cover", "Urban\ncover"),
                param = replace(param, param == "agri_cover", "Agriculture\ncover"),
                param = replace(param, param == "popdens_log", "Population\ndensity (log)")) %>%
  #dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", ordered=TRUE))) %>%
  dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", "Causal"), ordered=TRUE)) %>%
  ggplot() + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.6, alpha=0.7) +
  theme_classic() + 
  ylab("Posterior marginal slope (log odds)") + 
  xlab("") + 
  theme(legend.position = "top",
        axis.text = element_text(size=13), 
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13)) +
  scale_color_viridis_d(begin=0, end=0.8, direction=-1)

# uni + multi for colin's talk
params_plot = fx %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::filter(type != "Causal") %>%
  dplyr::mutate(param = replace(param, param == "healthtraveltime_log", "Travel time to\nhealth facility (log)"),
                param = replace(param, param == "livestock_log", "Livestock\ndensity (log)"),
                param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
                param = replace(param, param == "forest_cover", "Forest\ncover"),
                param = replace(param, param == "forest_loss", "Forest\nloss"),
                param = replace(param, param == "agri_expansion", "Agriculture\nexpansion"),
                param = replace(param, param == "social_vulnerability", "Social\nvulnerability"),
                param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
                param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
                param = replace(param, param == "urban_cover", "Urban\ncover"),
                param = replace(param, param == "agri_cover", "Agriculture\ncover"),
                param = replace(param, param == "popdens_log", "Population\ndensity (log)")) %>%
  dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", ordered=TRUE))) %>%
  #dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", "Causal"), ordered=TRUE)) %>%
  ggplot() + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.6, alpha=0.7) +
  theme_classic() + 
  ylab("Posterior marginal slope (log odds)") + 
  xlab("") + 
  theme(legend.position = "top",
        axis.text = element_text(size=13), 
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13)) +
  scale_color_viridis_d(begin=0, end=0.8, direction=-1)

pp = gridExtra::grid.arrange(spde_plot, params_plot, ncol=1, heights=c(1, 0.8))
ggsave(pp, file="./output/plots/powassan_dashboardtest_jun22.png", device="png", units="in", width=13.5, height=11, dpi=600, scale=0.9)
