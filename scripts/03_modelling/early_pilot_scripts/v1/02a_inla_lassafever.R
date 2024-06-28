



# ================ Lassa fever ===================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_setup_functions.R")
source("./scripts/03_modelling/00_covar_extraction_funcs.R")

# natural earth coastlines
study_area = rnaturalearth::ne_countries(scale=110, returnclass="sf")

# set seed for random points
set.seed(106)


# ------------- read in data and create study area boundary box -------------

# read disease data
disease_name = "lassa"
dz = read.csv("./output/spillovers_processed/spillovers_lassa.csv") %>%
  dplyr::filter(Year >= 1985) %>%
  distinct() %>%
  dplyr::mutate(record_id = 1:length(Disease))
load("./output/spillovers_processed/spillovers_lassa_shp.R")
dz_shp = shp

# West Africa shapefile
study_area <- st_read("./data/shapefiles/africa_shp.shp") %>% dplyr::filter(NAME != "Sao Tome and Principe")
bbox <- extent(c(-18, 18, 0, 20))
study_area <- st_crop(st_union(study_area), bbox)

# create buffered outline around points to clip study area
buf_sa = st_as_sf(dz, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_union() %>%
  st_convex_hull %>%
  st_buffer(150000) %>%
  rmapshaper::ms_simplify(keep = 0.03) %>%
  smoothr::smooth(method = "chaikin")
study_area = st_intersection(study_area, buf_sa)
study_area <- st_sf(id = 1, geometry = study_area)
plot(study_area$geometry)



# ------------ create standardised polygons for entire area --------------

# create buffer for points
pts = dz %>% dplyr::filter(DataType == "point")
pts$ADMcode = paste("Pts", 1:nrow(pts), sep="_")
pts_sf = st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326)
buf_size = 5000 # 5km buffer for points
buf = st_buffer(pts_sf, buf_size) %>%
  dplyr::left_join(pts %>% dplyr::select(record_id, Longitude, Latitude))

# polygons with buffer
polys = dz_shp %>%
  dplyr::select(ADMcode) %>%
  dplyr::full_join(dz[ dz$DataType == "polygon",])

# combine
dz = rbind(buf, polys)
dz$poly_area = as.vector(sf::st_area(dz)/10^6)

# filtering out all records with poly_area of > threshold (5000km, approx 70km2 cell size if grid)
area_thresh = 5000
dz = dz %>% dplyr::filter(poly_area <= area_thresh)

# ggplot() +
#   geom_sf(data=study_area, fill=NA, col="black") +
#   maptheme +
#   geom_sf(data=dz, col="red", pch=16)


# ---------- create pseudoabsences with same median buffer size -----------

# buffer size for pseudoabsences
med_area_size = median(dz$poly_area)
buf_size = sqrt(med_area_size / pi) * 1000
buf_size 

# # random pseudos
# npts = 3000
# pas <- sp::spsample(as_Spatial(study_area), n=npts, type="random") %>%
#   as.data.frame() %>%
#   dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
#   dplyr::mutate(record_id = paste("pa", 1: npts, sep="_"))
# pseudo = st_as_sf(pas, coords = c("Longitude", "Latitude"), crs = 4326)
# pbuf = st_buffer(pseudo, buf_size) %>%
#   dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude))

# using worldpop 
pd = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_1km_Aggregated.tif") %>%
  crop(study_area)
pd = raster::mask(pd, fasterize::fasterize(study_area, raster(pd), field = "id"))
pd = log(pd + 1)
vv = values(pd)
vv[ vv == 0 ] = median(vv[ vv != 0 & !is.na(vv)])
values(pd) = vv
#plot(pd, col=viridis::cividis(200))

npts = 2500
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




# ---------- combine into dataset for modeling ---------

dd = rbind(
  dz %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 1),
  pbuf %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 0)
)

# plot(study_area$geometry); points(dd$Longitude, dd$Latitude, pch=16, cex=0.25, col="grey50")
# points(dd$Longitude[ dd$presence == 1], dd$Latitude[dd$presence==1], pch=16, cex=0.9, col="red")

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
# ggsave(pdx, file="./lassa_weighted_background_wp.png", device="png", units="in", dpi=900, width=9, height=5.5)




# ============== covariate processing and extraction ==================

# full covariates df
covs_df = data.frame(
  cov = c("tmean_anomaly", "precip_anomaly", "forest_cover", "forest_loss", "crop_cover", "crop_expansion",
          "urban_cover", "urban_expansion", "evi_coefvar", "evi_dissimilarity", "livestock_all", "health_travel"),
  description = c("Mean standardised annual mean temperature anomaly 1991-2020 (compared to baseline 1950-1990) from ERA5-Land",
                  "Mean standardised annual precipitation anomaly 1991-2020 (compared to baseline 1950-1990) from ERA5-Land",
                  "Proportion forest cover 2015 from Coperncius 100m fractional land cover",
                  "Proportion area experiencing tree cover loss between 2000 and 2019 from Global Forest Change (Hansen)",
                  "Proportion cropland cover 2015 from Coperncius 100m fractional land cover",
                  "Proportion area experiencing net cropland expansion between 2000 and 2019 from GLAD (Global Cropland Expansion)",
                  "Proportion built-up cover 2015 from Coperncius 100m fractional land cover",
                  "Proportion area experiencing built-up land expansion 2000-2019 from ESA-CCI Land Cover",
                  "Landscape heterogeneity 1: Variance (coefficient of variation) in peak EVI (2001-2005) among all 250m pixels",
                  "Landscape heterogeneity 2: Second-order dissimilarity in peak EVI (2001-2005) among all 250m pixels",
                  "Mean density of livestock (cattle, sheep, goats, buffalo, pigs) in animals per km2 2010, from Gridded Livestock of the World v3",
                  "Mean motorised travel time to the nearest health centre")
)

# extract each covariate
covs = covs_df$cov
ddx = dd
for(c in covs){
  print(c)
  xx = extractCovariateVals(dd, covariate_name=c) # stored in '00_covar_extraction_funcs.R'
  ddx = cbind(ddx, xx)
}

# log transform any relevant variables and save
ddx = ddx %>%
  dplyr::mutate(health_travel_log = log(health_travel+1),
                livestock_log = log(livestock_all+1)) %>% 
  st_drop_geometry()
write.csv(ddx, "./output/model_df/lassa_df_wp.csv", row.names=FALSE)

# # view boxplot of covars
plot_covars = ddx %>%
  reshape2::melt(id.vars = 1:5) %>%
  #dplyr::filter(!variable %in% c("healthtraveltime", "popdens", "y")) %>%
  ggplot() +
  ggforce::geom_sina(aes(x = factor(presence), y=value, group=factor(presence)), width=1, size=0.1, col="grey60") +
  geom_boxplot(aes(x = factor(presence), y=value, group=factor(presence), fill=variable), width=0.25, alpha=0.8) +
  facet_wrap(~variable, scales="free_y") +
  theme_classic() +
  xlab("Presence") +
  ylab("Covariate value") +
  theme(legend.position="none",
        strip.background = element_blank(), strip.text = element_text(size=13))
ggsave(plot_covars, file="./lasv_covars.png", device="png", units="in", dpi=600, width=10, height=10)






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

# data, keep complete cases and scale predictors
dd = read.csv("./output/model_df/lassa_df_wp.csv") 
dd = dd[ complete.cases(dd), ]
dd[ , 6:ncol(dd)] = apply(dd[ , 6:ncol(dd)], 2, scale)

# # check for multicollinearity - urban cover and urban expansion are reasonably correlated
usdm::vifstep(dd[ , 6:ncol(dd) ] %>% dplyr::select(-urban_expansion, -health_travel, -livestock_all))
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
max.edge = 0.8
bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially

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
                           prior.range = c(30, 0.1), 
                           prior.sigma = c(5, 0.01))

# create projector matrix for point locations based on mesh
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

# response variable: cases, with total population as offset
dd$y = dd$presence

# error family
family1 = "binomial"

# create stack following haakon: https://haakonbakka.bitbucket.io/btopic107.html#3_stationary_models
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
predictors = covs_df$cov
predictors = c(predictors[ -which(predictors %in% c("health_travel", "livestock_all")) ], "health_travel_log", "livestock_log")

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

# multivariate model (all preds except highly collinear)
form = formula(paste(c(form_base, predictors[ -which(predictors %in% c("evi_coefvar", "urban_expansion")) ]), collapse=" + "))
m2 = fitINLAModel(formx=form, family=family1)

# multivariate model (causal subset)
preds_c = c("tmean_anomaly", "precip_anomaly", "crop_cover", 
            "crop_expansion", "urban_cover", "health_travel_log", "evi_dissimilarity")
form = formula(paste(c(form_base, preds_c), collapse=" + "))
m3 = fitINLAModel(formx=form, family=family1)




# ------------ view model results --------

# extract and view spatial field SPDE from full model
spde_plot = rasteriseSPDE(m3$summary.random$s$mean, mesh1, study_area) %>%
  as.data.frame(xy=TRUE) %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=field)) +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black") + 
  geom_point(data=dd[ dd$presence == 1, ], aes(Longitude, Latitude), col="red", size=0.8) + 
  theme(legend.position = c(0.98, 0.15), plot.title = element_text(hjust=0.5, size=16)) +
  ggtitle("Lassa fever")


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
  dplyr::mutate(Disease = "Lassa fever",
                Num_observations = sum(dd$presence))
row.names(fx) = c()
write.csv(fx, "./output/model_outputs/disease_params/lassa_params.csv", row.names=FALSE)


# plot parameters from 3 models
params_plot = fx %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::mutate(param = replace(param, param == "health_travel_log", "Travel time to\nhealth facility (log)"),
                param = replace(param, param == "livestock_log", "Livestock\ndensity (log)"),
                param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
                param = replace(param, param == "forest_cover", "Forest\ncover"),
                param = replace(param, param == "forest_loss", "Forest\nloss"),
                param = replace(param, param == "crop_expansion", "Cropland\nexpansion"),
                param = replace(param, param == "evi_dissimilarity", "EVI\ndissimilarity"),
                param = replace(param, param == "evi_coefvar", "EVI\nvariance"),
                param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
                param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
                param = replace(param, param == "urban_cover", "Urban\ncover"),
                param = replace(param, param == "crop_cover", "Cropland\ncover"),
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
  theme(legend.position = "none",
        axis.text = element_text(size=12), 
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=13), 
        strip.background = element_blank(),
        strip.text = element_text(size=13)) +
  scale_color_viridis_d(begin=0, end=0.8, direction=-1) + 
  facet_wrap(~type, nrow=3)

# uni + multi for colin's talk
# params_plot = fx %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::filter(type != "Causal") %>%
#   dplyr::mutate(param = replace(param, param == "healthtraveltime_log", "Travel time to\nhealth facility (log)"),
#                 param = replace(param, param == "livestock_log", "Livestock\ndensity (log)"),
#                 param = replace(param, param == "livestock_keepers", "Livestock\ndensity (log)"), # temporary
#                 param = replace(param, param == "tmean_anomaly", "Tmean\nanomaly"),
#                 param = replace(param, param == "forest_cover", "Forest\ncover"),
#                 param = replace(param, param == "forest_loss", "Forest\nloss"),
#                 param = replace(param, param == "agri_expansion", "Agriculture\nexpansion"),
#                 param = replace(param, param == "precip_anomaly", "Precip\nanomaly"),
#                 param = replace(param, param == "urban_expansion", "Urban\nexpansion"),
#                 param = replace(param, param == "urban_cover", "Urban\ncover"),
#                 param = replace(param, param == "agri_cover", "Agriculture\ncover"),
#                 param = replace(param, param == "popdens_log", "Population\ndensity (log)")) %>%
#   dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", ordered=TRUE))) %>%
#   #dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", "Causal"), ordered=TRUE)) %>%
#   ggplot() + 
#   geom_hline(yintercept=0, lty=2) +
#   geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) + 
#   geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.6, alpha=0.7) +
#   theme_classic() + 
#   ylab("Posterior marginal slope (log odds)") + 
#   xlab("") + 
#   theme(legend.position = "top",
#         axis.text = element_text(size=13), 
#         axis.title = element_text(size=13),
#         legend.title = element_blank(),
#         legend.text = element_text(size=13)) +
#   scale_color_viridis_d(begin=0, end=0.8, direction=-1)

pp = gridExtra::grid.arrange(spde_plot, params_plot, ncol=1, heights=c(1, 1.4))
ggsave(pp, file="./output/plots/lassa_dashboardtest_jun22.png", device="png", units="in", width=13, height=13, dpi=600, scale=0.9)




