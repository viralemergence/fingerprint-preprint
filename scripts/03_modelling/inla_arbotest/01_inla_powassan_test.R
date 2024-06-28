


# ================= Powassan ====================

# Comparison of results for sufficient surveillance data
# 1. Spatial model with case counts (nbinomial) across all polygons - transmission intensity
# 2. Fingerprint detection model



# ------------ housekeeping ----------------

# setup objects
setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# random seed and disease name
rseed = 68302
disease_name = "powassan"

# define important analysis variables for specific disease
occurrence_point_buffer_km = 5 # km buffer around each occurrence point for spatial uncertainty in infection event
filter_area_thresh = 5000 # maximum geographical area of occurrence polygons to include in km2 (i.e. exclude v imprecise records)
num_background_points = nrow(dz) * 5 # number of background points to generate
covar_coverage_thresh = 0.9 # only keep covariates that are non-NA for at least X% of points
zeroes_coverage_thresh = 0.95 # only keep covariates that are non-zero for at least X% of points


# ---------------- Disease data ---------------

# read disease data
dz = read.csv("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet.csv") %>%
  dplyr::filter(Disease == "Powassan virus disease" & Presentation == "All cases") %>%
  dplyr::mutate(record_id = 1:length(Disease), 
                ADMcode = as.character(as.integer(ADMcode)))
load("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet_shp.R")
dz_shp = shp %>%
  dplyr::mutate(ADMcode = as.character(as.integer(ADMcode))) %>%
  dplyr::filter(ADMcode %in% dz$ADMcode) %>%
  st_transform(crs = 4326)
rm(shp)

# create study area polygon with US border boundary
bounds = st_read("./data/shapefiles/cb_2018_us_state_5m.shp") %>% 
  dplyr::filter(!NAME %in% c("Hawaii", "Guam", "American Samoa", "Puerto Rico", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands", "Alaska")) %>%
  st_transform(crs = 4326)
bounds = st_union(bounds)
bounds = st_sf(id = 1, geometry = bounds)
study_area = createStudyAreaPolygon(lat=dz$Latitude, lon=dz$Longitude, custom_boundary=bounds)

# top ranked hypotheses
hyp_strict = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv") %>%
  dplyr::filter(covar_name != "") %>%
  dplyr::filter(dz_name == dz$Disease[1]) %>%
  dplyr::filter(TestCon == "Test") %>%
  dplyr::select(covar_name)
hyp_strict = as.vector(hyp_strict$covar_name)
hyp_strict = c("urban_cover", hyp_strict)

# covariates for USA polygons and set standard ADMcode
covars = read.csv("./output/model_df/usa_admin2_covars.csv")
covars$ADMcode = paste(
  sprintf("%02d", as.numeric( unlist(lapply(strsplit(covars$ADMcode, "_"), "[", 1)) )),
  sprintf("%03d", as.numeric( unlist(lapply(strsplit(covars$ADMcode, "_"), "[", 2)) )),
  sep=""
)

# log covars and save for later
covars$health_travel_log = log(covars$health_travel + 1)
covars$livestock_log = log(covars$livestock_all +  1)
covs2 = covars

# scale covars
covars[ , 5:ncol(covars) ] = apply(covars[ , 5:ncol(covars) ], 2, scale) 

# keep only hypothesised covars and remove highly collinear (urb exp)
covars = covars[ , which(names(covars) %in% c("ADMcode", hyp_strict)) ]
hyp_strict = hyp_strict[ ! hyp_strict %in% c("urban_expansion", "biodiv_intact") ]





# ------------ Create model dataframe ---------------

# all counties
shp = sf::st_read("./data/shapefiles/usa_counties.shp")
shp = sf::st_transform(shp, crs=sf::st_crs(study_area))

# exclude non continental counties
to_exclude = sf::st_read("./data/shapefiles/cb_2018_us_state_5m.shp") %>% 
  dplyr::filter(NAME %in% c("Hawaii", "Guam", "American Samoa", "Puerto Rico", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands", "Alaska")) %>%
  st_drop_geometry()
shp = shp %>%
  dplyr::filter(!STATEFP %in% to_exclude$STATEFP)

# identify counties that intersect study area
ii = sf::st_intersects(study_area, shp)
ll = unlist(ii)
shp = shp[ ll, ]

# exclude outsized to match other pipeline
shp$area = as.vector(sf::st_area(shp)) / 10^6
shp = shp %>% dplyr::filter(area < filter_area_thresh)

# all matched on covars?
all(shp$ADMcode %in% covars$ADMcode)

# set covars and shp admcodes to integer for matching
shp$ADMcode = as.character(as.integer(shp$ADMcode))
covars$ADMcode = as.character(as.integer(covars$ADMcode))

# create dataframe of all years by locations (for binomial model)
df = expand.grid(shp$ADMcode, 2004:2020) %>%
  as.data.frame() %>%
  dplyr::rename(
    "ADMcode"=1, "Year"=2
  ) %>%
  dplyr::left_join(
    dz[ , c("ADMcode", "Year", "NumCases") ]
  ) %>%
  dplyr::mutate(NumCases = replace(NumCases, is.na(NumCases), 0),
                Outbreak = ifelse(NumCases>0, 1, 0)) %>%
  dplyr::left_join(covars) %>%
  dplyr::left_join(
    shp %>% sf::st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)
  )

# create dataframe of total cases per county (for numcases model)
df1 = data.frame( ADMcode = shp$ADMcode ) %>%
  dplyr::left_join(
    dz %>% dplyr::group_by(ADMcode) %>% dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE))
  ) %>%
  dplyr::mutate(NumCases = replace(NumCases, is.na(NumCases), 0),
                Outbreak = ifelse(NumCases>0, 1, 0)) %>%
  dplyr::left_join(covars) %>%
  dplyr::left_join(
    shp %>% sf::st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)
  )

# extract population for offset in incidence model
pop = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_1km_Aggregated.tif")
pop_ext = exactextractr::exact_extract(pop, shp, fun="sum")
pop_df = data.frame(ADMcode = shp$ADMcode, population = pop_ext, logpop = log(pop_ext / 100000))
df = dplyr::left_join(df, pop_df)
df1 = dplyr::left_join(df1, pop_df)

# code year
df$Year = as.integer(as.factor(df$Year))



# ----------- set up modelling objects for all models ---------

# priors
hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=1, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for fixed effects


# ------------- check data ---------------

# shp %>%
#   dplyr::left_join(
#     df1[ c("ADMcode", "NumCases")]
#   ) %>%
#   dplyr::mutate(NumCases = replace(NumCases, NumCases == 0, NA)) %>%
#   ggplot() +
#   geom_sf(aes(fill=log(NumCases+1)), color=NA) +
#   theme_void() +
#   scale_fill_viridis_c(option="magma", direction=-1, na.value="white")
# 
# df %>%
#   dplyr::select(ADMcode, Year, NumCases) %>%
#   dplyr::group_by(ADMcode) %>%
#   dplyr::summarise(n_outbreaks = sum(NumCases>0)) %>%
#   dplyr::mutate(n_outbreaks = replace(n_outbreaks, n_outbreaks == 0, NA)) %>%
#   dplyr::full_join(shp) %>%
#   sf::st_as_sf() %>%
#   ggplot() +
#   geom_sf(aes(fill=n_outbreaks), color=NA) +
#   theme_void() +
#   scale_fill_viridis_c(option="magma", direction=-1, na.value="white")



# ----------- MODEL 1: Annual case incidence ---------

# build mesh using spatial locs
spatial_locs = df1[ , c("Longitude", "Latitude")]
max.edge = 0.7
bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
plot(mesh1)

# create projector matrix for point locations based on mesh
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

# define pc matern spde model
spde = inla.spde2.pcmatern(mesh1, 
                           prior.range = c(10, 0.1), 
                           prior.sigma = c(5, 0.01))

# build stack 
df1$y = df1$NumCases
stack1 = inla.stack(data=list(y=df1$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, df1[ , -which(names(df1) %in% c("y"))]),
                                 iidx = 1:nrow(df1)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')

form = paste(c("y ~ -1 + Intercept + offset(logpop) + f(s, model=spde)", hyp_strict), collapse=" + ")
m1 = fitINLAModel(formx = formula(form), family="zeroinflatednbinomial1", stack=stack1, spde=spde, inla.mode="experimental")
save(m1, file="./output/model_outputs/wnv_testmodels/powassan_spde_mod1.R")

# viz
df1$fitted = exp(m1$summary.linear.predictor$mean[1:nrow(df1)] + df1$logpop)
plot(df1$fitted, df1$NumCases)

# predicted spatial field
spde_pred = rasteriseSPDE(m1$summary.random$s$mean, mesh1, study_area)
spde_pred %>%
  as.data.frame(xy=TRUE) %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=field)) + 
  scale_fill_gradient2() + 
  theme_void()



# # ------------ MODEL 2: Annual presence absence all data -------------
# 
# # build mesh using spatial locs
# spatial_locs = df[ , c("Longitude", "Latitude")]
# max.edge = 1.2
# bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
# mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
#                      loc.domain=spatial_locs,
#                      max.edge=c(1, 5)*max.edge,
#                      cutoff = max.edge/5,
#                      offset = c(max.edge, bound.outer))
# plot(mesh1)
# 
# # create projector matrix for point locations based on mesh
# A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))
# 
# # define pc matern spde model
# spde = inla.spde2.pcmatern(mesh1, 
#                            prior.range = c(20, 0.1), 
#                            prior.sigma = c(5, 0.01))
# 
# # build stack 
# # df = df %>%
# #   dplyr::left_join(shp %>% sf::st_drop_geometry() %>% dplyr::select(ADMcode, area_km2)) %>%
# #   dplyr::mutate(popdens = population / area_km2, 
# #                 logpopdens = log(popdens))
# df$y = df$Outbreak
# 
# stack1 = inla.stack(data=list(y=df$y),
#                     effects=list(s = 1:mesh1$n, # spatial
#                                  data.frame(Intercept=1, df[ , -which(names(df) %in% c("y"))]),
#                                  iidx = 1:nrow(df)),
#                     A=list(A, 1, 1),
#                     remove.unused = FALSE, tag='est')
# 
# form = paste(c("y ~ -1 + Intercept + f(Year, model='iid') + f(s, model=spde)", hyp_strict), collapse=" + ")
# m2 = fitINLAModel(formx = formula(form), family="binomial", stack=stack1, spde=spde)
# save(m2, file="./output/model_outputs/wnv_testmodels/wnv_spde_mod2.R")
# 
# # predicted spatial field
# spde_pred = rasteriseSPDE(m2$summary.random$s$mean, mesh1, study_area)
# spde_pred %>%
#   as.data.frame(xy=TRUE) %>%
#   ggplot() + 
#   geom_raster(aes(x, y, fill=field)) + 
#   scale_fill_gradient2() + 
#   theme_void()



# ------------ MODEL 2: Annual presence + pseudoabsences (i.e. detection model) ----------

filename = paste("./output/model_df/" , disease_name, "_df_wp.csv", sep="")
dd = read.csv(filename) %>%
  dplyr::select(-health_travel, -livestock_all, -evi_coefvar, -tmean_anomaly, -precip_anomaly)

# proportion NAs in each predictor: exclude predictors that do not have coverage for at least X% of data
nas = apply(dd[ dd$presence == 1, 6:ncol(dd)], 2, function(x){ sum(!is.na(x))/length(x) })
nas = names(nas[ nas < covar_coverage_thresh ])

# excluding zeroes/nas
dd = dd %>% dplyr::select(-all_of(nas))

# scale predictors
dd[ , 6:ncol(dd)] = apply(dd[ , 6:ncol(dd)], 2, scale)

# response variable: 1/0 outbreaks
dd$y = dd$presence
family1 = "binomial"

# priors
hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=1, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for fixed effects

# build mesh using spatial locs
spatial_locs = dd[ , c("Longitude", "Latitude")]
max.edge = 0.7
bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain=spatial_locs,
                     max.edge=c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))
plot(mesh1)

# create projector matrix for point locations based on mesh
A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))

# define pc matern spde model
spde = inla.spde2.pcmatern(mesh1, 
                           prior.range = c(10, 0.1), 
                           prior.sigma = c(5, 0.01))

# build stack 
stack1 = inla.stack(data=list(y=dd$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, dd[ , -which(names(dd) %in% c("y"))]),
                                 iidx = 1:nrow(dd)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')

form_base = "y ~ -1 + Intercept + f(s, model=spde)"
form = formula(paste(c(form_base, hyp_strict), collapse=" + "))
m2 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde)
save(m2, file="./output/model_outputs/wnv_testmodels/powassan_spde_mod2.R")





# # ----------- examine params -----------
# 
# load(file="./output/model_outputs/wnv_testmodels/powassan_spde_mod1.R")
# load(file="./output/model_outputs/wnv_testmodels/powassan_spde_mod2.R")
# 
# results = do.call(
#   rbind.data.frame,
#   list(extractFixedINLA(m1, model_name = "Incidence"),
#        extractFixedINLA(m2, model_name = "Annual outbreak (pseudoabsences)")
#   )) %>%
#   dplyr::filter(param != "Intercept")
# 
# ggplot(results) +
#   geom_point(aes(param, mean, group=model, color=model), position=position_dodge(width=0.3)) +
#   geom_linerange(aes(param, ymin=lower, ymax=upper, group=model, color=model), position=position_dodge(width=0.3)) +
#   theme_classic() +
#   geom_hline(yintercept=0,lty=2)
# 








# 
# 
# 
# 
# 
# # ----------- MODEL 1: Case incidence by county  ----------------
# 
# # priors
# hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
# control.fixed1 = list(mean.intercept=0, # prior mean for intercept
#                       prec.intercept=1, # prior precision for intercept
#                       mean=0, # prior mean for fixed effects
#                       prec=1)  # prior precision for fixed effects
# 
# # modelling objects
# df$yearx = as.integer(as.factor(df$Year))
# 
# # # model 1: total case incidence (negative binomial) with space modelled as a CAR model
# # form = paste(c("y ~ 1 + offset(logpop) + f(polyid, model='bym2', graph='adjmatrix_uscounties')", hyp_strict), collapse=" + ")
# # df1$y = df1$NumCases
# # m1 = inla(formula(form),
# #           verbose = FALSE,
# #           data = df1,
# #           family = "nbinomial",
# #           control.fixed = control.fixed1, 
# #           control.predictor=list(compute=TRUE, link=1),
# #           control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=FALSE, return.marginals=TRUE),
# #           control.inla = list(strategy='adaptive', # adaptive gaussian
# #                               cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
# #           inla.mode = "experimental")
# # save(m1, file="./output/model_outputs/wnv_testmodels/wnv_test_model1.R")
# 
# # model 1b: annual case incidence (negative binomial) with "true" zeroes
# form = paste(c("y ~ 1 + offset(logpop) + f(yearx, model='rw1', hyper=hyper1.rw) + f(polyid, model='bym2', graph='adjmatrix_uscounties')", hyp_strict), collapse=" + ")
# df$y = df$NumCases
# m1b = inla(formula(form),
#           verbose = FALSE,
#           data = df,
#           family = "nbinomial",
#           control.fixed = control.fixed1, 
#           control.predictor=list(compute=TRUE, link=1),
#           control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=FALSE, return.marginals=TRUE),
#           control.inla = list(strategy='adaptive', # adaptive gaussian
#                               cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
#           inla.mode = "experimental")
# save(m1b, file="./output/model_outputs/wnv_testmodels/wnv_test_model1b.R")
# 
# 
# 
# # -------------- MODEL 2: Annual presence/absence --------------
# 
# # model 1b: annual case incidence (negative binomial) with "true" zeroes
# form = paste(c("y ~ 1 + f(yearx, model='rw1', hyper=hyper1.rw) + f(polyid, model='bym2', graph='adjmatrix_uscounties')", hyp_strict), collapse=" + ")
# df$y = ifelse(df$NumCases>0, 1, 0)
# m2 = inla(formula(form),
#            verbose = FALSE,
#            data = df,
#            family = "binomial",
#            control.fixed = control.fixed1, 
#            control.predictor=list(compute=TRUE, link=1),
#            control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=FALSE, return.marginals=TRUE),
#            control.inla = list(strategy='adaptive', # adaptive gaussian
#                                cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
#            inla.mode = "experimental")
# save(m2, file="./output/model_outputs/wnv_testmodels/wnv_test_model2.R")
# 
# 
# 
# # -------------- MODEL 3: Annual outbreak indicator with polygon pseudoabsences ---------------
# 
# # set up dataframe 
# set.seed(2048)
# df2 = rbind(
#   df[ df$Outbreak == 1, ],
#   df %>% dplyr::sample_n(9000) %>% dplyr::mutate(Outbreak = 0)
# ) %>%
#   dplyr::left_join(
#     shp %>% sf::st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)
#   )
# 
# # build mesh using spatial locs
# spatial_locs = df2[ , c("Longitude", "Latitude")]
# max.edge = 1.25
# bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
# mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
#                      loc.domain=spatial_locs,
#                      max.edge=c(1, 5)*max.edge,
#                      cutoff = max.edge/5,
#                      offset = c(max.edge, bound.outer))
# plot(mesh1)
# 
# # create projector matrix for point locations based on mesh
# A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))
# 
# # define pc matern spde model
# spde = inla.spde2.pcmatern(mesh1, 
#                            prior.range = c(30, 0.1), 
#                            prior.sigma = c(5, 0.01))
# 
# # build stack 
# df2$y = df2$Outbreak
# stack1 = inla.stack(data=list(y=df2$y),
#                     effects=list(s = 1:mesh1$n, # spatial
#                                  data.frame(Intercept=1, df2[ , -which(names(df2) %in% c("y"))]),
#                                  iidx = 1:nrow(df2)),
#                     A=list(A, 1, 1),
#                     remove.unused = FALSE, tag='est')
# 
# form = paste(c("y ~ -1 + Intercept + f(s, model=spde)", hyp_strict), collapse=" + ")
# m3 = fitINLAModel(formx = formula(form), family="binomial", stack=stack1, spde=spde)
# save(m3, file="./output/model_outputs/wnv_testmodels/wnv_test_model3.R")
# 
# 
# 
# # ------------- MODEL 4: Annual outbreak indicator with pw pseudoabsences -------------
# 
# filename = paste("./output/model_df/" , disease_name, "_df_wp.csv", sep="")
# dd = read.csv(filename) %>%
#   dplyr::select(-health_travel, -livestock_all, -evi_coefvar, -tmean_anomaly, -precip_anomaly)
# 
# # proportion NAs in each predictor: exclude predictors that do not have coverage for at least X% of data
# nas = apply(dd[ dd$presence == 1, 6:ncol(dd)], 2, function(x){ sum(!is.na(x))/length(x) })
# nas = names(nas[ nas < covar_coverage_thresh ])
# 
# # excluding zeroes/nas
# dd = dd %>% dplyr::select(-all_of(nas))
# 
# # scale predictors
# dd[ , 6:ncol(dd)] = apply(dd[ , 6:ncol(dd)], 2, scale)
# 
# # response variable: 1/0 outbreaks
# dd$y = dd$presence
# family1 = "binomial"
# 
# # priors
# hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
# control.fixed1 = list(mean.intercept=0, # prior mean for intercept
#                       prec.intercept=1, # prior precision for intercept
#                       mean=0, # prior mean for fixed effects
#                       prec=1)  # prior precision for fixed effects
# 
# # build mesh using spatial locs
# spatial_locs = dd[ , c("Longitude", "Latitude")]
# max.edge = 1.25
# bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
# mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
#                      loc.domain=spatial_locs,
#                      max.edge=c(1, 5)*max.edge,
#                      cutoff = max.edge/5,
#                      offset = c(max.edge, bound.outer))
# plot(mesh1)
# 
# # create projector matrix for point locations based on mesh
# A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))
# 
# # define pc matern spde model
# # prior median range at approximately 0.5*diff(range(spatial_locs$Longitude)) - roughly extent of study area
# # prior probability of marginal sigma of 3 or more is 0.01
# # (iterative tuning to ensure SPDE fits reasonably)
# spde = inla.spde2.pcmatern(mesh1, 
#                            prior.range = c(30, 0.1), 
#                            prior.sigma = c(5, 0.01))
# 
# # build stack 
# stack1 = inla.stack(data=list(y=dd$y),
#                     effects=list(s = 1:mesh1$n, # spatial
#                                  data.frame(Intercept=1, dd[ , -which(names(dd) %in% c("y"))]),
#                                  iidx = 1:nrow(dd)),
#                     A=list(A, 1, 1),
#                     remove.unused = FALSE, tag='est')
# 
# form_base = "y ~ -1 + Intercept + f(s, model=spde)"
# form = formula(paste(c(form_base, hyp_strict), collapse=" + "))
# m4 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde)
# save(m4, file="./output/model_outputs/wnv_testmodels/wnv_test_model4.R")
# 
# 
# 
# # ----------- examine params -----------
# 
# results = do.call(
#   rbind.data.frame,
#   list(extractFixedINLA(m1b, model_name = "Incidence"),
#        extractFixedINLA(m2, model_name = "Outbreak ('true' absences)"),
#        extractFixedINLA(m3, model_name = "Outbreak (pseudoabsence polygons)"),
#        extractFixedINLA(m4, model_name = "Outbreak (pseudoabsence points)")
#   )) %>%
#   dplyr::filter(param != "Intercept")
# 
# ggplot(results) +
#   geom_point(aes(param, mean, group=model, color=model), position=position_dodge(width=0.3)) +
#   geom_linerange(aes(param, ymin=lower, ymax=upper, group=model, color=model), position=position_dodge(width=0.3)) +
#   theme_classic() +
#   geom_hline(yintercept=0,lty=2)
# 
# 
# 
# 
# # # --------------- Generate background point absences (pop weighted and random) ----------------
# # 
# # # generate population weighted point+buffers with same median polygon area
# # poly_area = as.vector(sf::st_area(dz_shp) / 10^6)
# # med_area_size = median(poly_area)
# # buf_rad = sqrt(med_area_size / pi) * 1000
# # 
# # # weighted points, extract covars, save
# # bg_weighted = generateBackgroundPolygons(num_points = num_background_points,
# #                                          study_area = study_area,
# #                                          buffer_radius = buf_rad,
# #                                          method = "popweight",
# #                                          seed = rseed)
# # 
# # hyp2 = hyp_strict
# # hyp2 = replace(hyp2, hyp2 == "livestock_log", "livestock_all")
# # hyp2 = replace(hyp2, hyp2 == "health_travel_log", "health_travel")
# # for(cc in hyp2){
# #   print(cc)
# #   extr_x = extractCovariateVals(bg_weighted, covariate_name=cc)
# #   bg_weighted = cbind(bg_weighted, extr_x)
# # }
# # 
# # write.csv(
# #   bg_weighted %>% sf::st_drop_geometry(),
# #   "./output/model_outputs/wnv_testmodels/covars_bg_popweighted.csv",
# #   row.names=FALSE
# # )
# # 
# # 
# # # random points, same
# # bg_random = generateBackgroundPolygons(num_points = num_background_points,
# #                                          study_area = study_area,
# #                                          buffer_radius = buf_rad,
# #                                          method = "random",
# #                                          seed = rseed)
# # 
# # for(cc in hyp2){
# #   print(cc)
# #   extr_x = extractCovariateVals(bg_random, covariate_name=cc)
# #   bg_random = cbind(bg_random, extr_x)
# # }
# # 
# # save
# # write.csv(
# #   bg_random %>% sf::st_drop_geometry(),
# #   "./output/model_outputs/wnv_testmodels/covars_bg_random.csv",
# #   row.names=FALSE
# # )
# # 
# # 
# # 
# # # -------------- Background point models with outbreak indicator ----------------
# # 
# # dfw = read.csv("./output/model_outputs/wnv_testmodels/covars_bg_popweighted.csv")
# # 
# # # add covars to dz
# # dz1 = dz %>% 
# #   dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% 
# #   dplyr::mutate(presence = 1) %>%
# #   dplyr::left_join(
# #     covs2[ , which(names(covs2) %in% names(dfw))] %>%
# #       dplyr::select(-Longitude, -Latitude)
# #   )
# # 
# # # weighted data frame and log/scale covars
# # dfw = read.csv("./output/model_outputs/wnv_testmodels/covars_bg_popweighted.csv") %>%
# #   dplyr::mutate(presence = 0) %>%
# #   #sample_n(nrow(dz1)) %>%
# #   rbind(dz1) %>%
# #   dplyr::mutate(health_travel_log = log(health_travel+1))
# # dfw[ , which(names(dfw) %in% causal_vars) ] = apply(dfw[ , which(names(dfw) %in% causal_vars) ], 2, scale)
# # 
# # # random data frame
# # dfr = read.csv("./output/model_outputs/wnv_testmodels/covars_bg_random.csv") %>%
# #   dplyr::mutate(presence = 0) %>%
# #   #sample_n(nrow(dz1)) %>%
# #   rbind(dz1) %>%
# #   dplyr::mutate(health_travel_log = log(health_travel+1))
# # dfr[ , which(names(dfr) %in% causal_vars) ] = apply(dfr[ , which(names(dfr) %in% causal_vars) ], 2, scale)
# # 
# # 
# # 
# # # ------------ model 3: population weighted background points --------------
# #   
# # # y
# # dfw$y = dfw$presence
# # 
# # # build mesh using spatial locs
# # spatial_locs = dfw[ , c("Longitude", "Latitude")]
# # max.edge = 1.25
# # bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
# # mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
# #                      loc.domain=spatial_locs,
# #                      max.edge=c(1, 5)*max.edge,
# #                      cutoff = max.edge/5,
# #                      offset = c(max.edge, bound.outer))
# # plot(mesh1)
# # 
# # # create projector matrix for point locations based on mesh
# # A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))
# # 
# # # define pc matern spde model
# # spde = inla.spde2.pcmatern(mesh1, 
# #                            prior.range = c(30, 0.1), 
# #                            prior.sigma = c(5, 0.01))
# # 
# # # build stack 
# # stack1 = inla.stack(data=list(y=dfw$y),
# #                     effects=list(s = 1:mesh1$n, # spatial
# #                                  data.frame(Intercept=1, dfw[ , -which(names(dfw) %in% c("y"))]),
# #                                  iidx = 1:nrow(dfw)),
# #                     A=list(A, 1, 1),
# #                     remove.unused = FALSE, tag='est')
# # 
# # form = paste(c("y ~ -1 + Intercept + f(s, model=spde)", causal_vars), collapse=" + ")
# # m3 = fitINLAModel(formx = formula(form), family="binomial", stack=stack1, spde=spde)
# # save(m3, file="./output/model_outputs/wnv_testmodels/wnv_test_model3.R")
# # 
# # 
# # 
# # # ------------ model 4: random background points --------------
# # 
# # # y
# # dfr$y = dfr$presence
# # 
# # # build mesh using spatial locs
# # spatial_locs = dfr[ , c("Longitude", "Latitude")]
# # max.edge = 1.25
# # bound.outer = diff(range(spatial_locs$Longitude))/5 # outer bound: 1/3 of range diff initially
# # mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
# #                      loc.domain=spatial_locs,
# #                      max.edge=c(1, 5)*max.edge,
# #                      cutoff = max.edge/5,
# #                      offset = c(max.edge, bound.outer))
# # plot(mesh1)
# # 
# # # create projector matrix for point locations based on mesh
# # A = inla.spde.make.A(mesh1, loc=as.matrix(spatial_locs))
# # 
# # # define pc matern spde model
# # spde = inla.spde2.pcmatern(mesh1, 
# #                            prior.range = c(30, 0.1), 
# #                            prior.sigma = c(5, 0.01))
# # 
# # # build stack 
# # stack1 = inla.stack(data=list(y=dfr$y),
# #                     effects=list(s = 1:mesh1$n, # spatial
# #                                  data.frame(Intercept=1, dfr[ , -which(names(dfr) %in% c("y"))]),
# #                                  iidx = 1:nrow(dfr)),
# #                     A=list(A, 1, 1),
# #                     remove.unused = FALSE, tag='est')
# # 
# # form = paste(c("y ~ -1 + Intercept + f(s, model=spde)", causal_vars), collapse=" + ")
# # m4 = fitINLAModel(formx = formula(form), family="binomial", stack=stack1, spde=spde)
# # save(m4, file="./output/model_outputs/wnv_testmodels/wnv_test_model4.R")
# # 
# # 
# # 
# # 
# # 
# # # ----------- examine params -----------
# # 
# # results = do.call(
# #   rbind.data.frame,
# #   list(extractFixedINLA(m1.1, model_name = "Incidence"),
# #        extractFixedINLA(m2, model_name = "Outbreak ('true' absences)"),
# #        extractFixedINLA(m3, model_name = "Outbreak (pseudoabsences pop-weighted)"),
# #        extractFixedINLA(m4, model_name = "Outbreak (pseudoabsences random)")
# #   )) %>%
# #   dplyr::filter(param != "Intercept")
# # 
# # ggplot(results) +
# #   geom_point(aes(param, mean, group=model, color=model), position=position_dodge(width=0.3)) + 
# #   geom_linerange(aes(param, ymin=lower, ymax=upper, group=model, color=model), position=position_dodge(width=0.3)) + 
# #   theme_classic() + 
# #   geom_hline(yintercept=0,lty=2) 
# # 
# # 
# # 
# # # ------- maps of sampling design -----------
# # 
# # # model 1
# # map1 = shp %>% 
# #   dplyr::filter(ADMcode %in% df1$ADMcode) %>%
# #   dplyr::left_join(df1[ c("ADMcode", "NumCases", "population")]) %>%
# #   dplyr::mutate(incidence = NumCases / population * 100000) %>%
# #   ggplot() + 
# #   geom_sf(aes(fill=log(incidence+1)), color=NA) + 
# #   maptheme + 
# #   scale_fill_gradientn(colors=rev(viridisLite::mako(200)), name="Incidence\n(log)") + 
# #   theme(legend.position=c(0.9, 0.25)) + 
# #   ggtitle("Incidence (cases per 100,000 inhabitants)")
# # 
# # # model 2
# # map2 = shp %>% 
# #   dplyr::filter(ADMcode %in% df1$ADMcode) %>%
# #   dplyr::left_join(
# #     df %>%
# #       dplyr::group_by(ADMcode) %>%
# #       dplyr::summarise(n_outbreaks = sum(Outbreak))
# #   ) %>%
# #   #dplyr::mutate(incidence = NumCases / population * 100000) %>%
# #   ggplot() + 
# #   geom_sf(aes(fill=n_outbreaks), color=NA) + 
# #   maptheme + 
# #   scale_fill_gradientn(colors=rev(viridisLite::mako(200)), name="Num\noutbreak\nyears") + 
# #   theme(legend.position=c(0.9, 0.25)) +   
# #   ggtitle("Outbreaks ('true' absences)")
# # 
# # # model 3
# # map3 = shp %>% 
# #   dplyr::filter(ADMcode %in% df1$ADMcode) %>%
# #   dplyr::left_join(
# #     df %>%
# #       dplyr::group_by(ADMcode) %>%
# #       dplyr::summarise(n_outbreaks = sum(Outbreak))
# #   ) %>%
# #   dplyr::filter(n_outbreaks > 0) %>%
# #   ggplot() + 
# #   geom_point(data = dfw[ dfw$Longitude > -140, ], aes(Longitude, Latitude), size=0.1, color="grey50", alpha=0.7) + 
# #   geom_sf(aes(fill=n_outbreaks), color=NA, alpha=0.6) + 
# #   maptheme + 
# #   scale_fill_gradientn(colors=rev(viridisLite::mako(200)[50:200]), name="Num\noutbreak\nyears") +
# #   theme(legend.position=c(0.9, 0.25)) +   
# #   ggtitle("Outbreaks (pseudoabsences population weighted)")
# # 
# # # model 4
# # map4 = shp %>% 
# #   dplyr::filter(ADMcode %in% df1$ADMcode) %>%
# #   dplyr::left_join(
# #     df %>%
# #       dplyr::group_by(ADMcode) %>%
# #       dplyr::summarise(n_outbreaks = sum(Outbreak))
# #   ) %>%
# #   dplyr::filter(n_outbreaks > 0) %>%
# #   ggplot() + 
# #   geom_point(data = dfr[ dfr$Longitude > -140, ], aes(Longitude, Latitude), size=0.1, color="grey50", alpha=0.7) + 
# #   geom_sf(aes(fill=n_outbreaks), color=NA, alpha=0.6) + 
# #   maptheme + 
# #   scale_fill_gradientn(colors=rev(viridisLite::mako(200)[50:200]), name="Num\noutbreak\nyears") +
# #   theme(legend.position=c(0.9, 0.25)) +   
# #   ggtitle("Outbreaks (pseudoabsences random)")
# # 
# # # combine
# # maps = gridExtra::grid.arrange(grobs=list(map1, map2, map3, map4), ncol=2, nrow=2)
