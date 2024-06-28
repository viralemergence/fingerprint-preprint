

# ================ MPOX ===================


# ------------ housekeeping ----------------

# setup objects
setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")


# ---------------- Disease data ---------------

# random seed and disease name
rseed = 1985
disease_name = "mpx"

# disease data and shapefile post-1985
dz = read.csv("./output/spillovers_processed/spillovers_mpx.csv") %>%
  dplyr::filter(Year >= 1985) %>%
  distinct() %>%
  dplyr::mutate(record_id = 1:length(Disease))
load("./output/spillovers_processed/spillovers_mpx_shp.R")
dz_shp = shp %>% dplyr::filter(ADMcode %in% dz$ADMcode); rm(shp)

# create study area polygon
study_area = createStudyAreaPolygon(lat=dz$Latitude, lon=dz$Longitude)



# ------------ Analysis variables ------------

# define important analysis variables for specific disease
occurrence_point_buffer_km = 5 # km buffer around each occurrence point for spatial uncertainty in infection event
filter_area_thresh = 5000 # maximum geographical area of occurrence polygons to include in km2 (i.e. exclude v imprecise records)
num_background_points = nrow(dz) * 3 # number of background points to generate
covar_coverage_thresh = 0.9 # only keep covariates that are non-NA for at least X% of points
zeroes_coverage_thresh = 0.95 # only keep covariates that are non-zero for at least X% of points

# a priori confounders
confound = c("health_travel_log", "urban_cover")

# hypotheses liberal (from participatory exercise)
hyp_lib = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv") %>%
  dplyr::filter(covar_name != "") %>%
  dplyr::filter(dz_name == dz$Disease[1]) %>%
  dplyr::filter(TestLib == "Test") %>%
  dplyr::select(covar_name) 
hyp_lib = as.vector(hyp_lib$covar_name)

# hypotheses strict (from participatory exercise)
hyp_strict = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv") %>%
  dplyr::filter(covar_name != "") %>%
  dplyr::filter(dz_name == dz$Disease[1]) %>%
  dplyr::filter(TestCon == "Test") %>%
  dplyr::select(covar_name)
hyp_strict = as.vector(hyp_strict$covar_name)

# add a priori confounds
if("forest_loss" %in% hyp_lib ){ hyp_lib = c(hyp_lib, "forest_cover") }
if("forest_loss" %in% hyp_strict ){ hyp_strict = c(hyp_strict, "forest_cover") }
hyp_lib = unique(c(hyp_lib, confound)); hyp_lib = hyp_lib[ hyp_lib != "" ]
hyp_strict = unique(c(hyp_strict, confound))

# # replace with health_travel
hyp_lib = replace(hyp_lib, hyp_lib == "health_travel_log", "health_travel")
hyp_strict = replace(hyp_strict, hyp_strict == "health_travel_log", "health_travel")



# ------------- 3. INLA spatial regression model ---------------

# ------- a. read and prepare data ------

# read data and scale predictors (include log livestock and log health trav)
filename = paste("./output/model_df/" , disease_name, "_df_wp.csv", sep="")
dd = read.csv(filename) %>%
  dplyr::select(-health_travel_log,-livestock_all, -evi_coefvar, -tmean, -precip, -tmean_anomaly, -precip_anomaly)

# transform to hours
dd$health_travel = dd$health_travel/60

# calculate scaling factor
sc = sd(dd$health_travel, na.rm=TRUE)

# proportion NAs in each predictor: exclude predictors that do not have coverage for at least X% of data
nas = apply(dd[ dd$presence == 1, 6:ncol(dd)], 2, function(x){ sum(!is.na(x))/length(x) })
nas = names(nas[ nas < covar_coverage_thresh ])

# proportion 0s in each predictor: exclude predictors with zeroes for more than X% of data
zeroes = apply(dd[ , 6:ncol(dd)], 2, function(x){ sum(!is.na(x) & x == 0)/length(x) })
zeroes = names(zeroes[ zeroes > zeroes_coverage_thresh ])

# excluding zeroes/nas
dd = dd %>% dplyr::select(-all_of(nas), -all_of(zeroes))

# scale predictors
pred_names =names(dd)[6:ncol(dd)]
pred_names = pred_names[ -which(pred_names == "health_travel") ]
dd[ , which(names(dd) %in% pred_names) ] = apply(dd[ , which(names(dd) %in% pred_names) ], 2, scale)

# check for multicollinearity
# usdm::vifstep(dd[ , 6:ncol(dd) ])
# corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_lib))))
# corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_strict))))

# specify which multicollinear vars to exclude in later analyses if any
mc_vars_to_exclude = c("urban_expansion")

# examine VIFs/corplots for multivar and causal var subsets
# corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(-all_of(mc_vars_to_exclude))))
# usdm::vifstep(dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_strict)))



# ------ b. build model objects --------

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
max.edge = 1.3
bound.outer = diff(range(spatial_locs$Longitude))/7 # outer bound: 1/3 of range diff initially
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
                           prior.range = c(50, 0.1), 
                           prior.sigma = c(15, 0.1))

# build stack 
stack1 = inla.stack(data=list(y=dd$y),
                    effects=list(s = 1:mesh1$n, # spatial
                                 data.frame(Intercept=1, dd[ , -which(names(dd) %in% c("y"))]),
                                 iidx = 1:nrow(dd)),
                    A=list(A, 1, 1),
                    remove.unused = FALSE, tag='est')



# ----- c. inference -----

# model 0: random effects only baseline
form_base = "y ~ -1 + Intercept + f(s, model=spde)"
#m0 = fitINLAModel(formx = formula(form_base), family=family1, stack=stack1, spde=spde)

# model 1: univariate models including each covariate individually
predictors = names(dd)[ 6:(ncol(dd)-1)]
#m1 = vector("list", length=length(predictors))
fm1 = data.frame()

# for(i in 1:length(predictors)){
#   
#   pred = predictors[i]
#   print(paste("Fitting", pred, sep=" "))
#   form_i = formula(paste(form_base, pred, sep = " + "))
#   m1_i = fitINLAModel(formx = form_i, family=family1, stack=stack1, spde=spde)
#   fm1 = rbind(fm1, getIC(m1_i, mod_name = pred))
#   m1[[i]] = m1_i
#   
# }

# model 2: hypotheses liberally defined (any proposed association)
preds_m = predictors[ !predictors %in% mc_vars_to_exclude] 
preds_m = preds_m[ preds_m %in% hyp_lib ]
form = formula(paste(c(form_base, preds_m), collapse=" + "))
m2 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde)

# model 3: hypotheses strictly defined (top ranked)
preds_c = preds_m[ preds_m %in% hyp_strict ]
form = formula(paste(c(form_base, preds_c), collapse=" + "))
m3 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde)



# ----- e. extract fixed effects ------

# f1 = do.call(
#   rbind.data.frame, 
#   lapply(m1, function(x) extractFixedINLA(x))
# ) %>%
#   dplyr::mutate(type = "Univariate")

f2 = extractFixedINLA(m2) %>%
  dplyr::mutate(type = "Causal (broad)")

f3 = extractFixedINLA(m3) %>%
  dplyr::mutate(type = "Causal (strict)")

fx = do.call(rbind.data.frame, list(f2, f3)) %>% 
  dplyr::mutate(Disease = dz$Disease[1],
                Num_observations = sum(dd$presence))
row.names(fx) = c()

write.csv(fx, paste("./output/model_outputs/disease_params_htt/", disease_name,"_params.csv", sep=""), row.names=FALSE)

