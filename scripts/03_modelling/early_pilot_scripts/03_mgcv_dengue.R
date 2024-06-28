

# ================ DENGUE (GLOBAL) ===================


# ------------ housekeeping ----------------

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# ---------------- Disease data ---------------

# random seed and disease name
rseed = 90273
disease_name = "dengue"

# disease data and shapefile post-1985
dz = read.csv("./output/spillovers_processed/spillovers_dengue.csv") %>%
  distinct() %>%
  dplyr::mutate(record_id = 1:length(Disease)) %>%
  dplyr::filter(!is.na(Longitude)) %>%
  dplyr::filter(Year >= 1985)

# remov polygons wo admcode
dz = dz[ -which(dz$DataType=="polygon" & is.na(dz$ADMcode) ), ]

# load("./output/spillovers_processed/spillovers_dengue_shp1.R")
# load("./output/spillovers_processed/spillovers_dengue_shp2.R")
# shp = rbind(shp1, shp2)
# dz_shp = shp %>%
#   dplyr::filter(ADMcode %in% dz$ADMcode)
# 
# # create study area polygon with country/coastal boundaries
# study_area = createStudyAreaPolygon(lat=dz$Latitude, lon=dz$Longitude, buffer=FALSE, extent_border = 8)
# 

# ------------ Analysis variables ------------

# define important analysis variables for specific disease
occurrence_point_buffer_km = 10 # km buffer around each occurrence point for spatial uncertainty in infection event
filter_area_thresh = 5000 # maximum geographical area of occurrence polygons to include in km2 (i.e. exclude v imprecise records)
num_background_points = nrow(dz) # number of background points to generate
covar_coverage_thresh = 0.9 # only keep covariates that are non-NA for at least X% of points
zeroes_coverage_thresh = 0.95 # only keep covariates that are non-zero for at least X% of points

# a priori confounders
confound = c("health_travel_log", "urban_cover")

# hypotheses liberal (from participatory exercise)
hyp_lib = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv") %>%
  dplyr::filter(dz_name == dz$Disease[1]) %>%
  dplyr::filter(TestLib == "Test") %>%
  dplyr::select(covar_name) 
hyp_lib = as.vector(hyp_lib$covar_name)

# hypotheses strict (from participatory exercise)
hyp_strict = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv") %>%
  dplyr::filter(dz_name == dz$Disease[1]) %>%
  dplyr::filter(TestCon == "Test") %>%
  dplyr::select(covar_name)
hyp_strict = as.vector(hyp_strict$covar_name)

# add a priori confounds
if("forest_loss" %in% hyp_lib ){ hyp_lib = c(hyp_lib, "forest_cover") }
if("forest_loss" %in% hyp_strict ){ hyp_strict = c(hyp_strict, "forest_cover") }
hyp_lib = unique(c(hyp_lib, confound))
hyp_strict = unique(c(hyp_strict, confound))


# ------- a. read and prepare data ------

# read data and scale predictors (include log livestock and log health trav)
filename = paste("./output/model_df/" , disease_name, "_df_wp.csv", sep="")
dd = read.csv(filename) %>%
  dplyr::select(-health_travel, -livestock_all, -evi_coefvar, -tmean, -precip, -tmean_anomaly, -precip_anomaly)

# proportion NAs in each predictor: exclude predictors that do not have coverage for at least X% of data
nas = apply(dd[ dd$presence == 1, 6:ncol(dd)], 2, function(x){ sum(!is.na(x))/length(x) })
nas = names(nas[ nas < covar_coverage_thresh ])

# proportion 0s in each predictor: exclude predictors with zeroes for more than X% of data
zeroes = apply(dd[ , 6:ncol(dd)], 2, function(x){ sum(!is.na(x) & x == 0)/length(x) })
zeroes = names(zeroes[ zeroes > zeroes_coverage_thresh ])

# excluding zeroes/nas
dd = dd %>% dplyr::select(-all_of(nas), -all_of(zeroes))

# scale predictors
dd[ , 6:ncol(dd)] = apply(dd[ , 6:ncol(dd)], 2, scale)

# check for multicollinearity
usdm::vifstep(dd[ , 6:ncol(dd) ])
corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::filter(complete.cases(.))))

# specify which multicollinear vars to exclude in later analyses if any
mc_vars_to_exclude = c()

# examine VIFs/corplots for multivar and causal var subsets
corrplot::corrplot( cor( dd[ , 6:ncol(dd) ] %>% dplyr::select(-all_of(mc_vars_to_exclude)) %>% dplyr::filter(complete.cases(.) )))
usdm::vifstep( dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_strict)) %>% dplyr::filter(complete.cases(.)) )

# predictors
predictors = names(dd)[ 6:(ncol(dd)-1)]
preds_m = predictors[ !predictors %in% mc_vars_to_exclude] 
preds_m = preds_m[ preds_m %in% hyp_lib ]
preds_c = preds_m[ preds_m %in% hyp_strict ]







# =============== mgcv test; run these for SI ================

library(mgcv)

# response variable: 1/0 outbreaks
dd$y = dd$presence
family1 = "binomial"

formx = paste(c("y ~ 1 + s(Longitude, Latitude, k=29)", paste("s(", preds_m, ", k=6)", sep="")), collapse=" + ")

m1 = mgcv::gam(formula(formx), data=dd, family=binomial(link="logit"))

# draw the 2D tensor smooth
gratia::draw(m1, 
             select=1, 
             continuous_fill = scale_fill_distiller(palette = "Spectral", type = "div"))

# draw the splines
gratia::draw(m1, select=-1) & theme_classic()



