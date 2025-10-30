

# ================ MARBURG (CENTRAL & WEST AFRICA) ===================


# ------------ housekeeping ----------------

# setup objects
# setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
# library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
# library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# defaults to extracting covars; can be switched to FALSE for rerunning models only
if(!exists("generate_dataset")){
  generate_dataset = TRUE
}


# ---------------- Disease data ---------------

# random seed and disease name
rseed = 4582
disease_name = "marburg"

# disease data and shapefile post-1980
dz = read.csv("./output/spillovers_processed/spillovers_marv.csv") %>%
  dplyr::filter(Year >= 1980) %>%
  dplyr::mutate(record_id = 1:length(Disease)) %>% 
  dplyr::mutate(dup = duplicated(paste(LocationName, Year))) %>%
  dplyr::filter(dup == FALSE) %>%
  dplyr::select(-dup)
load("./output/spillovers_processed/spillovers_marv_shp.R")
dz_shp = shp %>% dplyr::filter(ADMcode %in% dz$ADMcode); rm(shp)

# create study area polygon
study_area = createStudyAreaPolygon(lat=dz$Latitude, lon=dz$Longitude)


# ------------ Analysis variables ------------

# define important analysis variables for specific disease
occurrence_point_buffer_km = 10 # km buffer around each occurrence point for spatial uncertainty in infection event
filter_area_thresh = 10000 # maximum geographical area of occurrence polygons to include in km2 (i.e. exclude v imprecise records)
num_background_points = nrow(dz) * 10 # number of background points to generate
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

# hypotheses weakest 
hyp_weak = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv") %>%
  dplyr::filter(covar_name != "") %>%
  dplyr::filter(dz_name == dz$Disease[1]) %>%
  dplyr::filter(TestAny == "Test") %>%
  dplyr::select(covar_name)
hyp_weak = as.vector(hyp_weak$covar_name)

# add a priori confounds
if("forest_loss" %in% hyp_lib ){ hyp_lib = c(hyp_lib, "forest_cover") }
if("forest_loss" %in% hyp_strict ){ hyp_strict = c(hyp_strict, "forest_cover") }
if("forest_loss" %in% hyp_weak ){ hyp_weak = c(hyp_weak, "forest_cover") }
hyp_weak = unique(c(hyp_weak, confound)); hyp_weak = hyp_weak[ hyp_weak != "" ]
hyp_lib = unique(c(hyp_lib, confound)); hyp_lib = hyp_lib[ hyp_lib != "" ]
hyp_strict = unique(c(hyp_strict, confound)); hyp_strict = hyp_strict[ hyp_strict != "" ]





# ----------- 1. Build objects for modelling (dependencies: 00_analysis_funcs.R) -------------

if(generate_dataset == TRUE){
  
  # harmonise points and polygons
  # creating X km radius buffer around all point locations, and filtering out any areas >threshold size
  sf::sf_use_s2(FALSE)
  dz = harmoniseDataTypes(dz, dz_shp, point_buffer_km=occurrence_point_buffer_km, filter_area_thresh=filter_area_thresh)
  sf::sf_use_s2(TRUE)
  
  # background points: population weighted point+buffers with same median polygon area
  med_area_size = median(dz$poly_area)
  buf_rad = sqrt(med_area_size / pi) * 1000 # calculates radius in km
  bg = generateBackgroundPolygons(num_points = num_background_points,
                                  study_area = study_area,
                                  buffer_radius = buf_rad,
                                  method = "popweight",
                                  seed = rseed)
  
  # combine into dataset
  dd = rbind(
    dz %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 1),
    bg %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 0)
  )
  
  # extract covariates
  print("Extracting covariate data")
  covs_df = read.csv("./scripts/03_modelling/covariate_lookup.csv")
  ddx = dd
  
  for(cc in covs_df$cov){
    print(cc)
    if(cc == "forest_loss"){ sf::sf_use_s2(FALSE) }
    extr_x = extractCovariateVals(dd, covariate_name=cc) # function stored in '00_covar_extraction_funcs.R'
    sf::sf_use_s2(TRUE)
    ddx = cbind(ddx, extr_x)
  }
  
  # log transform any relevant variables, drop geometry and save
  ddx = ddx %>%
    dplyr::mutate(health_travel_log = log(health_travel+1),
                  livestock_log = log(livestock_all+1))
  filename = paste("./output/model_df/" , disease_name, "_df_wp.csv", sep="")
  write.csv(sf::st_drop_geometry(ddx), filename, row.names=FALSE)
  
  
  
  # # ----------- 2. Initial visualisation: check everything has worked ------------
  # 
  # # p1: study area + polygons
  # p1 = ggplot() +
  #   geom_sf(data=study_area, fill="grey95", col="grey80") +
  #   maptheme +
  #   geom_sf(data=dz, col="red", fill=NA) + ggtitle(paste(dz$Disease[1], ", n=", nrow(dz), sep=""))
  # 
  # # p1 = study area + background
  # p2 = ggplot() +
  #   geom_sf(data=study_area, fill="grey95", col="grey80") +
  #   maptheme +
  #   geom_point(data=ddx[ ddx$presence==0,], aes(Longitude, Latitude), size=0.25, col="black") +
  #   geom_point(data=ddx[ ddx$presence==1,], aes(Longitude, Latitude), size=1, col="red")
  # 
  # # extracted covariates
  # p3 = ddx %>%
  #   reshape2::melt(id.vars = 1:5) %>%
  #   dplyr::filter(!variable %in% c("livestock_all", "health_travel")) %>%
  #   ggplot() +
  #   ggforce::geom_sina(aes(x = factor(presence), y=value, group=factor(presence)), width=1, size=0.1, col="grey60") +
  #   geom_boxplot(aes(x = factor(presence), y=value, group=factor(presence), fill=variable), width=0.25, alpha=0.8) +
  #   facet_wrap(~variable, scales="free_y") +
  #   theme_classic() +
  #   xlab("Presence") +
  #   ylab("Covariate value") +
  #   theme(legend.position="none",
  #         strip.background = element_blank(), strip.text = element_text(size=13))
  # 
  # pc1 = gridExtra::grid.arrange(p1, p2, nrow=2)
  # pc = gridExtra::grid.arrange(pc1, p3, nrow=1)
  
}


# ------------- 3. INLA spatial regression model ---------------

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
# usdm::vifstep(dd[ , 6:ncol(dd) ])
# corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_lib)) %>% dplyr::filter(complete.cases(.))))
# corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_strict)) %>% dplyr::filter(complete.cases(.))))

# specify which multicollinear vars to exclude in later analyses if any
mc_vars_to_exclude = c("urban_expansion", "crop_expansion", "hunting")

# examine VIFs/corplots for multivar and causal var subsets
# corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(-all_of(mc_vars_to_exclude)) %>% dplyr::filter(complete.cases(.))))
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
max.edge = 1.2
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
# prior median range at approximately 0.5*diff(range(spatial_locs$Longitude)) - roughly extent of study area
# prior probability of marginal sigma of 3 or more is 0.01
spde = inla.spde2.pcmatern(mesh1, 
                           prior.range = c(15, 0.1), 
                           prior.sigma = c(10, 0.01))

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
m0 = fitINLAModel(formx = formula(form_base), family=family1, stack=stack1, spde=spde)

# model 1: univariate models including each covariate individually
predictors = names(dd)[ 6:(ncol(dd)-1)]
m1 = vector("list", length=length(predictors))
#fm1 = data.frame()

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


# remove multicollinear vars
preds_m = predictors[ !predictors %in% mc_vars_to_exclude] 

# model 2: majority rule
preds_l = preds_m[ preds_m %in% hyp_lib ]
form = formula(paste(c(form_base, preds_l), collapse=" + "))
m2 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde, config=TRUE)

# model 3: top ranked
preds_c = preds_m[ preds_m %in% hyp_strict ]
form = formula(paste(c(form_base, preds_c), collapse=" + "))
m3 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde, config=TRUE)

# model 4: any author
preds_w = preds_m[ preds_m %in% hyp_weak ]
form = formula(paste(c(form_base, preds_w), collapse=" + "))
m4 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde, config=TRUE)
