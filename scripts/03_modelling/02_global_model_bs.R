

# ================ Runs global model (disease agnostic) with bootstrap to even sampling across diseases  ===================

# This script fits the global models to infer drivers of outbreak risk (including all diseases)
# To ensure even representation of diseases with different quantities of data, this applies 
# a bootstrap approach, subsampling to 100 points per disease
# submodels are fitted 100 times, posterior samples are extracted and these are used to build 
# up the combined posterior distribution across all submodels
# The results are visualised in the script "04_results/02_global_model.R"

# ------------ housekeeping ----------------

# setup objects
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")



# ---------------- 1. compile disease data ---------------

# reads extracted covariates and geo data for all diseases
# individual datasets generated in each individual disease script ("./scripts/03_modelling/inla_per_disease")

ll = list.files("./output/model_df/", pattern=".csv", full.names=TRUE)

ll = ll[ -grep("brazil|usa|filo", ll)]

rr = read.csv(ll[1])
nn = names(rr)
nn = nn[ - which(nn %in% c("livestock_all", "tmean", "tmean_anomaly", "precip", "precip_anomaly", "evi_coefvar", "health_travel")) ]

dd = do.call(rbind.data.frame, lapply(ll, function(x) read.csv(x) %>% dplyr::select(all_of(nn)) %>% dplyr::mutate(file=x)))

dd$Disease = unlist(
  lapply(
    strsplit(
      unlist(lapply(strsplit(dd$file, "/"), "[", 4)), 
      "_"
    ),
    "[", 1
  )
)

dd = cbind(dd[ , (ncol(dd)-1) : ncol(dd), drop=FALSE], dd[ , 1:(ncol(dd)-2)])

dd = dd %>% dplyr::select(-ADMcode, -record_id)

dd = dd %>% dplyr::filter(presence == 1)
dd = dd %>% dplyr::filter(Longitude > -130)
dd = dd %>% dplyr::filter(Latitude < 50)

# create study area polygon with country/coastal boundaries
study_area = createStudyAreaPolygon(lat=dd$Latitude, lon=dd$Longitude, buffer=FALSE, extent_border = 1)



# ------------ analysis variables ------------

# define important analysis variables for specific disease
occurrence_point_buffer_km = 5 # km buffer around each occurrence point for spatial uncertainty in infection event
num_background_points = 50000 # number of background points to generate
covar_coverage_thresh = 0.95 # only keep covariates that are non-NA for at least X% of points
zeroes_coverage_thresh = 0.95 # only keep covariates that are non-zero for at least X% of points

# seed for background point generation
#rseed = 61947


# --------- 2. extract covariates at background points across globe -----------

## THIS IS RUN ONCE, THEN COMMENTED OUT AS TIME CONSUMING TO RUN

# background points: population weighted point+buffers with same median polygon area
# buf_rad = 10000
# bg = generateBackgroundPolygons(num_points = num_background_points,
#                                 study_area = study_area,
#                                 buffer_radius = buf_rad,
#                                 method = "popweight",
#                                 seed = rseed,
#                                 pop_res = "10km")
#
# bg = bg %>% dplyr::select(record_id, Longitude, Latitude, ADMcode) %>% dplyr::mutate(presence = 0)
# sf::st_write(bg, "C:/Users/roryj/Dropbox/Research/fingerprint/global_background/global_background_points.shp")

#bg = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/global_background/global_background_points.shp")

# bg %>%
#   st_drop_geometry() %>%
#   ggplot() + 
#   geom_point(aes(Longitude, Latitude), size=0.1) + 
#   coord_fixed()

# # extract covariates
# print("Extracting covariate data")
# covs_df = read.csv("./scripts/03_modelling/covariate_lookup.csv")
# covs_df = covs_df %>% dplyr::filter(!cov %in% c("precip_anomaly", "tmean_anomaly"))
# ddx = bg
# disease_name = "global_background"
# for(cc in covs_df$cov){
#   print(cc)
#   extr_x = extractCovariateVals(bg, covariate_name=cc) # function stored in '00_covar_extraction_funcs.R'
#   write.csv(extr_x, paste("./output/model_df/global_background/bg_", cc, ".csv", sep=""), row.names=FALSE)
#   ddx = cbind(ddx, extr_x)
# }
# 
# # log transform any relevant variables, drop geometry and save
# ddx = ddx %>%
#   dplyr::mutate(health_travel_log = log(health_travel+1),
#                 livestock_log = log(livestock_all+1))
# filename = paste("./output/model_df/global_background/globalbackgroundpoints_df_wp.csv", sep="")
# write.csv(sf::st_drop_geometry(ddx), filename, row.names=FALSE)



# --------- 3. INLA geospatial regression models ---------------


# ------- read and prepare data ------

# background points
bg = read.csv("./output/model_df/global_background/globalbackgroundpoints_df_wp.csv")
bg = bg[ , c(2:3, 5:ncol(bg))] %>%
  dplyr::mutate(Disease = "Background") %>%
  dplyr::select(-all_of(c("tmean", "precip", "health_travel", "evi_coefvar", "livestock_all")))

dd = dd %>% dplyr::select(-file)

# combine with disease outbreaks
dd = rbind(dd, bg)

# proportion NAs in each predictor: exclude predictors that do not have coverage for at least X% of data
nas = apply(dd[ dd$presence == 1, 5:ncol(dd)], 2, function(x){ sum(!is.na(x))/length(x) })
nas = names(nas[ nas < covar_coverage_thresh ])

# proportion 0s in each predictor: exclude predictors with zeroes for more than X% of data
zeroes = apply(dd[ , 5:ncol(dd)], 2, function(x){ sum(!is.na(x) & x == 0)/length(x) })
zeroes = names(zeroes[ zeroes > zeroes_coverage_thresh ])

# excluding zeroes/nas
dd = dd %>% dplyr::select(-all_of(nas), -all_of(zeroes))

# scale predictors
dd[ , 5:ncol(dd)] = apply(dd[ , 5:ncol(dd)], 2, scale)

# check for multicollinearity
# usdm::vifstep(dd[ , 5:ncol(dd) ])
# corrplot::corrplot(cor(dd[ , 5:ncol(dd) ] %>% dplyr::filter(complete.cases(.))))

# specify which multicollinear vars to exclude in later analyses if any
mc_vars_to_exclude = c("urban_expansion")
# 
# # examine VIFs/corplots for multivar and causal var subsets
# corrplot::corrplot( cor( dd[ , 6:ncol(dd) ] %>% dplyr::filter(complete.cases(.) )))
# usdm::vifstep( dd[ , 6:ncol(dd) ] %>% dplyr::select(all_of(hyp_strict)) %>% dplyr::filter(complete.cases(.)) )

# predictors
predictors = names(dd)[ 5:(ncol(dd))]
preds_m = predictors[ !predictors %in% mc_vars_to_exclude] 



# ------ setup objects for submodels by principal human exposure route --------

dd$transmission = "ZVB"
dd$transmission[ dd$Disease %in% c("melioidosis") ] = "ENV"
dd$transmission[ dd$Disease %in% c("zika", "dengue", "chikv") ] = "VB"
dd$transmission[ dd$Disease %in% c("anthrax", "ebola", "h5n1", "hantavirus", "junin", "lassa", "mers", "mpx", "nipah", "marburg") ] = "Z"
dd$transmission[ dd$Disease == "Background" ] = "BK"

# study areas for specific transmission groups
zvb_studyarea = createStudyAreaPolygon(lat=dd$Latitude[ dd$transmission == "ZVB"], lon=dd$Longitude[ dd$transmission == "ZVB"], buffer=FALSE, extent_border = 1)
vb_studyarea = createStudyAreaPolygon(lat=dd$Latitude[ dd$transmission == "VB"], lon=dd$Longitude[ dd$transmission == "VB"], buffer=FALSE, extent_border = 1)

# identify intersecting background points 

# zoonotic and vector borne
bb = st_as_sf(dd, coords = c("Longitude", "Latitude")); st_crs(bb) = crs(study_area) 
sf_use_s2(FALSE)
ii = sf::st_intersects(bb, zvb_studyarea)
ii = unlist(lapply(ii, function(x) ifelse(length(x)==0, 0, x)))
dd$ZVB = ii
dd$ZVB[ dd$transmission %in% c("ENV", "Z", "VB") ] = 0

# # zoonotic (direct)
ii = sf::st_intersects(bb, zvb_studyarea)
ii = unlist(lapply(ii, function(x) ifelse(length(x)==0, 0, x)))
dd$Z = ii
dd$Z[ dd$transmission %in% c("ENV", "ZVB", "VB") ] = 0

# vector borne
ii = sf::st_intersects(bb, vb_studyarea)
ii = unlist(lapply(ii, function(x) ifelse(length(x)==0, 0, x)))
dd$VB = ii
dd$VB[ dd$transmission %in% c("ENV", "ZVB", "Z") ] = 0

# ggplot() +
#   geom_sf(data=zvb_studyarea, fill="grey95", color=NA) +
#   geom_point(data=dd[ (dd$Z == 1 | dd$ZVB == 1) & dd$Disease != "Background", ], aes(Longitude, Latitude, color=factor(transmission)), size=0.2, alpha=0.4) +
#   maptheme

save(zvb_studyarea, file="./output/model_outputs/global_model_bs/zvb_studyarea.R")


# ------ model objects and function --------

# likelihood and key priors
family1 = "binomial"
hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                      prec.intercept=0.1, # prior precision for intercept
                      mean=0, # prior mean for fixed effects
                      prec=1)  # prior precision for fixed effects

# build mesh using spatial locs
# n.b. built outside the function so standardised across all bootstrap iterations
spatial_locs = dd[ , c("Longitude", "Latitude")]
max.edge = 2.5
bound.outer = diff(range(spatial_locs$Longitude))/10 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area), 
                     loc.domain = spatial_locs,
                     max.edge = c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))

# function to fit global multidisease model (modifiable for subsets)

fitGlobalModel = function(data, # model dataframe
                          study_shp, #study_area shapefile
                          predictors, # which covariates to include in the model
                          type="Global", # giving a name to label the model
                          mesh_edge=2.5, # max edge parameter for mesh (controls how fine)
                          prange=20, # prior range for SPDE
                          psigma=10, # prior sigma for SPDE
                          fit_spde=TRUE,
                          verbose = FALSE){ # whether to include geospatial random effect
  
  # response variable: 1/0 outbreaks
  data$y = data$presence
  
  # create projector matrix for point locations based on mesh
  sl = data[ , c("Longitude", "Latitude")]
  A = inla.spde.make.A(mesh1, loc=as.matrix(sl))
  
  # define pc matern spde model
  spde = inla.spde2.pcmatern(mesh1, 
                             prior.range = c(prange, 0.1), 
                             prior.sigma = c(psigma, 0.01))
  
  # build stack 
  stack1 = inla.stack(data=list(y=data$y),
                      effects=list(s = 1:mesh1$n, # spatial
                                   data.frame(Intercept=1, data[ , -which(names(data) %in% c("y"))]),
                                   iidx = 1:nrow(data)),
                      A=list(A, 1, 1),
                      remove.unused = FALSE, tag='est')
  
  # baseline
  if(fit_spde == TRUE){
    form_base = "y ~ -1 + Intercept + f(s, model=spde)"
  } else{
    form_base = "y ~ -1 + Intercept"
  }

  # with predictors
  form = paste(c(form_base, predictors), collapse=" + ")
  m1 = fitINLAModel(formx = formula(form), family=family1, stack=stack1, spde=spde, verbose = verbose, config = TRUE)
  
  fx = extractFixedINLA(m1) %>% 
    dplyr::mutate(type = type, paramname = param) %>%
    dplyr::mutate(param = replace(param, param == "health_travel_log", "Travel time to health facility (log)"),
                  param = replace(param, param == "livestock_log", "Livestock density (log)"),
                  param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
                  param = replace(param, param == "forest_cover", "Forest cover"),
                  param = replace(param, param == "forest_loss", "Forest loss"),
                  param = replace(param, param == "crop_expansion", "Cropland expansion"),
                  param = replace(param, param == "evi_dissimilarity", "EVI dissimilarity"),
                  param = replace(param, param == "evi_coefvar", "EVI variance"),
                  param = replace(param, param == "precip_anomaly", "Precip anomaly"),
                  param = replace(param, param == "urban_expansion", "Urban expansion"),
                  param = replace(param, param == "urban_cover", "Urban cover"),
                  param = replace(param, param == "hunting", "Hunting pressure"),
                  param = replace(param, param == "crop_cover", "Cropland cover"),
                  param = replace(param, param == "mining", "Mining coverage"),
                  param = replace(param, param == "tmean_change", "Tmean change"),
                  param = replace(param, param == "precip_change", "Precip change"),
                  param = replace(param, param == "social_vulnerability", "Social vulnerability"),
                  param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
                  param = replace(param, param == "protected_areas", "Protected area coverage"),
                  param = replace(param, param == "popdens_log", "Population density (log)"))
  
  fx %>% 
    dplyr::filter(paramname != "Intercept") %>%
    ggplot() + 
    geom_hline(yintercept=0, lty=2) +
    geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=2) + 
    geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.6, alpha=0.7) +
    theme_classic() + 
    ylab("Posterior marginal slope (log odds)") + 
    xlab("") + 
    theme(legend.position = "none",
          plot.title = element_text(size=16, hjust=0.5),
          axis.text = element_text(size=12), 
          axis.title = element_text(size=13),
          legend.title = element_blank(),
          legend.text = element_text(size=13), 
          strip.background = element_blank(),
          strip.text = element_text(size=14)) +
    scale_color_viridis_d(begin=0, end=0.8, direction=-1) + 
    coord_flip() +
    facet_wrap(~type) 
  
  if(fit_spde){
    
    rasteriseSPDE = function(fitted_mesh_x, mesh_x, shp_x){
      
      proj = inla.mesh.projector(mesh_x, dims=c(2000, 2000))
      full.proj = expand.grid(x = proj$x, y=proj$y)
      full.proj$field = c(inla.mesh.project(proj, fitted_mesh_x))
      pras = rasterFromXYZ(full.proj)
      mask = fasterize::fasterize(shp_x, pras)
      pras = mask(pras, mask)
      pras = crop(pras, shp_x)
      return(pras)
    }
    
    spde_ras = rasteriseSPDE(m1$summary.random$s$mean, mesh1, study_shp)
    
    rasteriseSPDE(m1$summary.random$s$mean, mesh1, study_shp) %>%
      as.data.frame(xy=TRUE) %>%
      ggplot() + 
      geom_raster(aes(x, y, fill=field)) +
      scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
      maptheme +
      geom_sf(data=study_shp, fill=NA, col="black") + 
      geom_point(data=data[ data$presence == 1, ], aes(Longitude, Latitude), col="red", size=0.2, alpha=0.2) + 
      theme(legend.position = c(0.98, 0.15), plot.title = element_text(hjust=0.5, size=16))
  } else{
    spde_ras = "none"
  }
  
  m_obj = list(
    type = type,
    shp = study_shp,
    spde = spde_ras,
    fx = fx,
    mesh_edge = mesh_edge,
    spde_priorrange = prange,
    spde_priorsigma = psigma, 
    model = m1
  )
  return(m_obj)
  
}


# ===================== 1. Overall model across all diseases ===========================

# fit models: compare inference between models with the following attributes:
# 1. non-spatial non-detection
# 2. spatial but with no detection variables
# 3. geospatial with detection variables

# modelling has the issue of variable sampling effort/data quantity across diseases
# to address this, fit the model to the same number of points for each disease (100)
# fit N iterations of the model, each fitted to randomly drawn subsamples of 100 per disease
# draw posterior samples of fixed effects, recover SPDEs, and save each iteration
# then combine these outputs later

# ---- 1. Non spatial model with no detection covariates

for(i in 1:100){
  
  print(paste("fitting model", i))
  ts = Sys.time()
  
  # subset to 100 points per disease
  npts = dd %>% dplyr::group_by(Disease) %>% dplyr::summarise(npts = length(Disease))
  dd_s = dd %>% 
    dplyr::filter(Disease != "Background") %>%
    dplyr::filter(Disease %in% npts$Disease[ npts$npts >= 100 ]) %>%
    dplyr::group_by(Disease) %>%
    dplyr::sample_n(100, replace=FALSE) %>%
    rbind(
      dd %>% dplyr::filter(Disease %in% npts$Disease[ npts$npts < 100 ])
    )
  dd_s = dd_s %>%
    rbind(
      dd %>% dplyr::filter(Disease == "Background") %>%
        dplyr::sample_n(nrow(dd_s)*2, replace=FALSE)
    )
  
  # fit model
  m_s = fitGlobalModel(data=dd_s, study_shp=study_area, 
                       predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
                       type="Nonspatial_nondetection", 
                       mesh_edge=2.5, prange=20, psigma=10, fit_spde = FALSE)
  
  # draw posterior samples
  samp = INLA::inla.posterior.sample(n=1000, m_s$model)
  
  # extract fixefs from samples
  samp_extr = function(x){
    s = samp[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(names(dd_s), ":1", sep="")), ]
    t(as.data.frame(s))
  }
  fx_samp = do.call(rbind.data.frame, lapply(1:length(samp), samp_extr)) %>%
    dplyr::mutate(type = "Nonspatial_nondetection", id=sample(1:10^6, 1))
  
  # extract spde
  #spde = m_s$spde
  
  # save
  write.csv(fx_samp, paste("./output/model_outputs/global_model_bs/outputs/all_diseases/", fx_samp$type[1], "_", fx_samp$id[1], "_fixefs.csv", sep=""), row.names=FALSE)
  #save(spde, file=paste("./output/model_outputs/global_model_bs/outputs/all_diseases/", fx_samp$type[1], "_", fx_samp$id[1], "_spde.R", sep=""))
  
  te = Sys.time()
  print(te - ts)
}


# ---- 2. Spatial model with no detection covariates

for(i in 1:100){
  
  print(paste("fitting model", i))
  ts = Sys.time()
  
  # subset to 100 points per disease
  npts = dd %>% dplyr::group_by(Disease) %>% dplyr::summarise(npts = length(Disease))
  dd_s = dd %>% 
    dplyr::filter(Disease != "Background") %>%
    dplyr::filter(Disease %in% npts$Disease[ npts$npts >= 100 ]) %>%
    dplyr::group_by(Disease) %>%
    dplyr::sample_n(100, replace=FALSE) %>%
    rbind(
      dd %>% dplyr::filter(Disease %in% npts$Disease[ npts$npts < 100 ])
    )
  dd_s = dd_s %>%
    rbind(
      dd %>% dplyr::filter(Disease == "Background") %>%
        dplyr::sample_n(nrow(dd_s)*2, replace=FALSE)
    )
  
  # fit model
  m_s = fitGlobalModel(data=dd_s, study_shp=study_area, 
                       predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
                       type="Spatial_nondetection", 
                       mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
  
  # draw posterior samples
  samp = INLA::inla.posterior.sample(n=1000, m_s$model)
  
  # extract fixefs from samples
  samp_extr = function(x){
    s = samp[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(names(dd_s), ":1", sep="")), ]
    t(as.data.frame(s))
  }
  fx_samp = do.call(rbind.data.frame, lapply(1:length(samp), samp_extr)) %>%
    dplyr::mutate(type = "Spatial_nondetection", id=sample(1:10^6, 1))
  
  # extract spde
  spde = m_s$spde
  
  # save
  write.csv(fx_samp, paste("./output/model_outputs/global_model_bs/outputs/all_diseases/", fx_samp$type[1], "_", fx_samp$id[1], "_fixefs.csv", sep=""), row.names=FALSE)
  save(spde, file=paste("./output/model_outputs/global_model_bs/outputs/all_diseases/", fx_samp$type[1], "_", fx_samp$id[1], "_spde.R", sep=""))
  
  te = Sys.time()
  print(te - ts)
}


# --- 3. Spatial model with detection covariates

for(i in 1:50){
  
  print(paste("fitting model", i))
  ts = Sys.time()
  
  # subset to 100 points per disease
  npts = dd %>% dplyr::group_by(Disease) %>% dplyr::summarise(npts = length(Disease))
  dd_s = dd %>% 
    dplyr::filter(Disease != "Background") %>%
    dplyr::filter(Disease %in% npts$Disease[ npts$npts >= 100 ]) %>%
    dplyr::group_by(Disease) %>%
    dplyr::sample_n(100, replace=FALSE) %>%
    rbind(
      dd %>% dplyr::filter(Disease %in% npts$Disease[ npts$npts < 100 ])
    )
  dd_s = dd_s %>%
    rbind(
      dd %>% dplyr::filter(Disease == "Background") %>%
        dplyr::sample_n(nrow(dd_s)*2, replace=FALSE)
    )
  
  # fit model
  m_s = fitGlobalModel(data=dd_s, study_shp=study_area, 
                       predictors=preds_m,
                       type="Spatial_detection", 
                       mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE,
                       verbose = FALSE)
  
  # draw posterior samples
  samp = INLA::inla.posterior.sample(n=1000, m_s$model)
  
  # extract fixefs from samples
  samp_extr = function(x){
    s = samp[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(names(dd_s), ":1", sep="")), ]
    t(as.data.frame(s))
  }
  fx_samp = do.call(rbind.data.frame, lapply(1:length(samp), samp_extr)) %>%
    dplyr::mutate(type = "Spatial_detection", id=sample(1:10^6, 1))
  
  # extract spde
  spde = m_s$spde
  
  # save
  write.csv(fx_samp, paste("./output/model_outputs/global_model_bs/outputs/all_diseases/", fx_samp$type[1], "_", fx_samp$id[1], "_fixefs.csv", sep=""), row.names=FALSE)
  save(spde, file=paste("./output/model_outputs/global_model_bs/outputs/all_diseases/", fx_samp$type[1], "_", fx_samp$id[1], "_spde.R", sep=""))
  
  te = Sys.time()
  print(te - ts)
}




# ====================== breaking down models by key transmission mode ============================

# build mesh
spatial_locs = dd[ dd$ZVB == 1 | dd$Z == 1 | dd$VB == 1 , c("Longitude", "Latitude")]
max.edge = 2.5
bound.outer = diff(range(spatial_locs$Longitude))/10 # outer bound: 1/3 of range diff initially
mesh1 = inla.mesh.2d(boundary = as_Spatial(study_area),
                     loc.domain = spatial_locs,
                     max.edge = c(1, 5)*max.edge,
                     cutoff = max.edge/5,
                     offset = c(max.edge, bound.outer))

# -------- 1. zoonotic (any route) ----------- 

for(i in 1:100){
  
  print(paste("fitting model", i))
  ts = Sys.time()
  
  # subset to 100 points per disease
  dd_s = dd[ dd$ZVB == 1 | dd$Z == 1, ]
  npts = dd_s %>% dplyr::group_by(Disease) %>% dplyr::summarise(npts = length(Disease))
  dd_s = dd_s %>% 
    dplyr::filter(Disease != "Background") %>%
    dplyr::filter(Disease %in% npts$Disease[ npts$npts >= 100 ]) %>%
    dplyr::group_by(Disease) %>%
    dplyr::sample_n(100, replace=FALSE) %>%
    rbind(
      dd %>% dplyr::filter(Disease %in% npts$Disease[ npts$npts < 100 ])
    )
  dd_s = dd_s %>%
    rbind(
      dd %>% dplyr::filter(Disease == "Background") %>%
        dplyr::sample_n(nrow(dd_s)*2, replace=FALSE)
    )
  
  # fit model
  m_s = fitGlobalModel(data=dd_s, 
                       study_shp=study_area, 
                       predictors=preds_m,
                       type="Spatial_detection", 
                       mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE,
                       verbose=FALSE)
  
  # draw posterior samples
  samp = INLA::inla.posterior.sample(n=1000, m_s$model)
  
  # extract fixefs from samples
  samp_extr = function(x){
    s = samp[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(names(dd_s), ":1", sep="")), ]
    t(as.data.frame(s))
  }
  fx_samp = do.call(rbind.data.frame, lapply(1:length(samp), samp_extr)) %>%
    dplyr::mutate(type = "Spatial_detection", id=sample(1:10^6, 1))
  
  # extract spde
  spde = m_s$spde
  
  # save
  write.csv(fx_samp, paste("./output/model_outputs/global_model_bs/outputs/zoonotic/", fx_samp$type[1], "_", fx_samp$id[1], "_fixefs.csv", sep=""), row.names=FALSE)
  save(spde, file=paste("./output/model_outputs/global_model_bs/outputs/zoonotic/", fx_samp$type[1], "_", fx_samp$id[1], "_spde.R", sep=""))
  
  te = Sys.time()
  print(te - ts)
}

# -------- 2. vector-borne (any host) ----------- 

for(i in 1:100){
  
  print(paste("fitting model", i))
  ts = Sys.time()
  
  # subset to 100 points per disease
  dd_s = dd[ dd$ZVB == 1 | dd$VB == 1, ]
  npts = dd_s %>% dplyr::group_by(Disease) %>% dplyr::summarise(npts = length(Disease))
  dd_s = dd_s %>% 
    dplyr::filter(Disease != "Background") %>%
    dplyr::filter(Disease %in% npts$Disease[ npts$npts >= 100 ]) %>%
    dplyr::group_by(Disease) %>%
    dplyr::sample_n(100, replace=FALSE) %>%
    rbind(
      dd %>% dplyr::filter(Disease %in% npts$Disease[ npts$npts < 100 ])
    )
  dd_s = dd_s %>%
    rbind(
      dd %>% dplyr::filter(Disease == "Background") %>%
        dplyr::sample_n(nrow(dd_s)*2, replace=FALSE)
    )
  
  # fit model
  m_s = fitGlobalModel(data=dd_s, 
                       study_shp=zvb_studyarea, 
                       predictors=preds_m,
                       type="Spatial_detection", 
                       mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE,
                       verbose=FALSE)
  
  # draw posterior samples
  samp = INLA::inla.posterior.sample(n=1000, m_s$model)
  
  # extract fixefs from samples
  samp_extr = function(x){
    s = samp[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(names(dd_s), ":1", sep="")), ]
    t(as.data.frame(s))
  }
  fx_samp = do.call(rbind.data.frame, lapply(1:length(samp), samp_extr)) %>%
    dplyr::mutate(type = "Spatial_detection", id=sample(1:10^6, 1))
  
  # extract spde
  spde = m_s$spde
  
  # save
  write.csv(fx_samp, paste("./output/model_outputs/global_model_bs/outputs/vectorborne/", fx_samp$type[1], "_", fx_samp$id[1], "_fixefs.csv", sep=""), row.names=FALSE)
  save(spde, file=paste("./output/model_outputs/global_model_bs/outputs/vectorborne/", fx_samp$type[1], "_", fx_samp$id[1], "_spde.R", sep=""))
  
  te = Sys.time()
  print(te - ts)
}








# # ===================== 1. Overall model across all diseases ===========================
# 
# # fit models: compare inference between models with the following attributes:
# # 1. non-spatial non-detection
# # 2. spatial but with no detection variables
# # 3. geospatial with detection variables
# 
# m_nonspatial1 = fitGlobalModel(data=dd, study_shp=study_area, 
#                                predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
#                                type="Nonspatial_nondetection", 
#                                mesh_edge=2.5, prange=20, psigma=10, fit_spde = FALSE)
# save(m_nonspatial1, file="./output/model_outputs/global_models/global_m1_nonspatial.R")
# 
# m_spatial1 = fitGlobalModel(data=dd, study_shp=study_area, 
#                                predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
#                                type="Spatial_nondetection", 
#                                mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
# save(m_spatial1, file="./output/model_outputs/global_models/global_m2_spatial.R")
# 
# m_spatial2 = fitGlobalModel(data=dd, study_shp=study_area, 
#                             predictors=preds_m,
#                             type="Spatial_detection", 
#                             mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
# save(m_spatial2, file="./output/model_outputs/global_models/global_m3_spatial.R")
# 
# 
# # ===================== 2. Only for zoonotic reservoir ===========================
# 
# # Same approach as above, but only for diseases where transmission is from a wildlife reservoir i.e. explicitly ecosystem-linked
# # Excludes predominantly environmental and anthroponotic transmission
# 
# mz_nonspatial1 = fitGlobalModel(data=dd[ dd$ZVB == 1 | dd$Z == 1, ], 
#                                 study_shp=zvb_studyarea, 
#                                 predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
#                                 type="Nonspatial_nondetection", 
#                                 mesh_edge=2.5, prange=20, psigma=10, fit_spde = FALSE)
# save(mz_nonspatial1, file="./output/model_outputs/global_models/global_m1_nonspatial_zzvb.R")
# 
# mz_spatial1 = fitGlobalModel(data=dd[ dd$ZVB == 1 | dd$Z == 1, ], 
#                                 study_shp=zvb_studyarea, 
#                                 predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
#                                 type="Spatial_nondetection", 
#                                 mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
# save(mz_spatial1, file="./output/model_outputs/global_models/global_m2_spatial_zzvb.R")
# 
# mz_spatial2 = fitGlobalModel(data=dd[ dd$ZVB == 1 | dd$Z == 1, ], 
#                             study_shp=zvb_studyarea, 
#                             predictors=preds_m,
#                             type="Spatial_detection", 
#                             mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
# save(mz_spatial2, file="./output/model_outputs/global_models/global_m3_spatial_zzvb.R")
# 

# # ===================== 3. Only for vector-borne transmission ===========================
# 
# # includes diseases with any vectors involved in transmission
# 
# mz_nonspatial1 = fitGlobalModel(data=dd[ dd$ZVB == 1 | dd$VB == 1, ], 
#                                 study_shp=zvb_studyarea, 
#                                 predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
#                                 type="Nonspatial_nondetection", 
#                                 mesh_edge=2.5, prange=20, psigma=10, fit_spde = FALSE)
# save(mz_nonspatial1, file="./output/model_outputs/global_models/global_m1_nonspatial_vb.R")
# 
# mz_spatial1 = fitGlobalModel(data=dd[ dd$ZVB == 1 | dd$VB == 1, ], 
#                              study_shp=zvb_studyarea, 
#                              predictors=preds_m[ !preds_m %in% c("health_travel_log", "urban_cover")],
#                              type="Spatial_nondetection", 
#                              mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
# save(mz_spatial1, file="./output/model_outputs/global_models/global_m2_spatial_vb.R")
# 
# mz_spatial2 = fitGlobalModel(data=dd[ dd$ZVB == 1 | dd$VB == 1, ], 
#                              study_shp=zvb_studyarea, 
#                              predictors=preds_m,
#                              type="Spatial_detection", 
#                              mesh_edge=2.5, prange=20, psigma=10, fit_spde = TRUE)
# save(mz_spatial2, file="./output/model_outputs/global_models/global_m3_spatial_vb.R")
# 
