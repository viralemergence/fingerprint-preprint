

# ============== Fits the global disease-agnostic model for different regions ================


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



# ---------- identify regions ---------------

# get region
rr = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
ddx = st_as_sf(dd, coords = c("Longitude", "Latitude"), crs=crs(rr))
ii = sf::st_intersects(ddx, rr)
ii = unlist(lapply(ii, function(x) ifelse(length(x)==0, NA, x)))
dd$iso3 = rr$iso3[ ii ]
dd$region = countrycode::countrycode(dd$iso3, origin="iso3c", destination = "region")
dd$region[ dd$iso3 == "ANT" ] = "Latin America & Caribbean"
dd$region[ dd$iso3 == "ESH" ] = "Middle East & North Africa"
dd$region[ dd$iso3 %in% c("MYT", "REU") ] = "Sub-Saharan Africa"



# ------------- fit model for each region ---------------

rgs = c("East Asia & Pacific",
        "Latin America & Caribbean",
        "North America",
        "Sub-Saharan Africa",
        "South Asia")

for(rg in rgs){
  
  # ------- subset data to just region and recreate study area ----------
  print(rg)
  dd_r = dd %>% dplyr::filter(region == rg)
  study_area = createStudyAreaPolygon(lat=dd_r$Latitude, lon=dd_r$Longitude, buffer=FALSE, extent_border = 1)
  
  if(rg == "East Asia & Pacific"){ rgcode = "EAP" }
  if(rg == "Latin America & Caribbean"){  rgcode = "LAC" }
  if(rg == "North America"){  rgcode = "NAM" }
  if(rg == "Sub-Saharan Africa"){  rgcode = "SSA" }
  if(rg == "South Asia"){  rgcode = "SAS" }
  
  
  # ------- check correlation matrix and exclude highly collinear vars  ---------
  
  # corrplot::corrplot(cor(dd_r[ , 5:19 ] %>% dplyr::filter(complete.cases(.))))
  # usdm::vifstep(dd_r[ , 5:19 ])
  if(rg == "East Asia & Pacific"){ preds_rg = preds_m[ !preds_m == "social_vulnerability"] }
  if(rg == "Latin America & Caribbean"){ preds_rg = preds_m }
  if(rg == "North America"){ preds_rg = preds_m[ !preds_m %in% c("social_vulnerability", "health_travel_log")] }
  if(rg == "Sub-Saharan Africa"){ preds_rg = preds_m }
  if(rg == "South Asia"){ preds_rg = preds_m[ !preds_m == "social_vulnerability"] }
  
  
  # ------ model objects and function --------
  
  # likelihood and key priors
  family1 = "binomial"
  hyper1.iid = list(theta = list(prior="pc.prec", param=c(1,0.01)))
  control.fixed1 = list(mean.intercept=0, # prior mean for intercept
                        prec.intercept=0.1, # prior precision for intercept
                        mean=0, # prior mean for fixed effects
                        prec=1)  # prior precision for fixed effects
  
  # build mesh using spatial locs
  spatial_locs = dd_r[ , c("Longitude", "Latitude")]
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
    
    print(
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
    )
    
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
      region = type,
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
  
  
  # -------------- fit model 100 times to bootstrap subsets recovering samples for each -------------
  
  for(i in 1:100){
    
    print(paste("fitting model", i))
    ts = Sys.time()
    
    # subset to 100 points per disease
    npts = dd_r %>% dplyr::group_by(Disease) %>% dplyr::summarise(npts = length(Disease))
    dd_s = dd_r %>% 
      dplyr::filter(Disease != "Background") %>%
      dplyr::filter(Disease %in% npts$Disease[ npts$npts >= 100 ]) %>%
      dplyr::group_by(Disease) %>%
      dplyr::sample_n(100, replace=FALSE) %>%
      rbind(
        dd_r %>% dplyr::filter(Disease %in% npts$Disease[ npts$npts < 100 ])
      )
    dd_s = dd_s %>%
      rbind(
        dd_r %>% dplyr::filter(Disease == "Background") %>%
          dplyr::sample_n(nrow(dd_s)*2, replace=FALSE)
      )
    
    # fit model
    m_s = fitGlobalModel(data=dd_s, 
                         study_shp=study_area, 
                         predictors=preds_rg,
                         type=rg, 
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
      dplyr::mutate(region = rg, id=sample(1:10^6, 1))
    
    # extract spde
    spde = m_s$spde
    
    # save
    write.csv(fx_samp, paste("./output/model_outputs/global_model_bs/outputs/regions/", rgcode, "_", fx_samp$id[1], "_fixefs.csv", sep=""), row.names=FALSE)
    save(spde, file=paste("./output/model_outputs/global_model_bs/outputs/regions/", rgcode, "_", fx_samp$id[1], "_spde.R", sep=""))
    
    te = Sys.time()
    print(te - ts)
    
  } # end of model bootstrapping code
  
  
} # end of multi region loop









