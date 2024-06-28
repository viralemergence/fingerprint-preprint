

# ================ JUNIN VIRUS (Argentina) ===================


# ------------ housekeeping ----------------

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# random seed and disease name
rseed = 7859
disease_name = "junin"

# define important analysis variables for specific disease
occurrence_point_buffer_km = 10 # km buffer around each occurrence point for spatial uncertainty in spillover event
filter_area_thresh = 5000 # maximum geographical area of occurrence polygons to include in km2 (i.e. exclude v imprecise records)
num_background_points = 1000 # number of background points to generate
covar_coverage_thresh = 0.9 # only keep covariates that are non-NA for at least X% of points (baseline 80%)

# variables with strong hypothesised association with this disease
causal_vars = c("tmean_anomaly", 
                "precip_anomaly",
                "evi_dissimilarity", 
                "urban_cover",
                "crop_cover",
                "mining",
                "health_travel_log",
                "social_vulnerability",
                "biodiv_intact")




# ---------------- 1. Build dataset for modelling (dependencies: 00_analysis_funcs.R) ---------------

# read disease data
dz = read.csv("./output/spillovers_processed/spillovers_junin.csv") %>%
  dplyr::filter(Year >= 1985) %>%
  distinct() %>%
  dplyr::mutate(record_id = 1:length(Disease))

# create study area polygon with US border boundary
bounds = sf::st_read("./data/spillovers/Junin/gadm41_ARG_0.shp")
study_area = createStudyAreaPolygon(lat=dz$Latitude, lon=dz$Longitude, custom_boundary=bounds)

# harmonise points and polygons
# creating X km radius buffer around all point locations, and filtering out any areas >5000km2 (approx 70x70 km square)
dz = harmoniseDataTypes(dz, dz_shp, point_buffer_km=occurrence_point_buffer_km, filter_area_thresh=filter_area_thresh)

# generate population weighted point+buffers with same median polygon area
med_area_size = median(dz$poly_area)
buf_rad = sqrt(med_area_size / pi) * 1000
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
  extr_x = extractCovariateVals(dd, covariate_name=cc) # function stored in '00_covar_extraction_funcs.R'
  ddx = cbind(ddx, extr_x)
}

# log transform any relevant variables, drop geometry and save
ddx = ddx %>%
  dplyr::mutate(health_travel_log = log(health_travel+1),
                livestock_log = log(livestock_all+1)) %>% 
  st_drop_geometry()
filename = paste("./output/model_df/" , disease_name, "_df_wp_10k.csv", sep="")
write.csv(ddx, filename, row.names=FALSE)





# ----------- 2. Initial visualisation: check everything has worked ------------

# p1: study area + polygons
p1 = ggplot() +
  geom_sf(data=study_area, fill="grey95", col="grey80") +
  maptheme +
  geom_sf(data=dz, col="red", fill=NA) + ggtitle(paste(dz$Disease[1], ", n=", nrow(dz), sep=""))

# p1 = study area + background
p2 = ggplot() +
  geom_sf(data=study_area, fill="grey95", col="grey80") +
  maptheme +
  geom_point(data=ddx[ ddx$presence==0,], aes(Longitude, Latitude), size=0.25, col="black") +
  geom_point(data=ddx[ ddx$presence==1,], aes(Longitude, Latitude), size=1, col="red")

# extracted covariates
p3 = ddx %>%
  reshape2::melt(id.vars = 1:5) %>%
  dplyr::filter(!variable %in% c("livestock_all", "health_travel")) %>%
  ggplot() +
  ggforce::geom_sina(aes(x = factor(presence), y=value, group=factor(presence)), width=1, size=0.1, col="grey60") +
  geom_boxplot(aes(x = factor(presence), y=value, group=factor(presence), fill=variable), width=0.25, alpha=0.8) +
  facet_wrap(~variable, scales="free_y") +
  theme_classic() +
  xlab("Presence") +
  ylab("Covariate value") +
  theme(legend.position="none",
        strip.background = element_blank(), strip.text = element_text(size=13))

pc1 = gridExtra::grid.arrange(p1, p2, nrow=2)
pc = gridExtra::grid.arrange(pc1, p3, nrow=1)




# ------------- 3. INLA spatial regression model ---------------

# ------- a. read and prepare data ------

# read data and scale predictors (include log livestock and log health trav)
filename = paste("./output/model_df/" , disease_name, "_df_wp.csv", sep="")
dd = read.csv(filename) %>%
  dplyr::select(-health_travel, -livestock_all, -evi_coefvar)
dd[ , 6:ncol(dd)] = apply(dd[ , 6:ncol(dd)], 2, scale)

# proportion NAs in each predictor: exclude predictors that do not have coverage for at least 80% of data
nas = apply(dd[ dd$presence == 1, 6:ncol(dd)], 2, function(x){ sum(!is.na(x))/length(x) })
nas = names(nas[ nas < covar_coverage_thresh ])

# # check for multicollinearity: lots in this area
# urban correlated with urban expansion and evi coefvar and health travel time (keep urban)
# crop and forest cover both correlated with each other and with biodiversity intactness (keep latter)
usdm::vifstep(dd[ , 6:ncol(dd) ])
corrplot::corrplot(cor(dd[ , 6:ncol(dd) ]))

# specify which to exclude in later analyses if any
mc_vars_to_exclude = c("urban_expansion",  nas)

# examine VIFs/corplots for multivar and causal var subsets
corrplot::corrplot(cor(dd[ , 6:ncol(dd) ] %>% dplyr::select(-mc_vars_to_exclude)))
usdm::vifstep(dd[ , 6:ncol(dd) ] %>% dplyr::select(causal_vars))


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
max.edge = 0.4
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
                           prior.range = c(5, 0.1), 
                           prior.sigma = c(5, 0.01))

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
fm1 = data.frame()
for(i in 1:length(predictors)){
  pred = predictors[i]
  print(paste("Fitting", pred, sep=" "))
  form_i = formula(paste(form_base, pred, sep = " + "))
  m1_i = fitINLAModel(formx = form_i, family=family1, stack=stack1, spde=spde)
  fm1 = rbind(fm1, getIC(m1_i, mod_name = pred))
  m1[[i]] = m1_i
}

# model 2: multivariate model including all predictors except multicollinear
preds_m = predictors[ !predictors %in% mc_vars_to_exclude] 
form = formula(paste(c(form_base, preds_m), collapse=" + "))
m2 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde)

# model 3: hypothesis led (causal) model including only covars with strong hypotheses
preds_c = preds_m[ preds_m %in% causal_vars ]
form = formula(paste(c(form_base, preds_c), collapse=" + "))
m3 = fitINLAModel(formx=form, family=family1, stack=stack1, spde=spde)



# ----- d. save info criteria from all models ----

info_crit = do.call(
  rbind.data.frame,
  list(
    getIC(m0, mod_name = "Baseline"),
    fm1,
    getIC(m2, mod_name = "Multivariate"),
    getIC(m3, mod_name = "Causal")
  )
) %>%
  dplyr::mutate(Disease = dz$Disease[1])
write.csv(info_crit, paste("./output/model_outputs/disease_params/", disease_name,"_infocriteria.csv", sep=""), row.names=FALSE)


# ----- e. extract fixed effects ------

f1 = do.call(
  rbind.data.frame, 
  lapply(m1, function(x) extractFixedINLA(x))
) %>%
  dplyr::mutate(type = "Univariate")

f2 = extractFixedINLA(m2) %>%
  dplyr::mutate(type = "Multivariate")

f3 = extractFixedINLA(m3) %>%
  dplyr::mutate(type = "Causal")

fx = do.call(rbind.data.frame, list(f1, f2, f3)) %>% 
  dplyr::mutate(Disease = dz$Disease[1],
                Num_observations = sum(dd$presence))
row.names(fx) = c()
write.csv(fx, paste("./output/model_outputs/disease_params/", disease_name,"_params.csv", sep=""), row.names=FALSE)


# ----- f. plots ---------

# view spde
rasteriseSPDE(m3$summary.random$s$mean, mesh1, study_area) %>%
  as.data.frame(xy=TRUE) %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=field)) +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black") + 
  geom_point(data=dd[ dd$presence == 1, ], aes(Longitude, Latitude), col="red", size=0.8) + 
  theme(legend.position = c(0.98, 0.15), plot.title = element_text(hjust=0.5, size=16)) +
  ggtitle(dz$Disease[1])

# plot parameters from 3 models
params_plot = fx %>%
  dplyr::filter(param != "Intercept") %>%
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
                param = replace(param, param == "social_vulnerability", "Social vulnerability"),
                param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
                param = replace(param, param == "protected_areas", "Protected area coverage"),
                param = replace(param, param == "popdens_log", "Population density (log)")) %>%
  #dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", ordered=TRUE))) %>%
  dplyr::mutate(type = factor(type, levels=c("Univariate", "Multivariate", "Causal"), ordered=TRUE),
                param = factor(param, levels=rev(unique(param)[order(unique(param))]), ordered=TRUE)) %>%
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
  facet_wrap(~type, ncol=3) + 
  ggtitle(dz$Disease[1])
ggsave(params_plot, file=paste("./output/plots/dz_dashboards/", disease_name, "_dashboard.png", sep=""), device="png", units="in", width=15, height=5.5, dpi=600, scale=0.9)



