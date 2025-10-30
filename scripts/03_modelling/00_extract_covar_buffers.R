

# ================ EXTRACT COVARIATES IN DIFFERENT SIZED BUFFERS ===================

# ------------ housekeeping ----------------

# setup objects
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)
library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")



# ---------------- 1. create global study area polygon ---------------

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



# --------- 2. extract covariates at background points across globe -----------

num_points = 1000
study_area = study_area
method = "popweight"
pop_res = "10km"
seed = 1000

# generate points
print(paste("Generating", num_points, "background points using method:", method, sep=" "))

# set random seed for reproducibility
if(!is.na(seed)){ set.seed(seed) }
if(!method %in% c("popweight", "random")) return("Error: method must be 'popweight' or 'random")
if(!pop_res %in% c("1km", "10km")) return("Error: pop_res must be 1km or 10km")

if(method == "popweight"){
  
  # create log population raster
  if(pop_res == "1km"){
    pd = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_1km_Aggregated.tif") %>% crop(study_area)
  }
  if(pop_res == "10km"){
    pd = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_10km_agg.tif") %>% crop(study_area)
  }
  pd = raster::mask(pd, fasterize::fasterize(study_area, raster(pd), field = "id"))
  pd = log(pd + 1)
  vv = values(pd)
  vv[ vv == 0 ] = median(vv[ vv != 0 & !is.na(vv)])
  values(pd) = vv
  
  # sample random points (num_points) across log pop raster (could have enmSdm dependency removed)
  pas = enmSdm::sampleRast(pd, n=num_points, replace=TRUE, prob=TRUE)
  if(is.null(pas)){ pas = enmSdm::sampleRast(pd, n=npts, replace=TRUE, prob=TRUE) } # replicate to deal with buggy
  pas = as.data.frame(pas) %>%
    dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
    dplyr::mutate(record_id = paste("bg", 1:num_points, sep="_"))
}

# generate buffers of different sizes 
print("Generating background point buffers")

bg = st_as_sf(pas, coords = c("Longitude", "Latitude"), crs = 4326)
bg2500 = st_buffer(bg, 2500) %>%
  dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
  dplyr::mutate(ADMcode = paste("bg", 1:num_points, sep="_"))

bg5000 = st_buffer(bg, 5000) %>%
  dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
  dplyr::mutate(ADMcode = paste("bg", 1:num_points, sep="_"))

bg10000 = st_buffer(bg, 10000) %>%
  dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
  dplyr::mutate(ADMcode = paste("bg", 1:num_points, sep="_"))

bg20000 = st_buffer(bg, 20000) %>%
  dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
  dplyr::mutate(ADMcode = paste("bg", 1:num_points, sep="_"))

# extract covariates at different scales
print("Extracting covariate data")
covs_df = read.csv("./scripts/03_modelling/covariate_lookup.csv")
covs_df = covs_df %>% dplyr::filter(!cov %in% c("precip_anomaly", "tmean_anomaly",
                                                "evi_coefvar", "forest_loss", 
                                                "hunting", "precip", "tmean"))

# combine into one big sf for extraction
ddx = 
  bg2500 %>% dplyr::mutate(buffer = "2500m") %>%
  rbind(bg5000 %>% dplyr::mutate(buffer = "5000m")) %>%
  rbind(bg10000 %>% dplyr::mutate(buffer = "10000m")) %>%
  rbind(bg20000 %>% dplyr::mutate(buffer = "20000m"))
  
disease_name = "buffer_test"
for(cc in covs_df$cov){
  print(cc)
  extr_x = extractCovariateVals(ddx, covariate_name=cc) # function stored in '00_covar_extraction_funcs.R'
  write.csv(extr_x, paste("./output/model_df/buffer_comparison/bg_", cc, ".csv", sep=""), row.names=FALSE)
  ddx = cbind(ddx, extr_x)
}

write.csv(
  ddx %>% sf::st_drop_geometry(),
  paste("./output/model_df/buffer_comparison/all_covars_multiplebuffers.csv", sep=""),
  row.names = FALSE
)

# view correlations between covars at different buffer sizes

# pivot longer then pivot wider for easy comparison
ll = ddx %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-record_id, -Longitude, -Latitude) %>%
  reshape2::melt(id.vars = 1:2)

cor_m = data.frame()
for(i in unique(ll$variable)){
  a = ll %>% 
    dplyr::filter(variable == i) %>%
    dplyr::mutate(buffer = paste("X", buffer, sep="")) %>%
    tidyr::pivot_wider(
      names_from=buffer, 
      values_from=value
    )
  cor_m = rbind(cor_m, a)
}

# correlations
cors = cor_m %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    cor_5000_2500 = as.vector(cor.test(X5000m, X2500m)$estimate),
    cor_5000_10000 = as.vector(cor.test(X5000m, X10000m)$estimate),
    cor_5000_20000 = as.vector(cor.test(X5000m, X20000m)$estimate)
  )

p1 = cors %>%
  dplyr::rename("covar"=1) %>%
  reshape2::melt(id.vars=1) %>%
  dplyr::mutate(
    variable = as.character(variable),
    variable = replace(variable, variable == "cor_5000_2500", "2.5km"),
    variable = replace(variable, variable == "cor_5000_10000", "10km"),
    variable = replace(variable, variable == "cor_5000_20000", "20km"),
    variable = factor(variable, levels=c("2.5km", "10km", "20km"), ordered=TRUE),
    covar = factor(covar, levels=rev(unique(covar)), ordered=TRUE) 
  ) %>%
  ggplot() + 
  geom_tile(aes(variable, covar, fill=value)) +
  theme_bw() + 
  scale_fill_gradientn(colors=viridisLite::mako(200), name="Correlation\nto 5km\nbuffer") +
  ylab("Covariate") + 
  xlab("Buffer size")

ggsave(p1, file="./output/plots/ED_buffer_correlations.jpg", device="jpg",
       units="in", width=4, height=3, dpi=600)

write.csv(cor_m, "./output/buffer_correlations.csv", row.names=FALSE)
