

# ================ Runs individual disease models using INLA geospatial regression ===================

# Runs region-specific models for dengue and yellow fever for supplement
# Individual disease scripts stored in "./scripts/03_modelling/inla_arboregions/"

# setup objects
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr) 
library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)#; library(enmSdm)

# generate covariate dataset param
# TRUE = re-extract covariates for each disease
# FALSE = just re-run models (if covars already extracted?)
generate_dataset = FALSE

# scripts per disease
scr = list.files("./scripts/03_modelling/inla_arboregions/", full.names=TRUE, pattern=".R")

# run each script in turn, print details, clear workspace after each
for(scr_i in scr){
  
  print(scr_i)
  source(scr_i)
  
  
  # --------- sample posterior for fixed effects ------------
  
  print("Extracting samples")
  
  for(m in 1:length(models)){
    
    mm = models[[m]]
  
    # samples
    s2 = inla.posterior.sample(n=2000, mm$model)
    samp_extr = function(x){
      s = s2[[x]]$latent
      s = s[ which(dimnames(s)[[1]] %in% paste(predictors, ":1", sep="")), ]
      t(as.data.frame(s))
    }
    s2_samp = do.call(rbind.data.frame, lapply(1:length(s2), samp_extr)) %>%
      dplyr::mutate(region = mm$region, id=sample(1:10^6, 1)) %>%
      dplyr::mutate(Disease = dz$Disease[1])
    
    write.csv(s2_samp, file=paste("./output/model_outputs/disease_pooled_arboregions/", disease_name, "_", mm$region, "_majorityrule_samp.csv", sep=""), row.names=FALSE)
    
    # fixed effects summary
    f2 = extractFixedINLA(mm$model) %>%
      dplyr::mutate(region = mm$region, Disease = dz$Disease[1])
    write.csv(f2, paste("./output/model_outputs/disease_pooled_arboregions/", disease_name, "_", mm$region, "_params_summary.csv", sep=""), row.names=FALSE)
    
  }
  
  # # ---- model objects for plotting and visualisation -----
  # 
  # # predicted spatial field from "majority rule" model
  # spde_pred = rasteriseSPDE(m2$summary.random$s$mean, mesh1, study_area)
  # names(spde_pred) = disease_name
  # crs(spde_pred) = st_crs(study_area)
  # 
  # model_obj = list(data = dd %>% dplyr::select(Longitude, Latitude, presence) %>% dplyr::mutate(Disease = dz$Disease[1]),
  #                  spde = spde_pred, 
  #                  study_area = study_area %>% dplyr::mutate(Disease = dz$Disease[1]))
  # save(model_obj, file=paste("./output/model_outputs/disease_pooled/", disease_name, "_objects.R", sep=""))
  # 
  # 
  # -----------------------
  
  # clear workspace (keeping generate_dataset and scr)
  wspace = ls()
  wspace = wspace[ -which(wspace %in% c("scr", "generate_dataset"))]
  rm(list=wspace)
  
}

# ends