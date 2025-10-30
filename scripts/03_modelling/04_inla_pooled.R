

# ================ Runs individual disease models using INLA geospatial regression ===================

# Wrapper script to run multivariable models for all individual diseases (geospatial logistic regression using R-INLA)
# Re-runs fully adjusted models and extracts fixed effects posterior samples - used for revised MS analyses
# These samples are pooled later for Figure 4a, to estimate the mean effect of covariates pooled 
# across all diseases
# Individual disease scripts are stored in "./scripts/03_modelling/inla_pooled/"

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
scr = list.files("./scripts/03_modelling/inla_pooled/", full.names=TRUE)

# run each script in turn, print details, clear workspace after each
for(scr_i in scr){
  
  print(scr_i)
  source(scr_i)
  
  
  # --------- sample posterior for fixed effects ------------
  
  print("Extracting samples")
  
  # model 2: majority rule
  s2 = inla.posterior.sample(n=2000, m2)
  samp_extr = function(x){
    s = s2[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(predictors, ":1", sep="")), ]
    t(as.data.frame(s))
  }
  s2_samp = do.call(rbind.data.frame, lapply(1:length(s2), samp_extr)) %>%
    dplyr::mutate(type = "Majority rule", id=sample(1:10^6, 1)) %>%
    dplyr::mutate(Disease = dz$Disease[1])
  
  # model 3: top ranked
  s3 = inla.posterior.sample(n=2000, m3)
  samp_extr = function(x){
    s = s3[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(predictors, ":1", sep="")), ]
    t(as.data.frame(s))
  }
  s3_samp = do.call(rbind.data.frame, lapply(1:length(s3), samp_extr)) %>%
    dplyr::mutate(type = "Top-ranked", id=sample(1:10^6, 1)) %>%
    dplyr::mutate(Disease = dz$Disease[1])
  
  # model 4: any author
  s4 = inla.posterior.sample(n=2000, m4)
  samp_extr = function(x){
    s = s4[[x]]$latent
    s = s[ which(dimnames(s)[[1]] %in% paste(predictors, ":1", sep="")), ]
    t(as.data.frame(s))
  }
  s4_samp = do.call(rbind.data.frame, lapply(1:length(s4), samp_extr)) %>%
    dplyr::mutate(type = "Any", id=sample(1:10^6, 1)) %>%
    dplyr::mutate(Disease = dz$Disease[1])
  
  # save
  write.csv(s2_samp, file=paste("./output/model_outputs/disease_pooled/", disease_name, "_majorityrule_samp.csv", sep=""), row.names=FALSE)
  write.csv(s3_samp, file=paste("./output/model_outputs/disease_pooled/", disease_name, "_topranked_samp.csv", sep=""), row.names=FALSE)
  write.csv(s4_samp, file=paste("./output/model_outputs/disease_pooled/", disease_name, "_any_samp.csv", sep=""), row.names=FALSE)

  
  # ------------- summary fixed effects estimates -----------
  
  f2 = extractFixedINLA(m2) %>%
    dplyr::mutate(type = "Majority rule")
  f3 = extractFixedINLA(m3) %>%
    dplyr::mutate(type = "Top-ranked")
  f4 = extractFixedINLA(m4) %>%
    dplyr::mutate(type = "Any")
  fx = do.call(rbind.data.frame, list(f2, f3, f4)) %>%
    dplyr::mutate(Disease = dz$Disease[1],
                  Num_observations = sum(dd$presence))
  row.names(fx) = c()
  write.csv(fx, paste("./output/model_outputs/disease_pooled/", disease_name,"_params_summary.csv", sep=""), row.names=FALSE)

  print(
    fx %>%
    dplyr::filter(param != "Intercept") %>%
    dplyr::mutate(type = factor(type, levels=c("Any", "Majority rule", "Top-ranked"), ordered=TRUE)) %>%
    ggplot() + 
    geom_point(aes(x=median, y=param, color=type), position=position_dodge(width=0.75)) + 
    geom_linerange(aes(xmin=lower, xmax=upper, y=param, color=type), position=position_dodge(width=0.75)) + 
    theme_bw() +
    scale_color_viridis_d(end=0.9) + 
    geom_vline(xintercept=0, lty=2) + 
    ggtitle(fx$Disease[1])
  )
  
  
  # ----------- info criteria from all models ------------
  
  info_crit = do.call(
    rbind.data.frame,
    list(
      getIC(m0, mod_name = "Baseline"),
      getIC(m2, mod_name = "Majority rule"),
      getIC(m3, mod_name = "Top-ranked"),
      getIC(m4, mod_name = "Any")
    )
  ) %>%
    dplyr::mutate(Disease = dz$Disease[1])
  write.csv(info_crit, paste("./output/model_outputs/disease_pooled/", disease_name,"_infocriteria.csv", sep=""), row.names=FALSE)
  
  
  # ---- model objects for plotting and visualisation -----
  
  # predicted spatial field from "majority rule" model
  spde_pred = rasteriseSPDE(m2$summary.random$s$mean, mesh1, study_area)
  names(spde_pred) = disease_name
  crs(spde_pred) = st_crs(study_area)
  
  model_obj = list(data = dd %>% dplyr::select(Longitude, Latitude, presence) %>% dplyr::mutate(Disease = dz$Disease[1]),
                   spde = spde_pred, 
                   study_area = study_area %>% dplyr::mutate(Disease = dz$Disease[1]))
  save(model_obj, file=paste("./output/model_outputs/disease_pooled/", disease_name, "_objects.R", sep=""))
  
  
  # -----------------------
  
  # clear workspace (keeping generate_dataset and scr)
  wspace = ls()
  wspace = wspace[ -which(wspace %in% c("scr", "generate_dataset"))]
  rm(list=wspace)
  
}

# ends