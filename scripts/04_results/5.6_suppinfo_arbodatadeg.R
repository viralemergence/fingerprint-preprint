



# setup objects
setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")


# read model objects
load("./output/model_outputs/wnv_testmodels/jac_spde_mod1.R"); m1_jac = m1
load("./output/model_outputs/wnv_testmodels/jac_spde_mod2.R"); m2_jac = m2

load("./output/model_outputs/wnv_testmodels/wnv_spde_mod1.R"); m1_wnv = m1
load("./output/model_outputs/wnv_testmodels/wnv_spde_mod2.R"); m2_wnv = m2

load("./output/model_outputs/wnv_testmodels/lacrosse_spde_mod1.R"); m1_lac = m1
load("./output/model_outputs/wnv_testmodels/lacrosse_spde_mod2.R"); m2_lac = m2

load("./output/model_outputs/wnv_testmodels/powassan_spde_mod1.R"); m1_pow = m1
load("./output/model_outputs/wnv_testmodels/powassan_spde_mod2.R"); m2_pow = m2

# filter fitted effects
results = do.call(
  rbind.data.frame,
  list(
    extractFixedINLA(m1_jac, model_name = "Case incidence") %>% dplyr::mutate(dz = "Jamestown Canyon encephalitis"),
    extractFixedINLA(m2_jac, model_name = "Outbreak event risk") %>% dplyr::mutate(dz = "Jamestown Canyon encephalitis"),
    
    extractFixedINLA(m1_wnv, model_name = "Case incidence") %>% dplyr::mutate(dz = "West Nile fever"),
    extractFixedINLA(m2_wnv, model_name = "Outbreak event risk") %>% dplyr::mutate(dz = "West Nile fever"),
    
    extractFixedINLA(m1_lac, model_name = "Case incidence") %>% dplyr::mutate(dz = "LaCrosse encephalitis"),
    extractFixedINLA(m2_lac, model_name = "Outbreak event risk") %>% dplyr::mutate(dz = "LaCrosse encephalitis"),
    
    extractFixedINLA(m1_pow, model_name = "Case incidence") %>% dplyr::mutate(dz = "Powassan encephalitis"),
    extractFixedINLA(m2_pow, model_name = "Outbreak event risk") %>% dplyr::mutate(dz = "Powassan encephalitis")
  )) %>%
  dplyr::filter(param != "Intercept")

param_lookup = data.frame(param = unique(results$param))
param_lookup$p2 = NA
param_lookup$p2[ param_lookup$param == "urban_cover"] = "Urban cover"
param_lookup$p2[ param_lookup$param == "forest_cover"] = "Forest cover"
param_lookup$p2[ param_lookup$param == "forest_loss"] = "Forest loss"
param_lookup$p2[ param_lookup$param == "evi_dissimilarity"] = "Vegetation heterogeneity"
param_lookup$p2[ param_lookup$param == "precip_change"] = "Precipitation change"
param_lookup$p2[ param_lookup$param == "tmean_change"] = "Temperature change"
param_lookup$p2[ param_lookup$param == "biodiv_intact"] = "Biodiversity Intactness index"
param_lookup$p2[ param_lookup$param == "social_vulnerability"] = "Social vulnerability"
param_lookup$p2[ param_lookup$param == "livestock_log"] = "Livestock density (log)"
param_lookup$p2 = factor(param_lookup$p2, 
                         levels = c("Urban cover", "Social vulnerability", "Livestock density (log)",
                                    "Forest cover", "Forest loss", "Vegetation heterogeneity", "Biodiversity Intactness index",
                                    "Precipitation change", "Temperature change"))

p1 = results %>%
  dplyr::mutate(dz = factor(dz,
                            levels = c("West Nile fever", "LaCrosse encephalitis", "Jamestown Canyon encephalitis", "Powassan encephalitis"),
                            ordered=TRUE)) %>%
  dplyr::left_join(param_lookup) %>%
  ggplot() + 
  geom_point(aes(p2, mean, group=model, color=model), position=position_dodge(width=0.4)) +
  geom_linerange(aes(p2, ymin=lower, ymax=upper, group=model, color=model), position=position_dodge(width=0.4)) +
  theme_classic() +
  geom_hline(yintercept=0,lty=2) +
  facet_wrap(~dz, scales="free") +
  coord_flip() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size=13)) +
  ylab("Slope parameter estimate")+
  xlab("Driver") +
  scale_color_viridis_d(begin=0.3, end=0.7, name="Model type") +
  theme(legend.position ="bottom",
        legend.text = element_text(size=11), legend.title = element_text(size=12))

ggsave(p1, file="./output/plots/ExtendedData_ModelSpecificationExercise.jpg", 
       device="jpg", units = "in", dpi=600, width=9, height=5.5)
