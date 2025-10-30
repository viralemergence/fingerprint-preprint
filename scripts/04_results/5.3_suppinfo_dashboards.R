
# ====================== SUPP FIGURE: PARAMETER DASHBOARDS FOR ALL DISEASES =========================

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA);
source("./scripts/00_plot_themes.R")

# result location
loc = "./output/model_outputs/disease_pooled/"
ll = list.files(loc, pattern = "summary", full.names=TRUE)

# read in results
params = do.call(rbind.data.frame, lapply(ll, read.csv)) %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::mutate(param = replace(param, param == "health_travel_log", "Health travel"),
                param = replace(param, param == "livestock_log", "Livestock"),
                param = replace(param, param == "forest_cover", "Forest"),
                param = replace(param, param == "forest_loss", "Forest loss"),
                param = replace(param, param == "crop_expansion", "Crop expansion"),
                param = replace(param, param == "evi_dissimilarity", "Veg heterogeneity"),
                param = replace(param, param == "urban_expansion", "Urban expansion"),
                param = replace(param, param == "urban_cover", "Urban"),
                param = replace(param, param == "hunting", "HPI (Hunting)"),
                param = replace(param, param == "crop_cover", "Cropland"),
                param = replace(param, param == "mining", "Mining"),
                param = replace(param, param == "tmean_change", "Temp change"),
                param = replace(param, param == "precip_change", "Precip change"),
                param = replace(param, param == "social_vulnerability", "Vulnerability"),
                param = replace(param, param == "biodiv_intact", "BII (Biodiversity)"),
                param = replace(param, param == "protected_areas", "Protected area"),
                param = replace(param, param == "popdens_log", "Population density (log)")) %>%
  dplyr::filter(type != "Univariate") %>%
  dplyr::filter(type != "Consensus") %>%
  dplyr::mutate(#type = ifelse(type == "Causal (broad)", "Hypothesis 1 (majority rule)", "Hypothesis 2 (top-ranked)"), 
                #type = factor(type, levels=c("Causal (broad)", "Causal (strict)"), ordered=TRUE),
                param = factor(param, levels=rev(unique(param)[order(unique(param))]), ordered=TRUE),
                type = replace(type, type == "Any", "Any author"))

params = params %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
  dplyr::mutate(Disease = abbrev3) %>%
  dplyr::mutate(Disease = replace(Disease, Disease=="Eastern equine encephalitis", "E. equine encephalitis"))

# param order
params$param = factor(params$param,
                      levels = rev(c("Urban",
                                 "Health travel", 
                                 "Forest", 
                                 "Veg heterogeneity", 
                                 "Cropland", 
                                 "Protected area",
                                 "BII (Biodiversity)", 
                                 "Precip change", 
                                 "Livestock",
                                 "Temp change", 
                                 "Forest loss", 
                                 "Vulnerability", 
                                 "Crop expansion", 
                                 "HPI (Hunting)", 
                                 "Mining",
                                 "Urban expansion")),
                      ordered=TRUE)

px = params %>%
  ggplot() + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=1), size=1.1) + 
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.8), size=0.6, alpha=0.7) +
  theme_minimal() + 
  coord_flip() +
  facet_wrap(~Disease, scales="free_x", ncol=6) + 
  scale_colour_viridis_d(end=.7, direction=-1) +
  ylab("Scaled covariate effect on outbreak risk (log odds)") + xlab("Driver") +
  theme(axis.title = element_text(size=22), 
        axis.text.y = element_text(size=13),
        strip.text = element_text(size=18), 
        legend.position = c(0.75, 0.09),
        legend.text = element_text(size=35), 
        legend.title = element_blank()
  )

ggsave(px,
       file="./output/plots/SuppFigure_Dashboards.jpg", 
       device="jpg",
       units="in", 
       width=17, height=21, 
       dpi=600)
  
  