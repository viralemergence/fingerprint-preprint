
# ====================== SUPP FIGURE: MULTIPANEL EFFECTS SUMMARY ACROSS ALL DISEASES =========================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")

params = list.files("./output/model_outputs/disease_params/", pattern=".csv", full.names=TRUE)
params = params[ -grep("infocriteria", params)]
pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0

# add abbrevs
pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
  dplyr::mutate(Disease = abbrev6)

pp = pp %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
  dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
                param = replace(param, param == "livestock_log", "Livestock density (log)"),
                param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
                param = replace(param, param == "hunting", "Hunting Pressure Index"),
                param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
                param = replace(param, param == "forest_cover", "Forest cover"),
                param = replace(param, param == "forest_loss", "Forest loss"),
                param = replace(param, param == "crop_expansion", "Cropland expansion"),
                param = replace(param, param == "evi_dissimilarity", "Vegetation heterogeneity"),
                param = replace(param, param == "evi_coefvar", "EVI variance"),
                param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
                param = replace(param, param == "urban_expansion", "Urban expansion"),
                param = replace(param, param == "urban_cover", "Built-up land"),
                param = replace(param, param == "crop_cover", "Cropland cover"),
                param = replace(param, param == "mining", "Mining"),
                param = replace(param, param == "tmean_change", "Temperature change"),
                param = replace(param, param == "precip_change", "Precipitation change"),
                param = replace(param, param == "social_vulnerability", "Social vulnerability"),
                param = replace(param, param == "biodiv_intact", "Biodiversity Intactness Index"),
                param = replace(param, param == "protected_areas", "Protected area coverage"),
                param = replace(param, param == "popdens_log", "Population density (log)")) 

# order by n points
np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)

# ranked driver order from causal strict models
rr = pp %>%
  dplyr::filter(type == "Causal (strict)") %>%
  dplyr::group_by(param) %>%
  dplyr::rename("signif"=sig) %>%
  dplyr::summarise(
    sig = n_distinct(Disease[ signif == TRUE]),
    pos = n_distinct(Disease[ signif == TRUE & mean > 0 ]),
    neg = n_distinct(Disease[ signif == TRUE & mean < 0 ]),
    tested = n_distinct(Disease),
    prop_tested = sig / tested,
    null = tested - sig
  ) %>%
  dplyr::arrange(desc(sig)) %>%
  dplyr::mutate(
    rank = rank(-sig, ties.method = "min")
  )

rank_drivers =  rr$param
param_order = c(rank_drivers)

pp = pp %>%
  dplyr::mutate(param = factor(param, levels=rev(param_order), ordered=TRUE))

pp$type2 = "Univariate"
pp$type2[ pp$type == "Causal (broad)" ] = "Hypotheses 1 (majority rule)"
pp$type2[ pp$type == "Causal (strict)" ] = "Hypotheses 2 (top-ranked)"
pp$type2 = factor(pp$type2, levels = c("Univariate", "Hypotheses 1 (majority rule)", "Hypotheses 2 (top-ranked)"))

plot_width = 17.25

# multipanel plot
csx = MetBrewer::met.brewer("Benedictus")
figure2 = pp %>%
  dplyr::filter(!is.na(mean)) %>%
  ggplot() +
  #geom_tile(aes(Disease, param, fill=mean)) + 
  geom_point(aes(Disease, param, fill=mean), size=5, pch=21, stroke=0.25) + 
  geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(Disease, param), col="black", size=5, pch=21, fill=NA, stroke=1) + 
  #scale_fill_gradient2(low=scales::muted(csx[1]), high=scales::muted(csx[length(csx)]), na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
  scale_fill_gradient2(low=csx[length(csx)], high=csx[1], na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
  theme_classic() + 
  coord_fixed() + 
  theme(axis.text.x = element_text(size=12, angle=90, vjust=0.4), 
        #axis.text.x = element_blank(), 
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14),
        #axis.title = element_blank(),
        strip.text = element_text(size=14.5),
        panel.grid.major.x = element_line(size=0.2, color="grey78"),
        strip.background = element_blank()) + 
  scale_x_discrete(position="bottom") + 
  scale_y_discrete() +
  theme(legend.position="bottom") + 
  xlab("Disease") + ylab("Driver") +
  facet_wrap(~type2, nrow=3)
ggsave(figure2, file="./output/plots/SuppFigure_MultiDiseaseMatrix.png", device="png", units="in", dpi=600, width=8.5, height=13.5, scale=0.95)













# ====================== OLD COLOUR SCHEMES ===================
# 
# #csx = colorRampPalette(MetBrewer::met.brewer("Cassatt1"))(200)
# csx = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200))
# lims = c(-max(abs(pp$mean), na.rm=TRUE), max(abs(pp$mean), na.rm=TRUE))
# figure2 = pp %>%
#   dplyr::filter(!is.na(mean)) %>%
#   ggplot() +
#   #geom_tile(aes(Disease, param, fill=mean)) + 
#   geom_point(aes(Disease, param, fill=mean), size=5, pch=21, stroke=0.25) + 
#   geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(Disease, param), col="black", size=5, pch=21, fill=NA, stroke=1) + 
#   scale_fill_gradientn(colours = csx, limits=lims, na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   theme_classic() + 
#   coord_fixed() + 
#   theme(axis.text.x = element_text(size=12, angle=90, vjust=0.4), 
#         #axis.text.x = element_blank(), 
#         axis.text.y = element_text(size=12), 
#         axis.title = element_text(size=14),
#         #axis.title = element_blank(),
#         strip.text = element_text(size=14.5),
#         strip.background = element_blank()) + 
#   scale_x_discrete(position="bottom") + 
#   scale_y_discrete() +
#   theme(legend.position="top") + 
#   xlab("Disease") + ylab("Driver") +
#   facet_wrap(~type2, nrow=1)
# ggsave(figure2, file="./output/plots/Figure2_matrix_blobtest_cs3.png", device="png", units="in", dpi=600, width=plot_width, height=6.5, scale=0.95)
# 
# 
# # try setting hard boundaries to colour scale
# csx = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(200))
# lims = c(-3, 3)
# figure2 = pp %>%
#   dplyr::filter(!is.na(mean)) %>%
#   dplyr::mutate(mean = replace(mean, mean < -3, -3),
#                 mean = replace(mean, mean > 3, 3)) %>%
#   ggplot() +
#   #geom_tile(aes(Disease, param, fill=mean)) + 
#   geom_point(aes(Disease, param, fill=mean), size=5, pch=21, stroke=0.25) + 
#   geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(Disease, param), col="black", size=5, pch=21, fill=NA, stroke=1) + 
#   scale_fill_gradientn(colours = csx, limits=lims, na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   theme_classic() + 
#   coord_fixed() + 
#   theme(axis.text.x = element_text(size=12, angle=90, vjust=0.4), 
#         #axis.text.x = element_blank(), 
#         axis.text.y = element_text(size=12), 
#         axis.title = element_text(size=14),
#         #axis.title = element_blank(),
#         strip.text = element_text(size=14.5),
#         strip.background = element_blank(),
#         panel.grid.major.x = element_line(size=0.25, color="grey90"),
#         panel.grid.minor.x = element_line(size=0.25, color="grey90")) + 
#   scale_x_discrete(position="bottom") + 
#   scale_y_discrete() +
#   theme(legend.position="top") + 
#   xlab("Disease") + ylab("Driver") +
#   facet_wrap(~type2, nrow=1)
# ggsave(figure2, file="./output/plots/Figure2_matrix_blobtest_cs3_lim3.png", device="png", units="in", dpi=600, width=plot_width, height=6.5, scale=0.95)
# 
# 
# # other color scales
# csx = MetBrewer::met.brewer("Benedictus")
# figure2 = pp %>%
#   dplyr::filter(!is.na(mean)) %>%
#   dplyr::mutate(mean = replace(mean, mean < -3, -3),
#                 mean = replace(mean, mean > 3, 3)) %>%
#   ggplot() +
#   #geom_tile(aes(Disease, param, fill=mean)) + 
#   geom_point(aes(Disease, param, fill=mean), size=5, pch=21, stroke=0.25) + 
#   geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(Disease, param), col="black", size=5, pch=21, fill=NA, stroke=1) + 
#   #scale_fill_gradient2(low=scales::muted(csx[1]), high=scales::muted(csx[length(csx)]), na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   scale_fill_gradient2(low=csx[length(csx)], high=csx[1], na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   theme_classic() + 
#   coord_fixed() + 
#   theme(axis.text.x = element_text(size=12, angle=90, vjust=0.4), 
#         #axis.text.x = element_blank(), 
#         axis.text.y = element_text(size=12), 
#         axis.title = element_text(size=14),
#         #axis.title = element_blank(),
#         strip.text = element_text(size=14.5),
#         strip.background = element_blank(),
#         panel.grid.major.x = element_line(size=0.25, color="grey90"),
#         panel.grid.minor.x = element_line(size=0.25, color="grey90")) + 
#   scale_x_discrete(position="bottom") + 
#   scale_y_discrete() +
#   theme(legend.position="top") + 
#   xlab("Disease") + ylab("Driver") +
#   facet_wrap(~type2, nrow=1)
# ggsave(figure2, file="./output/plots/Figure2_matrix_blobtest_cs2_lim3.png", device="png", units="in", dpi=600, width=plot_width, height=6.5, scale=0.95)
# 
# 
# # 
# # 
# # # flip to vertical
# # csx = MetBrewer::met.brewer("Benedictus")
# # figure2 = pp %>%
# #   dplyr::filter(!is.na(mean)) %>%
# #   dplyr::mutate(mean = replace(mean, mean < -3, -3),
# #                 mean = replace(mean, mean > 3, 3)) %>%
# #   ggplot() +
# #   #geom_tile(aes(Disease, param, fill=mean)) + 
# #   geom_point(aes(Disease, param, fill=mean), size=5, pch=21, stroke=0.25) + 
# #   geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(Disease, param), col="black", size=5, pch=21, fill=NA, stroke=1) + 
# #   #scale_fill_gradient2(low=scales::muted(csx[1]), high=scales::muted(csx[length(csx)]), na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
# #   scale_fill_gradient2(low=csx[length(csx)], high=csx[1], na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
# #   theme_classic() + 
# #   coord_fixed() + 
# #   theme(axis.text.x = element_text(size=12, angle=90, vjust=0.4), 
# #         #axis.text.x = element_blank(), 
# #         axis.text.y = element_text(size=12), 
# #         axis.title = element_text(size=14),
# #         #axis.title = element_blank(),
# #         strip.text = element_text(size=14.5),
# #         strip.background = element_blank(),
# #         panel.grid.major.x = element_line(size=0.25, color="grey90"),
# #         panel.grid.minor.x = element_line(size=0.25, color="grey90")) + 
# #   scale_x_discrete(position="bottom") + 
# #   scale_y_discrete() +
# #   theme(legend.position="top") + 
# #   xlab("Disease") + ylab("Driver") +
# #   facet_wrap(~type2, ncol=1)
# # ggsave(figure2, file="./output/plots/Figure2_matrix_blobtest_cs2_lim3_vertical.png", device="png", units="in", dpi=600, width=7.5, height=16, scale=0.95)
