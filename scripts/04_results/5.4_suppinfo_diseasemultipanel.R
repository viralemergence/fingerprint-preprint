

# ====================== SUPP FIGURE: MULTIPANEL EFFECTS SUMMARY ACROSS ALL DISEASES =========================

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf);# library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")

# multivariate
loc = "./output/model_outputs/disease_pooled/"
ll = list.files(loc, pattern = "summary", full.names=TRUE)
pp = do.call(rbind.data.frame, lapply(ll, read.csv)) 

# univariate
params = list.files("./output/model_outputs/disease_params/", pattern="params.csv", full.names=TRUE)
pp2 = do.call(rbind.data.frame, lapply(params, read.csv)) 
pp2 = pp2 %>% dplyr::filter(type == "Univariate") %>% dplyr::mutate(type = "Univariable")

# combine
pp = rbind(pp, pp2)
pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0

# add abbrevs
pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
  dplyr::mutate(Disease = abbrev6)

pp = pp %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
  # dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
  #               param = replace(param, param == "livestock_log", "Livestock density (log)"),
  #               param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
  #               param = replace(param, param == "hunting", "Hunting Pressure Index"),
  #               param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
  #               param = replace(param, param == "forest_cover", "Forest cover"),
  #               param = replace(param, param == "forest_loss", "Forest loss"),
  #               param = replace(param, param == "crop_expansion", "Cropland expansion"),
  #               param = replace(param, param == "evi_dissimilarity", "Landscape heterogeneity"),
  #               param = replace(param, param == "evi_coefvar", "EVI variance"),
  #               param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
  #               param = replace(param, param == "urban_expansion", "Urban expansion"),
  #               param = replace(param, param == "urban_cover", "Urban cover"),
  #               param = replace(param, param == "crop_cover", "Cropland cover"),
  #               param = replace(param, param == "mining", "Mining"),
  #               param = replace(param, param == "tmean_change", "Temperature change"),
  #               param = replace(param, param == "precip_change", "Precipitation change"),
  #               param = replace(param, param == "social_vulnerability", "Social vulnerability"),
  #               param = replace(param, param == "biodiv_intact", "Biodiversity Intactness Index"),
  #               param = replace(param, param == "protected_areas", "Protected area coverage"),
  #               param = replace(param, param == "popdens_log", "Population density (log)")) %>%
  dplyr::mutate(param = replace(param, param == "health_travel_log", "Health travel"),
                param = replace(param, param == "livestock_log", "Livestock"),
                param = replace(param, param == "livestock_ruminants_log", "Livestock"),
                param = replace(param, param == "hunting", "HPI (Hunting)"),
                param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
                param = replace(param, param == "forest_cover", "Forest"),
                param = replace(param, param == "forest_loss", "Forest loss"),
                param = replace(param, param == "crop_expansion", "Crop expansion"),
                param = replace(param, param == "evi_dissimilarity", "Land heterogeneity"),
                param = replace(param, param == "evi_coefvar", "EVI variance"),
                param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
                param = replace(param, param == "urban_expansion", "Urb expansion"),
                param = replace(param, param == "urban_cover", "Urban"),
                param = replace(param, param == "crop_cover", "Cropland"),
                param = replace(param, param == "mining", "Mining"),
                param = replace(param, param == "tmean_change", "Temp. change"),
                param = replace(param, param == "precip_change", "Precip. change"),
                param = replace(param, param == "social_vulnerability", "Vulnerability"),
                param = replace(param, param == "biodiv_intact", "BII (Biodiversity)"),
                param = replace(param, param == "protected_areas", "Protected area"),
                param = replace(param, param == "popdens_log", "Population density (log)")) 


# order by n points
np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)

# ranked driver order
rr = pp %>%
  dplyr::filter(type == "Majority rule") %>%
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

pp$type2 = pp$type
pp$type2[ pp$type == "Any" ] = "Any author"
pp$type2[ pp$type == "Top-ranked" ] = "Top ranked"
pp$type2 = factor(pp$type2, levels = c("Univariable", "Any author", "Majority rule", "Top ranked"))

plot_width = 17.25

# rescale variables 
pp = pp %>%
  dplyr::mutate(
    mean2 = mean,
    mean2 = replace(mean2, mean2 >= 2, "> 2"),
    mean2 = replace(mean2, mean2 >= 1 & mean2 < 2, "1 to 2"),
    mean2 = replace(mean2, mean2 >= 0.5 & mean2 < 1, "0.5 to 1"),
    mean2 = replace(mean2, mean2 >= 0 & mean2 < 0.5, "0 to 0.5"),
    mean2 = replace(mean2, mean2 < 0 & mean2 < -0.5, "0 to -0.5"),
    mean2 = replace(mean2, mean2 <= -0.5 & mean2 > -1, "-0.5 to -1"),
    mean2 = replace(mean2, mean2 <= -1 & mean2 > -2, "-1 to -2"),
    mean2 = replace(mean2, mean2 <= -2, "< -2")
  )

# bins = data.frame(bin = c(-4, -2, -1, -0.5, 0, 0.5, 1, 2, 4)) %>% 
#   dplyr::mutate(
#     n = 1:length(bin),
#     binname = c("< -2", "-2 to -1", "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1", "1 to 2", "> 2", "erp")
#   )

# bin variables for categorical colour scheme
bins = data.frame(bin = c(-4, -1.5, -0.5, 0, 0.5, 1.5, 4)) %>% 
  dplyr::mutate(
    n = 1:length(bin),
    binname = c("< -1.5", "-1.5 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1.5", "> 1.5", "erp")
  )
                                                                         
pp$n = .bincode(pp$mean, breaks = bins$bin)
pp = pp %>% dplyr::left_join(bins)
#pp$binname = factor(pp$binname, levels=c("< -2", "-2 to -1", "-1 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1", "1 to 2", "> 2"), ordered=TRUE)
pp$binname = factor(pp$binname, levels=c("< -1.5", "-1.5 to -0.5", "-0.5 to 0", "0 to 0.5", "0.5 to 1.5", "> 1.5"), ordered=TRUE)

pp %>% ggplot() + geom_boxplot( aes(x = factor(binname), y = mean, group=binname) )
table(pp$binname)

# multipanel plot
#csx = MetBrewer::met.brewer("Benedictus")
# figure2 = pp %>%
#   dplyr::filter(!is.na(mean)) %>%
#   ggplot() +
#   #geom_tile(aes(Disease, param, fill=mean)) + 
#   geom_point(aes(Disease, param, fill=mean), size=5, pch=21, stroke=0.2) + 
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
#         panel.grid.major.x = element_line(size=0.2, color="grey78"),
#         strip.background = element_blank()) + 
#   scale_x_discrete(position="bottom") + 
#   scale_y_discrete() +
#   theme(legend.position="bottom") + 
#   xlab("Disease") + ylab("Driver") +
#   facet_wrap(~type2, nrow=3) 
# ggsave(figure2, file="./output/plots/SuppFigure_MultiDiseaseMatrix.png", device="png", units="in", dpi=600, width=8.5, height=13.5, scale=0.95)

pp = pp %>% dplyr::mutate(param = factor(param, levels=rev(param_order), ordered=TRUE))
pp$Disease = factor(pp$Disease, levels = rev(np$Disease), ordered=TRUE)

# figure2 = pp %>%
#   dplyr::filter(!is.na(mean)) %>%
#   ggplot() +
#   geom_point(aes(param, Disease, fill=mean), size=5, pch=22, stroke=0.2, color="grey20") + 
#   geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(param, Disease), col="grey20", size=5, pch=22, fill=NA, stroke=1) + 
#   #scale_fill_gradient2(low=scales::muted(csx[1]), high=scales::muted(csx[length(csx)]), na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   scale_fill_gradient2(low=csx[length(csx)], high=csx[1], na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   theme_classic() + 
#   coord_fixed() + 
#   theme(axis.text.x = element_text(size=13, angle=90, vjust=0.4), 
#         #axis.text.x = element_blank(), 
#         axis.text.y = element_text(size=13), 
#         axis.title = element_text(size=14),
#         #axis.title = element_blank(),
#         axis.line = element_blank(),
#         strip.text = element_text(size=16),
#         panel.background = element_rect(fill="grey85"),
#         panel.border = element_rect(color="grey20", fill=NA),
#         panel.grid.major.y = element_line(size=0.2, color="grey73"),
#         strip.background = element_blank()) + 
#   scale_x_discrete(position="bottom") + 
#   scale_y_discrete() +
#   theme(legend.position="bottom") + 
#   xlab("Driver") + ylab("Disease") +
#   facet_wrap(~type2, nrow=1) 
# ggsave(figure2, file="./output/plots/SuppFigure_MultiDiseaseMatrix_lg.png", device="png", units="in", dpi=600, width=20, height=10, scale=0.95)

# calculate % of possible dz-driver combs per type
poss = n_distinct(pp$param) * n_distinct(pp$Disease)
pc = as.data.frame(table(pp$type2) / poss) %>%
  dplyr::mutate(type3 = paste(Var1, " (", round(Freq, 2)*100, "%)", sep="")) %>%
  dplyr::select(1, 3) %>%
  dplyr::rename("type2"=1)
pp = left_join(pp, pc)
pp$type3 = factor(pp$type3, levels = c("Univariable (92%)", "Any author (77%)", "Majority rule (69%)", "Top ranked (50%)"))

# with a discrete colour scale

csx = as.vector(MetBrewer::met.brewer(name="Benedictus", n=n_distinct(pp$binname)))
#csx = RColorBrewer::brewer.pal(n=n_distinct(pp$binname), name = "BrBG")
csx = pals::coolwarm(n=n_distinct(pp$binname))

figure2 = pp %>%
  dplyr::filter(!is.na(mean)) %>%
  ggplot() +
  geom_point(aes(Disease, param, fill=binname), size=4.7, pch=22, stroke=0.1, color="grey55") + 
  geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(Disease, param), col="grey20", size=4.7, pch=22, fill=NA, stroke=1.2) + 
  #scale_fill_gradient2(low=scales::muted(csx[1]), high=scales::muted(csx[length(csx)]), na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
  scale_fill_manual(values=csx, na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
  theme_classic() + 
  coord_fixed() + 
  theme(axis.text.x = element_text(size=13, angle=90, vjust=0.4), 
        #axis.text.x = element_blank(), 
        axis.text.y = element_text(size=13), 
        axis.title = element_text(size=18),
        #axis.title = element_blank(),
        legend.text = element_text(size=13), 
        legend.title = element_text(size=13.5),
        axis.line = element_blank(),
        strip.text = element_text(size=18),
        panel.background = element_rect(fill="grey92"),
        panel.border = element_rect(color="grey20", fill=NA),
        panel.grid.major.x = element_line(size=0.2, color="grey78"),
        strip.background = element_blank()) + 
  scale_x_discrete(position="bottom") + 
  scale_y_discrete() +
  theme(legend.position="bottom") + 
  xlab("Disease") + ylab("Driver") +
  facet_wrap(~type3, ncol=1) +
  guides(
    fill = guide_legend(nrow = 2)
  ) 
ggsave(figure2, file="./output/plots/SuppFigure_MultiDiseaseMatrix_lg3.png", device="png", units="in", dpi=600, width=10, height=18, scale=0.95)


# figure2 = pp %>%
#   dplyr::filter(!is.na(mean)) %>%
#   ggplot() +
#   geom_point(aes(param, Disease, fill=binname), size=5, pch=21, stroke=0.35, color="grey30") + 
#   geom_point(data = pp[ pp$sig == TRUE & !is.na(pp$sig), ], aes(param, Disease), col="grey20", size=5, pch=21, fill=NA, stroke=1) + 
#   #scale_fill_gradient2(low=scales::muted(csx[1]), high=scales::muted(csx[length(csx)]), na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   scale_fill_manual(values=csx, na.value = "grey70", name="Effect on outbreak risk\n(log odds)") + 
#   theme_classic() + 
#   coord_fixed() + 
#   theme(axis.text.x = element_text(size=13, angle=90, vjust=0.4), 
#         #axis.text.x = element_blank(), 
#         axis.text.y = element_text(size=13), 
#         axis.title = element_text(size=14),
#         #axis.title = element_blank(),
#         legend.text = element_text(size=13), 
#         legend.title = element_text(size=13.5),
#         axis.line = element_blank(),
#         strip.text = element_text(size=16),
#         panel.background = element_rect(fill="grey90"),
#         panel.border = element_rect(color="grey20", fill=NA),
#         panel.grid.major.y = element_line(size=0.2, color="grey65"),
#         strip.background = element_blank()) + 
#   scale_x_discrete(position="bottom") + 
#   scale_y_discrete() +
#   theme(legend.position="bottom") + 
#   xlab("Driver") + ylab("Disease") +
#   facet_wrap(~type3, nrow=1) +
#   guides(
#     fill = guide_legend(nrow = 1)
#   ) 
# ggsave(figure2, file="./output/plots/SuppFigure_MultiDiseaseMatrix_lg2.png", device="png", units="in", dpi=600, width=20, height=10, scale=0.95)
# 






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
