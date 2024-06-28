

setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")

params = list.files("./output/model_outputs/disease_params/", pattern=".csv", full.names=TRUE)
params = params[ -grep("infocriteria", params)]
pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
pp$mean = round(pp$mean, 3)
pp$lower = round(pp$lower, 3)
pp$upper = round(pp$upper, 3)
pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0

# add abbrevs
pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
  dplyr::mutate(Disease = abbrev2)

pp = pp %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
  dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
                param = replace(param, param == "livestock_log", "Livestock density (log)"),
                param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
                param = replace(param, param == "hunting", "Hunting Pressure Index"),
                param = replace(param, param == "forest_cover", "Forest cover"),
                param = replace(param, param == "forest_loss", "Forest loss"),
                param = replace(param, param == "crop_expansion", "Cropland expansion"),
                param = replace(param, param == "evi_dissimilarity", "Vegetation heterogeneity"),
                param = replace(param, param == "urban_expansion", "Urban expansion"),
                param = replace(param, param == "urban_cover", "Urban cover"),
                param = replace(param, param == "crop_cover", "Cropland cover"),
                param = replace(param, param == "mining", "Mining cover"),
                param = replace(param, param == "tmean_change", "Temperature change"),
                param = replace(param, param == "precip_change", "Precipitation change"),
                param = replace(param, param == "social_vulnerability", "Social vulnerability"),
                param = replace(param, param == "biodiv_intact", "Biodiversity Intactness index"),
                param = replace(param, param == "protected_areas", "Protected area cover"),
                param = replace(param, param == "popdens_log", "Population density (log)")) 

# # add NAs for non-modelled vars
# pp = pp %>% dplyr::select(param, Disease, type, mean, sig, Num_observations)
# add_nas = expand.grid(unique(pp$param), unique(pp$Disease), unique(pp$type)) %>%
#   as.data.frame() %>%
#   dplyr::rename("param" = 1, "Disease"=2, "type"=3) %>%
#   dplyr::mutate(sig = NA, mean = NA) %>% 
#   dplyr::left_join(pp %>% dplyr::select(Disease, Num_observations) %>% distinct())
# add_nas = add_nas[ ! paste(add_nas$param, add_nas$Disease, add_nas$type, sep="_") %in% paste(pp$param, pp$Disease, pp$type, sep="_"), ]
# pp = rbind(pp, add_nas)


# # order by n points
# np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
# pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)


# ============== subset ===============

htt = pp %>% dplyr::filter() %>%
  dplyr::filter(type == "Causal (strict)") %>%
  dplyr::filter(param == "Healthcare travel time (log)")

htt %>%
  ggplot() + 
  geom_point(aes(Disease, mean)) + 
  geom_linerange(aes(Disease, ymin=lower, ymax=upper)) +
  theme_classic() + 
  geom_hline(lty=2, yintercept=0) + 
  coord_flip()



# ============ for each disease read in data =============

loc = list.files("./output/model_df/", pattern=".csv", full.names=TRUE)
dd = data.frame()

for(i in loc){
  
  cc = read.csv(i)

  # health travel data
  cc$health_travel_log = log(cc$health_travel+1)
  
  # get scaling factor from data
  sc = sd(cc$health_travel_log, na.rm=TRUE)
  
  # calculate mean travel time across background points
  mean_tt_bg = mean(cc$health_travel[ cc$presence == 0 ], na.rm=TRUE)
  mean_tt_pres = mean(cc$health_travel[ cc$presence == 1 ], na.rm=TRUE)
  
  # get ranges of distance
  meandist = mean(cc$health_travel, na.rm=TRUE)
  mindist = min(cc$health_travel, na.rm=TRUE)
  maxdist = quantile(cc$health_travel, 0.975, na.rm=TRUE)
  
  d_i = data.frame(
    dz = i,
    sc = sc,
    mindist = mindist,
    meandist = meandist,
    maxdist = maxdist,
    mean_tt_bg = mean_tt_bg,
    mean_tt_pres = mean_tt_pres
  )
  dd = rbind(dd, d_i)
  
}

dd = dd[ -grep("brazil|usa", dd$dz), ]
dd$Disease = c("Anthrax", "CCHF", "Chagas", "Chikungunya",
               "Dengue", "Ebola", "EEE", "H5N1", "Hantavirus",
               "Jamestown Canyon", "Japanese enceph.",
               "Junin", "LaCrosse", "Lassa", "Lyme", "Marburg", "Mayaro",
               "Melioidosis", "MERS", "Mpox", "Nipah", 
               "Oropouche", "P. knowlesi", "Plague", "Powassan", 
               "Brazil spotted", "Rift Valley", "SLE", 
               "West Nile", "Yellow fever", "Zika")

# add scaling factor
htt = left_join(htt, dd)

# rescale params
htt[ , c("mean", "median", "lower", "upper")] = htt[ , c("mean","median", "lower", "upper")] / htt$sc


# # ============== get data for all diseases =================
# 
# # load three submodels
# load(file="./output/model_outputs/global_models/global_m1_nonspatial.R")
# load(file="./output/model_outputs/global_models/global_m2_spatial.R")
# load(file="./output/model_outputs/global_models/global_m3_spatial.R")
# 
# # effects
# all_models = do.call(
#   rbind.data.frame,
#   list(m_nonspatial1$fx, m_spatial1$fx, m_spatial2$fx)
# ) %>%
#   dplyr::mutate(
#     model_type = ifelse(grepl("Nonspatial", type), "Non-spatial", "Spatial"),
#     detection = ifelse(grepl("nondetection", type), "Non-detection", "Detection"),
#     type2 = paste(model_type, detection, sep=" + "),
#     type3 = "No adjustment",
#     type3 = replace(type3, type2 == "Spatial + Non-detection", "Geospatial"),
#     type3 = replace(type3, type2 == "Spatial + Detection", "Geospatial + Detection"),
#     type2 = factor(type2, levels = rev(c("Non-spatial + Non-detection", "Non-spatial + Detection", "Spatial + Detection")), ordered=TRUE),
#     type3 = factor(type3, levels = c("No adjustment", "Geospatial", "Geospatial + Detection"), ordered=TRUE)
#   ) %>%
#   dplyr::filter(param == "Travel time to health facility (log)")
# 
# scg = read.csv("./output/htt_global_scaling.csv")
# 
# # rescale
# all_models = all_models %>%
#   dplyr::select(mean, lower, median, upper) %>%
#   dplyr::mutate(Disease = "Global (all diseases)")
# all_models[, 1:4] = all_models[, 1:4]  * scg$sc
# 
# 
# # plot 1 global
# p1 = all_models %>%
#   ggplot() + 
#   geom_point(aes(Disease, (exp(median)-1)*100), color=viridis::viridis(200)[80]) + 
#   geom_linerange(aes(Disease, ymin=(exp(lower)-1)*100, ymax=(exp(upper)-1)*100), color=viridis::viridis(200)[100]) +
#   theme_classic() + 
#   geom_hline(lty=2, yintercept=0) + 
#   theme(panel.grid.major = element_line(color="grey95")) +
#   ylim(-100, 90) + 
#   coord_flip() + 
#   ylab("% change in odds of outbreak event report\nper logarithmic unit of distance to clinic") +
#   xlab("Disease")
  
# descriptive
p1 = htt %>%
  dplyr::arrange(desc(mean)) %>%
  # dplyr::mutate(
  #   lower = replace(lower, Disease == "Marburg", NA),
  #   upper = replace(upper, Disease == "Marburg", NA)
  # )%>%
  dplyr::filter(Disease != "Marburg") %>% # uncertainty too high for viz
  dplyr::mutate(abbrev6 = factor(abbrev6, levels=abbrev6, ordered=TRUE)) %>% 
  dplyr::select(abbrev6, mean_tt_bg, mean_tt_pres) %>%
  dplyr::rename("Background"=2, "Outbreak"=3) %>%
  reshape2::melt(id.vars=1) %>%
  dplyr::mutate(variable = ifelse(variable == "Background", "Study area population", "Outbreak locations")) %>%
  ggplot() + 
  geom_bar(aes(abbrev6, value/60, group=variable, fill=variable), width=0.55, color=NA, stat="identity", position=position_dodge(width=0.65)) + 
  theme_classic() + 
  theme(panel.grid.major = element_line(color="grey95")) +
  coord_flip() + 
  ylab("Mean estimated travel time\nto health centre (hours)") +
  xlab("Disease") +
  theme(legend.position = c(0.78, 0.16), 
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=8.5),
        axis.title = element_text(size=11), 
        axis.text = element_text(size=10),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_fill_viridis_d(begin=0.3, end=0.7, direction=1, guide = guide_legend(reverse = TRUE))

# plot 
p2 = htt %>%
  dplyr::arrange(desc(mean)) %>%
  dplyr::filter(Disease != "Marburg") %>% # uncertainty too high for viz
  # dplyr::mutate(
  #   lower = replace(lower, Disease == "Marburg", NA),
  #   upper = replace(upper, Disease == "Marburg", NA)
  # ) %>%
  dplyr::mutate(abbrev6 = factor(abbrev6, levels=abbrev6, ordered=TRUE)) %>%
  ggplot() + 
  geom_linerange(aes(abbrev6, ymin=(exp(lower)-1)*100, ymax=(exp(upper)-1)*100), color=viridis::viridis(200)[100]) +
  geom_point(aes(abbrev6, (exp(median)-1)*100), color=viridis::viridis(200)[80]) + 
  theme_classic() + 
  geom_hline(lty=2, yintercept=0) + 
  theme(panel.grid.major = element_line(color="grey95")) +
  coord_flip() + 
  theme(axis.title = element_text(size=11), 
        axis.text = element_text(size=10)) +
  ylab("% change in odds of outbreak report\nper logarithmic unit of distance to health centre") +
  xlab("Disease")


library(patchwork)

pc = p2 + p1 + plot_layout(ncol=2, widths=c(1, 1))

ggsave(pc, file="./output/plots/ExtendedData_TravelTimeProj.jpg", device="jpg", units="in", width=10, height=4.2, dpi=600, scale=0.85)




# for each rescale by scaling factor to get on the log scale
htt_rs = data.frame()

for(i in 1:nrow(htt)){
  
  h_i = htt[i, ]
  
  dd_i = dd %>% dplyr::filter(Disease == h_i$Disease)
  
  h_i[ , c("mean", "median", "lower", "upper")] = h_i[ , c("mean","median", "lower", "upper")] * dd_i$sc
  h_i$scale = "rescaled"
  h_i$sc = 
  
  # dd_i
  # 
  # newdat = data.frame(
  #   health_travel = seq(dd_i$meandist, dd_i$maxdist, length.out = 250)
  # ) %>%
  #   dplyr::mutate(
  #     health_travel_log = log(health_travel+1),
  #     htl_scaled = health_travel_log / dd_i$sc # scale
  #   ) 
  # 
  # # project change in odds
  # newdat$log_odds = h_i$mean * newdat$htl_scaled
  # newdat$OR = exp(newdat$log_odds)
  # plot(newdat$health_travel, newdat$OR)
  
  # combine
  htt_rs = rbind(htt_rs, h_i)
  
}

# express slope as % of 50% decline on log scale
dec = log(0.5) * sc
h_i$mean / dec


dec_resc = dec * sc
exp(dec_resc) 

log(0.5)

100 * exp(-0.69)
5 + log(0.5)


# ============== specify criterion =============

library(patchwork)



# ========== create combined figure across all diseases =============


# --------- 1. Hypotheses 1: majority rule ----------

hyp = "Causal (broad)"

rr = pp %>%
  dplyr::filter(type == hyp) %>%
  dplyr::group_by(param) %>%
  dplyr::rename("signif"=sig) %>%
  dplyr::summarise(
    sig = n_distinct(Disease[ signif == TRUE]),
    pos = n_distinct(Disease[ signif == TRUE & mean > 0 ]),
    neg = n_distinct(Disease[ signif == TRUE & mean < 0 ]),
    tested = n_distinct(Disease),
    prop_tested = tested / n_distinct(pp$Disease),
    null = tested - sig
  ) %>%
  dplyr::arrange(desc(sig)) %>%
  dplyr::mutate(
    rank = rank(-sig, ties.method = "min")
  )

rank_drivers =  rr$param

# plot 0 = rank
p0 = rr %>%
  dplyr::mutate(grp = "a",
                param = factor(param, levels=rev(rank_drivers), ordered=TRUE)) %>%
  ggplot() + 
  geom_text(aes(grp, param, label=rank), fontface="bold") + 
  theme_classic() + 
  xlab("Rank") + 
  ylab("Driver") +
  theme(axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_line(color="grey75"),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(color="white", size=10.5),
        axis.line = element_blank(),
        panel.border = element_rect(color="grey75", fill=NA)) + 
  theme(panel.background = element_rect(fill="white"))

# plot 1 = directionality - including "not tested"
p1 = rr %>%
  dplyr::mutate(not_tested = n_distinct(pp$Disease) - tested) %>%
  tidyr::pivot_longer(cols=c("pos", "neg", "null", "not_tested")) %>%
  dplyr::mutate(param = factor(param, levels=rev(rank_drivers), ordered=TRUE),
                name = replace(name, name == "pos", "Pos"),
                name = replace(name, name == "neg", "Neg"),
                name = replace(name, name == "null", "None"),
                name = replace(name, name == "not_tested", "NT"),
                value = replace(value, value==0, NA),
                name = factor(name, levels=c("NT", "Neg", "None", "Pos"), ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(name, param, fill=name, size=value, group=name), pch=21, alpha=1) + 
  theme_classic() + 
  xlab("Relationship") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size=11.5),
        legend.text = element_text(size=11),
        panel.border = element_rect(color="grey75", fill=NA),
        axis.ticks.x = element_line(color="grey75"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(color="black", size=10.5)) +
  scale_fill_manual(values=c("Pos"="violetred3", "Neg"="steelblue2", "None"="white", "NT"="grey85"), name="Direction", guide="none") + 
  scale_size(name = "Num. diseases", breaks=c(1, 5, 10, 15, 20, 25)) +
  theme(panel.background = element_rect(fill="white"))

# plot 1 = directionality - without "Not tested"
# p1 = rr %>%
#   tidyr::pivot_longer(cols=c("pos", "neg", "null")) %>%
#   dplyr::mutate(param = factor(param, levels=rev(rank_drivers), ordered=TRUE),
#                 name = replace(name, name == "pos", "Pos"),
#                 name = replace(name, name == "neg", "Neg"),
#                 name = replace(name, name == "null", "None"),
#                 value = replace(value, value==0, NA)) %>%
#   ggplot() + 
#   geom_point(aes(name, param, fill=name, size=value, group=name), pch=21, alpha=1) + 
#   theme_classic() + 
#   xlab("Relationship") + 
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.line = element_blank(),
#         legend.title = element_text(size=11.5),
#         legend.text = element_text(size=11),
#         panel.border = element_rect(color="grey75", fill=NA),
#         axis.ticks.x = element_line(color="grey75"),
#         axis.ticks.y = element_blank(),
#         axis.title.x = element_text(size=11),
#         axis.text.x = element_text(color="black", size=10.5)) +
#   scale_fill_manual(values=c("Pos"="violetred3", "Neg"="steelblue2", "None"="grey80"), name="Direction", guide="none") + 
#   scale_size(name = "Num. diseases") +
#   theme(panel.background = element_rect(fill="white"))


# plot 2 = effect sizes
cc = pp %>%
  dplyr::filter(type == hyp)
cc$vartype = "Land use change/intensity"
cc$vartype[ cc$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
cc$vartype[ cc$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
cc$vartype[ cc$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
cc$vartype[ cc$param %in% c("Forest cover", "Cropland cover", "Vegetation heterogeneity", "Biodiversity Intactness index") ] = "Ecosystem structure"
cc$param = factor(cc$param, levels=rev(rank_drivers), ordered=TRUE)
cc$vartype = factor(cc$vartype, levels=c("Detection", "Socioeconomic", "Ecosystem structure", "Climate change", "Land use change/intensity"))

p2 = ggplot() + 
  geom_linerange(data=cc, aes(y=param, xmin=lower, xmax=upper, group=Disease, color=vartype), alpha=0.25, size=0.8) +
  geom_point(data=cc[ cc$sig == FALSE, ], aes(y=param, x=mean, group=Disease, color=vartype), fill=NA, pch=21, size=4, alpha=0.7) +
  geom_point(data=cc[ cc$sig == TRUE, ], aes(y=param, x=mean, group=Disease, fill=vartype), pch=21, size=4, alpha=0.7) + 
  geom_vline(xintercept=0, lty=2) + 
  theme_classic() + 
  xlim(-max(abs(range(c(cc$lower, cc$upper)))), max(abs(range(c(cc$lower, cc$upper))))) +
  xlab("Scaled covariate effect on outbreak risk (log odds)") + 
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type") +
  scale_color_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size=11.5),
        legend.text = element_text(size=11),
        panel.border = element_rect(color="grey75", fill=NA),
        axis.ticks.x = element_line(color="grey75"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(color="black", size=10.5)) +
  theme(panel.background = element_rect(fill="white"))
  
# create combined plot with shared legend
combined = p0 + p1 + p2 & theme(legend.position = "right")

# including/excluding "not tested"
combined = combined + plot_layout(guides = "collect", widths=c(0.2, 0.65, 2.5)) # incl
#combined = combined + plot_layout(guides = "collect", widths=c(0.2, 0.55, 2.5)) # excl

# save
ggsave(combined, file="./output/plots/Figure3_majorityrule.jpg", device="jpg", units="in", width=12.7, height=5, dpi=600)





# --------- 1. Hypotheses 2: Top ranked ----------

hyp = "Causal (strict)"

rr = pp %>%
  dplyr::filter(type == hyp) %>%
  dplyr::group_by(param) %>%
  dplyr::rename("signif"=sig) %>%
  dplyr::summarise(
    sig = n_distinct(Disease[ signif == TRUE]),
    pos = n_distinct(Disease[ signif == TRUE & mean > 0 ]),
    neg = n_distinct(Disease[ signif == TRUE & mean < 0 ]),
    tested = n_distinct(Disease),
    prop_tested = tested / n_distinct(pp$Disease),
    null = tested - sig
  ) %>%
  dplyr::arrange(desc(sig)) %>%
  dplyr::mutate(
    rank = rank(-sig, ties.method = "min")
  )

rank_drivers =  rr$param

# plot 0 = rank
p0 = rr %>%
  dplyr::mutate(grp = "a",
                param = factor(param, levels=rev(rank_drivers), ordered=TRUE)) %>%
  ggplot() + 
  geom_text(aes(grp, param, label=rank), fontface="bold") + 
  theme_classic() + 
  xlab("Rank") + 
  ylab("Driver") +
  theme(axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_line(color="grey75"),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(color="white", size=10.5),
        axis.line = element_blank(),
        panel.border = element_rect(color="grey75", fill=NA)) + 
  theme(panel.background = element_rect(fill="white"))

# plot 1 = directionality - including "not tested"
p1 = rr %>%
  dplyr::mutate(not_tested = n_distinct(pp$Disease) - tested) %>%
  tidyr::pivot_longer(cols=c("pos", "neg", "null", "not_tested")) %>%
  dplyr::mutate(param = factor(param, levels=rev(rank_drivers), ordered=TRUE),
                name = replace(name, name == "pos", "Pos"),
                name = replace(name, name == "neg", "Neg"),
                name = replace(name, name == "null", "None"),
                name = replace(name, name == "not_tested", "NT"),
                value = replace(value, value==0, NA),
                name = factor(name, levels=c("NT", "Neg", "None", "Pos"), ordered=TRUE)) %>%
  ggplot() + 
  geom_point(aes(name, param, fill=name, size=value, group=name), pch=21, alpha=1) + 
  theme_classic() + 
  xlab("Relationship") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size=11.5),
        legend.text = element_text(size=11),
        panel.border = element_rect(color="grey75", fill=NA),
        axis.ticks.x = element_line(color="grey75"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(color="black", size=10.5)) +
  scale_fill_manual(values=c("Pos"="violetred3", "Neg"="steelblue2", "None"="white", "NT"="grey85"), name="Direction", guide="none") + 
  scale_size(name = "Num. diseases", breaks=c(1, 5, 10, 15, 20, 25)) +
  theme(panel.background = element_rect(fill="white"))

# plot 1 = directionality - without "Not tested"
# p1 = rr %>%
#   tidyr::pivot_longer(cols=c("pos", "neg", "null")) %>%
#   dplyr::mutate(param = factor(param, levels=rev(rank_drivers), ordered=TRUE),
#                 name = replace(name, name == "pos", "Pos"),
#                 name = replace(name, name == "neg", "Neg"),
#                 name = replace(name, name == "null", "None"),
#                 value = replace(value, value==0, NA)) %>%
#   ggplot() + 
#   geom_point(aes(name, param, fill=name, size=value, group=name), pch=21, alpha=1) + 
#   theme_classic() + 
#   xlab("Relationship") + 
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.line = element_blank(),
#         legend.title = element_text(size=11.5),
#         legend.text = element_text(size=11),
#         panel.border = element_rect(color="grey75", fill=NA),
#         axis.ticks.x = element_line(color="grey75"),
#         axis.ticks.y = element_blank(),
#         axis.title.x = element_text(size=11),
#         axis.text.x = element_text(color="black", size=10.5)) +
#   scale_fill_manual(values=c("Pos"="violetred3", "Neg"="steelblue2", "None"="grey80"), name="Direction", guide="none") + 
#   scale_size(name = "Num. diseases") +
#   theme(panel.background = element_rect(fill="white"))


# plot 2 = effect sizes
cc = pp %>%
  dplyr::filter(type == hyp)
cc$vartype = "Land use change/intensity"
cc$vartype[ cc$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
cc$vartype[ cc$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
cc$vartype[ cc$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
cc$vartype[ cc$param %in% c("Forest cover", "Cropland cover", "Vegetation heterogeneity", "Biodiversity Intactness index") ] = "Ecosystem structure"
cc$param = factor(cc$param, levels=rev(rank_drivers), ordered=TRUE)
cc$vartype = factor(cc$vartype, levels=c("Detection", "Socioeconomic", "Ecosystem structure", "Climate change", "Land use change/intensity"))

p2 = ggplot() + 
  geom_linerange(data=cc, aes(y=param, xmin=lower, xmax=upper, group=Disease, color=vartype), alpha=0.25, size=0.8) +
  geom_point(data=cc[ cc$sig == FALSE, ], aes(y=param, x=mean, group=Disease, color=vartype), fill=NA, pch=21, size=4, alpha=0.7) +
  geom_point(data=cc[ cc$sig == TRUE, ], aes(y=param, x=mean, group=Disease, fill=vartype), pch=21, size=4, alpha=0.7) + 
  geom_vline(xintercept=0, lty=2) + 
  theme_classic() + 
  xlim(-max(abs(range(c(cc$lower, cc$upper)))), max(abs(range(c(cc$lower, cc$upper))))) +
  xlab("Scaled covariate effect on outbreak risk (log odds)") + 
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type") +
  scale_color_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        legend.title = element_text(size=11.5),
        legend.text = element_text(size=11),
        panel.border = element_rect(color="grey75", fill=NA),
        axis.ticks.x = element_line(color="grey75"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=11),
        axis.text.x = element_text(color="black", size=10.5)) +
  theme(panel.background = element_rect(fill="white"))

# create combined plot with shared legend
combined = p0 + p1 + p2 & theme(legend.position = "right")

# including/excluding "not tested"
combined = combined + plot_layout(guides = "collect", widths=c(0.2, 0.65, 2.5)) # incl
#combined = combined + plot_layout(guides = "collect", widths=c(0.2, 0.55, 2.5)) # excl

# save
ggsave(combined, file="./output/plots/Figure3_topranked.jpg", device="jpg", units="in", width=12.7, height=5, dpi=600)

