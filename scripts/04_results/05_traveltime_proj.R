

# ============ Examining effects of travel time to healthcare ===============

setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA)
source("./scripts/00_plot_themes.R")

params = list.files("./output/model_outputs/disease_params_htt/", pattern=".csv", full.names=TRUE)
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
                param = replace(param, param == "health_travel", "Healthcare travel time"),
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
  dplyr::filter(param == "Healthcare travel time")

# htt %>%
#   ggplot() + 
#   geom_point(aes(Disease, mean)) + 
#   geom_linerange(aes(Disease, ymin=lower, ymax=upper)) +
#   theme_classic() + 
#   geom_hline(lty=2, yintercept=0) + 
#   coord_flip()



# ============ for each disease read in data =============

loc = list.files("./output/model_df/", pattern=".csv", full.names=TRUE)
loc = loc[ -grep("usa|brazil", loc)]
dd = data.frame()
ht_ests = data.frame()

for(i in loc){
  
  cc = read.csv(i)

  # median travel time across background locations
  median_tt_bg = median(cc$health_travel[ cc$presence == 0 ], na.rm=TRUE)
  median_tt_pres = median(cc$health_travel[ cc$presence == 1 ], na.rm=TRUE)
  
  # get ranges of distance
  meandist = mean(cc$health_travel, na.rm=TRUE)
  mindist = min(cc$health_travel, na.rm=TRUE)
  maxdist = quantile(cc$health_travel, 0.975, na.rm=TRUE)
  
  # % locations > 1 or 2 hours away
  bg = cc$health_travel[ cc$presence == 0 ]
  percent_over_1hr = sum(bg > 60, na.rm=TRUE) / length(bg)
  percent_over_2hr = sum(bg > 120, na.rm=TRUE) / length(bg)
  
  d_i = data.frame(
    dz = i,
    sc = sc,
    mindist = mindist,
    meandist = meandist,
    maxdist = maxdist,
    median_tt_bg = median_tt_bg,
    median_tt_pres = median_tt_pres,
    percent_over_1hr = percent_over_1hr,
    percent_over_2hr = percent_over_2hr
  )
  dd = rbind(dd, d_i)
  
  # raw
  ht_ests = rbind(
    ht_ests, 
    cc %>% dplyr::select(presence, health_travel) %>% dplyr::mutate(dz = i)
  )
  
}

dd$Disease = c("Anthrax", "CCHF", "Chagas", "Chikungunya",
               "Dengue", "Ebola", "EEE", "H5N1", "Hantavirus",
               "Jamestown Canyon", "Japanese enceph.",
               "Junin", "LaCrosse", "Lassa", "Lyme", "Marburg", "Mayaro",
               "Melioidosis", "MERS", "Mpox", "Nipah", 
               "Oropouche", "P. knowlesi", "Plague", "Powassan", 
               "Brazil spotted", "Rift Valley", "SLE", 
               "West Nile", "Yellow fever", "Zika")

ht_ests = left_join(ht_ests, dd %>% dplyr::select(dz, Disease) %>% distinct()) %>%
  dplyr::left_join(
    read.csv("scripts/04_results/dz_abbrevs.csv") %>% dplyr::select(Disease, abbrev6)
  )

# add into df
htt = left_join(htt, dd)

# median and range stats
median((exp(htt$median)-1)*100)
quantile((exp(htt$median)-1)*100)

# read in proportion of population > X hours under motorized and non-motorized 
res = do.call(
  rbind.data.frame, 
  lapply(list.files("./output/model_outputs/tthc_per_disease/", pattern=".csv", full.names=TRUE), read.csv)
) %>%
  dplyr::filter(mins_tt == 60) %>%
  dplyr::left_join(
    read.csv("scripts/04_results/dz_abbrevs.csv") %>% dplyr::select(Disease, abbrev6)
  ) %>%
  dplyr::mutate(prop_pop = 1 - prop_pop)



# ==================== plots ========================

# remove SLE; uncertainty too high for viz
htt = htt %>% dplyr::filter(Disease != "SLE")

# order of dz points
#fac_order = htt %>% dplyr::arrange(desc(mean)) %>% dplyr::select(abbrev6) 
fac_order = htt %>% dplyr::arrange(mean) %>% dplyr::select(abbrev6) 

# parameter estimates
htt_params = htt %>%
  dplyr::mutate(abbrev6 = factor(abbrev6, levels=fac_order$abbrev6, ordered=TRUE)) %>%
  ggplot() + 
  geom_linerange(aes(abbrev6, ymin=(exp(lower)-1)*100, ymax=(exp(upper)-1)*100), color=viridis::viridis(200)[100]) +
  geom_point(aes(abbrev6, (exp(median)-1)*100), color=viridis::viridis(200)[80], size=2) + 
  theme_classic() + 
  geom_hline(lty=2, yintercept=0) + 
  theme(panel.grid.major = element_line(color="grey95")) +
  coord_flip() + 
  theme(axis.title = element_text(size=11), 
        axis.text = element_text(size=10)) +
  ylab("% change in odds of outbreak report per\nadditional hour's travel time to health facility") +
  xlab("Disease")

# descriptive
outbreak_vs_background = htt %>%
  dplyr::filter(Disease != "SLE") %>% # uncertainty too high for viz
  dplyr::mutate(abbrev6 = factor(abbrev6, levels=fac_order$abbrev6, ordered=TRUE)) %>% 
  dplyr::select(abbrev6, median_tt_bg, median_tt_pres) %>%
  dplyr::rename("Background"=2, "Outbreak"=3) %>%
  reshape2::melt(id.vars=1) %>%
  dplyr::mutate(variable = ifelse(variable == "Background", "Background\nlocations", "Outbreak\nlocations")) %>%
  ggplot() + 
  geom_bar(aes(abbrev6, value/60, group=variable, fill=variable), width=0.45, color=NA, stat="identity", position=position_dodge(width=0.58)) + 
  theme_classic() + 
  theme(panel.grid.major = element_line(color="grey95")) +
  coord_flip() + 
  ylab("Median motorized travel time\nto health facility (hours)") +
  xlab("Disease") +
  theme(#legend.position = c(0.82, 0.16),
        legend.position = c(0.85, 0.8),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7.5),
        axis.title = element_text(size=11), 
        axis.text = element_text(size=10),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_fill_viridis_d(begin=0.3, end=0.7, direction=-1, guide = guide_legend(reverse = TRUE, override.aes = list(size=2)))

# % background locations over 2 hours away
tt_ests = htt %>%
  dplyr::mutate(abbrev6 = factor(abbrev6, levels=fac_order$abbrev6, ordered=TRUE)) %>%
  ggplot() +
  geom_bar(aes(abbrev6, percent_over_2hr*100), width=0.35, color=NA, fill=pals::kovesi.rainbow(200)[150], stat="identity", position=position_dodge(width=0.75)) +
  #geom_bar(aes(abbrev6, prop_pop*100, group=type, fill=type), width=0.65, color=NA, stat="identity", position=position_dodge(width=0.75)) +
  theme_classic() +
  theme(panel.grid.major = element_line(color="grey95")) +
  coord_flip() +
  ylab("% locations\n>2 hours away") +
  xlab("Disease") +
  theme(#legend.position = c(0.82, 0.16),
        legend.position = c(0.84, 0.77),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=7.5),
        axis.title = element_text(size=11),
        axis.text = element_text(size=10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_viridis_d(begin=0.3, end=0.7, direction=1, guide = guide_legend(reverse = TRUE), option="magma")

# next two are % of population
# tt_ests1 has both motorized and nonmotorized
# tt_ests = res %>%
#   dplyr::filter(abbrev6 %in% htt$abbrev6) %>%
#   dplyr::mutate(abbrev6 = factor(abbrev6, levels=fac_order$abbrev6, ordered=TRUE)) %>% 
#   ggplot() + 
#   #geom_bar(aes(abbrev6, prop_pop*100), width=0.65, color=NA, fill="coral1", stat="identity", position=position_dodge(width=0.75)) + 
#   geom_bar(aes(abbrev6, prop_pop*100, group=type, fill=type), width=0.65, color=NA, stat="identity", position=position_dodge(width=0.75)) + 
#   theme_classic() + 
#   theme(panel.grid.major = element_line(color="grey95")) +
#   coord_flip() + 
#   ylab("% population\n>2 hours away") +
#   xlab("Disease") +
#   theme(legend.position = c(0.82, 0.16), 
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size=7.5),
#         axis.title = element_text(size=11), 
#         axis.text = element_text(size=10),
#         axis.text.y = element_blank(), 
#         axis.title.y = element_blank()) + 
#   scale_fill_viridis_d(begin=0.3, end=0.7, direction=1, guide = guide_legend(reverse = TRUE), option="magma")
# 
# # just motorized
# tt_ests = res %>%
#   dplyr::filter(type == "motorized") %>%
#   dplyr::filter(abbrev6 %in% htt$abbrev6) %>%
#   dplyr::mutate(abbrev6 = factor(abbrev6, levels=fac_order$abbrev6, ordered=TRUE)) %>% 
#   ggplot() + 
#   geom_bar(aes(abbrev6, prop_pop*100), width=0.35, color=NA, fill=pals::parula(200)[140], stat="identity", position=position_dodge(width=0.75)) + 
#   #geom_bar(aes(abbrev6, prop_pop*100, group=type, fill=type), width=0.65, color=NA, stat="identity", position=position_dodge(width=0.75)) + 
#   theme_classic() + 
#   theme(panel.grid.major = element_line(color="grey95")) +
#   coord_flip() + 
#   ylab("% population\n>1 hour away") +
#   xlab("Disease") +
#   theme(legend.position = c(0.82, 0.16), 
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(size=7.5),
#         axis.title = element_text(size=11), 
#         axis.text = element_text(size=10),
#         axis.text.y = element_blank(), 
#         axis.title.y = element_blank()) + 
#   scale_fill_viridis_d(begin=0.3, end=0.7, direction=1, guide = guide_legend(reverse = TRUE), option="magma")

pc = gridExtra::grid.arrange(htt_params, outbreak_vs_background, tt_ests, ncol=3, widths=c(1, 0.8, 0.3))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           fontface = "plain", size = 22, 
                           x = c(0.01, 0.45, 0.835), y = c(1, 1, 1))

ggsave(pc, file="./output/plots/Figure5_TravelTimeProj.jpg", device="jpg", units="in", width=11.5, height=4.4, dpi=900, scale=0.85)




