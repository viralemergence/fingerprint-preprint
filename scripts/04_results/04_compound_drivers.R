
# ============== CORRELATION MATRICES AMONG PREDICTORS AND DRIVERS =====================

setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(igraph); library(ggraph); library(tidygraph)
source("./scripts/00_plot_themes.R")





# ================ co-occurrence of significant disease drivers from models ==================

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
                param = replace(param, param == "hunting", "Hunting pressure"),
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
                param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
                param = replace(param, param == "protected_areas", "Protected area coverage"),
                param = replace(param, param == "popdens_log", "Population density (log)")) 

# order by n points
np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)

# fix for rvf
pp$human_infection_route[ pp$Disease == "Rift Valley" ] = "Wildlife + vector"

# # exclude any with no shared drivers (mining)
# pp = pp %>% 
#   dplyr::filter(param != "Mining")

build_graph = function(param_results){
  
  # significant effects only 
  sig_efs = param_results %>% 
    dplyr::filter(sig == TRUE)
  num_dz = n_distinct(sig_efs$Disease)
  
  # 1. nodes including all 
  nodes = data.frame(node = as.vector(unique(pp$param)))
  nodesx = sig_efs %>%
    dplyr::group_by(param) %>%
    dplyr::summarise(n_models = n_distinct(Disease)) %>%
    #dplyr::arrange(desc(n_models)) %>%
    dplyr::rename(node=1) %>%
    dplyr::mutate(node = as.vector(node))
  nodes = dplyr::left_join(nodes, nodesx) %>%
    dplyr::mutate(n_models = replace(n_models, is.na(n_models), 0),
                  prop_models = n_models / num_dz,
                  Occurrence = n_models)
  
  # add abbrevs for visualisation
  nodes = left_join(nodes, read.csv("./scripts/04_results/driver_abbrevs.csv"), by=c("node"="driver")) %>%
    dplyr::mutate(
      abbrev2 = replace(abbrev2, abbrev2 == "Cexp", "Crop\nexp"),
      abbrev2 = replace(abbrev2, abbrev2 == "Uexp", "Urban\nexp"),
      abbrev2 = replace(abbrev2, abbrev2 == "Vul", "Social\nvul"),
      abbrev2 = replace(abbrev2, abbrev2 == "VegHet", "Veg\nhet"),
      abbrev2 = replace(abbrev2, abbrev2 == "ForLoss", "Forest\nloss"),
      abbrev2 = replace(abbrev2, abbrev2 == "Precip", "Precip\nchange"),
      abbrev2 = replace(abbrev2, abbrev2 == "Temp", "Temp\nchange"),
      abbrev2 = replace(abbrev2, abbrev2 == "Health", "Health\naccess"),
      abbrev2 = replace(abbrev2, abbrev2 == "PA", "Prot\narea")
    )
  
  # nodes$abbrev2 = factor(nodes$abbrev2,
  #                        levels=order(nodes$abbrev2),
  #                        ordered=TRUE)
  #node_order = c("Uexp", "Urban", "Crop", "Cexp", "Biodiv", "Health", "Hunting", "Mining", "VegHet", "Livestock", "Vul", "PA", "Forest",  "ForLoss", "Precip", "Temp")
  node_order = c("Urban\nexp", "Urban", "Crop", "Crop\nexp", "Biodiv", "Health\naccess", "Hunting", "Mining", "Veg\nhet", "Livestock", "Social\nvul", "Prot\narea", "Forest",  "Forest\nloss", "Precip\nchange", "Temp\nchange")
  nodes = nodes[ match(node_order, nodes$abbrev2), ]
  
  # driver type
  nodes$driver_type = "Land use change/intensity"
  nodes$driver_type[ nodes$node %in% c("Biodiversity Intactness", "Cropland cover", "Forest cover", "Vegetation heterogeneity")] = "Ecosystem structure"
  nodes$driver_type[ nodes$node %in% c("Livestock density (log)", "Social vulnerability")] = "Socioeconomic"
  nodes$driver_type[ nodes$node %in% c("Temperature change", "Precipitation change")] = "Climate change"
  nodes$driver_type[ nodes$node %in% c("Built-up land", "Healthcare travel time (log)")] = "Detection"
  nodes$driver_type = factor(nodes$driver_type, levels=c("Detection", "Socioeconomic", "Ecosystem structure", "Climate change", "Land use change/intensity"), ordered=TRUE)
  
  # edges
  edges = sig_efs %>%
    dplyr::select(Disease, param) %>%
    dplyr::full_join(
      sig_efs %>% dplyr::select(Disease, param) %>% dplyr::rename(param2 = param)
    )
  edges = edges[ -which(edges$param == edges$param2), ]
  
  edges = edges %>%
    dplyr::group_by(param, param2) %>%
    dplyr::summarise(co_occ = n_distinct(Disease)) %>%
    #dplyr::arrange(desc(co_occ)) %>%
    dplyr::left_join(
      nodes %>% dplyr::mutate(from = 1:nrow(nodes)) %>% dplyr::rename(param = node) %>% dplyr::select(param, from)
    ) %>%
    dplyr::left_join(
      nodes %>% dplyr::mutate(to = 1:nrow(nodes)) %>% dplyr::rename(param2 = node) %>% dplyr::select(param2, to)
    ) %>%
    dplyr::mutate(`Driver\nco-occurrence\n(num. diseases)`=co_occ,
                  `Driver\nco-occurrence\n(prop. diseases)`=co_occ / num_dz,
                  `Co-occurrence`=co_occ)
  
  # one instance per edge
  ee = apply(edges[ , c("from", "to")], 1, function(x) paste(x[ order(x)], collapse=", "))
  edges = edges[ !duplicated(ee), ]
  
  # edges %>%
  #   ggplot() + 
  #   geom_tile(aes(x = param, y=param2, fill=co_occ)) + 
  #   theme_classic() + 
  #   scale_fill_viridis_c(option="magma")
  
  # combine into network
  graph = tbl_graph(nodes=nodes, edges=edges, directed = FALSE, node_key="node")
  return(graph)
  
}

# ------------ 1. causal strict -----------

# break down to causal strict
pp1 = pp %>%
  dplyr::filter(type == "Causal (strict)")

# all diseases
# graph = build_graph(pp1)
# pg1 = ggraph(graph, layout="circle") +
#   geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
#   geom_node_point(aes(size=Occurrence), pch=21, fill="grey98", color="grey75", alpha=0.95) +
#   geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
#   scale_color_gradientn(colors=pals::ocean.deep(n=30)[10:26], guide="none") +
#   scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[10:26]), guide="none") +
#   maptheme +
#   scale_edge_width(range = c(0.5, 3.5), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
#   scale_edge_alpha(range=c(0.3, 1), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
#   scale_size(range=c(3, 15), labels=c(0, 1, 5, 10, 15, 20), breaks=c(0, 1, 5, 10, 15, 20)) +
#   theme(legend.position="bottom", 
#         legend.text = element_text(size=10), 
#         legend.title = element_text(size=11),
#         #legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box = "vertical") +
#   expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
#   ggtitle("All diseases [n=31]") +
#   guides(size = guide_legend(order = 1))


# all diseases
graph = build_graph(pp1)
pg1 = ggraph(graph, layout="circle") +
  geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
  geom_node_point(aes(size=Occurrence), pch=21, fill="white", color="grey75", alpha=1) +
  geom_node_point(aes(size=Occurrence, fill=driver_type), pch=21, color="white", alpha=0.35) +
  geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
  scale_edge_colour_gradientn(colors=colorRampPalette(colors=c("grey70", "grey40"))(200), guide="none") +
  # scale_color_gradientn(colors=pals::ocean.deep(n=30)[10:26], guide="none") +
  # scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[10:26]), guide="none") +
  maptheme +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type", guide="none") + 
  scale_edge_width(range = c(0.5, 3.5), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
  scale_edge_alpha(range=c(0.3, 1), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
  scale_size(range=c(3, 15), labels=c(0, 1, 5, 10, 15, 20), breaks=c(0, 1, 5, 10, 15, 20)) +
  theme(legend.position="bottom", 
        legend.text = element_text(size=10), 
        legend.title = element_text(size=11),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical") +
  expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
  ggtitle("All diseases [n=31]") +
  guides(size = guide_legend(order = 1))

# zoonotic and vector-borne (n=17 dz)
graph = build_graph(
  pp1 %>% dplyr::filter(human_infection_route %in% c("Wildlife + vector")) %>% dplyr::mutate(Disease = as.vector(Disease))
)
# pg2 = ggraph(graph, layout="circle") +
#   geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
#   geom_node_point(aes(size=Occurrence), pch=21, fill="grey98", color="grey75", alpha=0.95) +
#   geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
#   scale_color_gradientn(colors=pals::ocean.deep(n=30)[10:26], guide="none") +
#   scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[10:26]), guide="none") +
#   maptheme +
#   scale_edge_width(range = c(0.5, 3.5), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
#   scale_edge_alpha(range=c(0.3, 1), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
#   scale_size(range=c(3, 15), labels=c(0, 1, 5, 10), breaks=c(0, 1, 5, 10)) +
#   theme(legend.position="bottom", 
#         legend.text = element_text(size=10), 
#         legend.title = element_blank(),
#         #legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box = "vertical") +
#   expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
#   ggtitle("Zoonotic (vector-borne) [n=16]") +
#   guides(size = guide_legend(order = 1))

pg2 = ggraph(graph, layout="circle") +
  geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
  geom_node_point(aes(size=Occurrence), pch=21, fill="white", color="grey75", alpha=1) +
  geom_node_point(aes(size=Occurrence, fill=driver_type), pch=21, color="white", alpha=0.35) +
  geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
  scale_edge_colour_gradientn(colors=colorRampPalette(colors=c("grey70", "grey40"))(200), guide="none") +
  maptheme +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type", guide="none") + 
  scale_edge_width(range = c(0.5, 3.5), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
  scale_edge_alpha(range=c(0.3, 1), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
  scale_size(range=c(3, 15), labels=c(0, 1, 5, 10), breaks=c(0, 1, 5, 10)) +
  theme(legend.position="bottom", 
        legend.text = element_text(size=10), 
        legend.title = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical") +
  expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
  ggtitle("Zoonotic (vector-borne) [n=17]") +
  guides(size = guide_legend(order = 1))


# zoonotic
graph = build_graph(
  pp1 %>% dplyr::filter(human_infection_route %in% c("Wildlife")) %>% dplyr::mutate(Disease = as.vector(Disease))
)
# pg3 = ggraph(graph, layout="circle") +
#   geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
#   geom_node_point(aes(size=Occurrence), pch=21, fill="grey98", color="grey75", alpha=0.95) +
#   geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
#   scale_color_gradientn(colors=pals::ocean.deep(n=30)[12:26], guide="none") +
#   scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[12:26]), guide="none") +
#   maptheme +
#   scale_edge_width(range = c(0.5, 2), labels=c(1,2), breaks=c(1,2)) +
#   scale_edge_alpha(range=c(0.38, 1), labels=c(1,2), breaks=c(1,2)) +
#   scale_size(range=c(3, 15), labels=c(0, 1, 2, 4), breaks=c(0, 1, 2, 4)) +
#   theme(legend.position="bottom", 
#         legend.text = element_text(size=10), 
#         legend.title = element_blank(),
#         #legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box = "vertical") +
#   expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
#   ggtitle("Zoonotic (direct) [n=9]") +
#   guides(size = guide_legend(order = 1))

pg3 = ggraph(graph, layout="circle") +
  geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
  geom_node_point(aes(size=Occurrence), pch=21, fill="white", color="grey75", alpha=1) +
  geom_node_point(aes(size=Occurrence, fill=driver_type), pch=21, color="white", alpha=0.35) +
  geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
  scale_edge_colour_gradientn(colors=colorRampPalette(colors=c("grey70", "grey40"))(200), guide="none") +
  maptheme +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type", guide="none") + 
  scale_edge_width(range = c(0.5, 2), labels=c(1,2,3), breaks=c(1,2,3)) +
  scale_edge_alpha(range=c(0.38, 1), labels=c(1,2,3), breaks=c(1,2,3)) +
  scale_size(range=c(3, 15), labels=c(0, 1, 2, 4), breaks=c(0, 1, 2, 4)) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical") +
  expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
  ggtitle("Zoonotic (direct) [n=10]") +
  guides(size = guide_legend(order = 1))

pc = gridExtra::grid.arrange(pg1, pg3, pg2, ncol=3)
# pc = ggpubr::as_ggplot(pc)  +
#   cowplot::draw_plot_label(label = c("a", "b", "c"), 
#                            fontface = "plain", size = 28, 
#                            x = c(0.01, 0.34, 0.69), y = c(0.96, 0.96, 0.96))

ggsave(pc, file="./output/plots/Figure4_CoOccurence_TopRanked_v3.jpg", device="jpg", units="in", width=16, height=7.2, dpi=600, scale=0.8)



# ------------ 1. causal broad -----------

pp1 = pp %>%
  dplyr::filter(type == "Causal (broad)")

# all diseases
# graph = build_graph(pp1)
# pg1 = ggraph(graph, layout="circle") +
#   geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
#   geom_node_point(aes(size=Occurrence), pch=21, fill="grey98", color="grey75", alpha=0.95) +
#   geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
#   scale_color_gradientn(colors=pals::ocean.deep(n=30)[10:26], guide="none") +
#   scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[10:26]), guide="none") +
#   maptheme +
#   scale_edge_width(range = c(0.5, 3.5), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
#   scale_edge_alpha(range=c(0.3, 1), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
#   scale_size(range=c(3, 15), labels=c(0, 1, 5, 10, 15, 20), breaks=c(0, 1, 5, 10, 15, 20)) +
#   theme(legend.position="bottom", 
#         legend.text = element_text(size=10), 
#         legend.title = element_text(size=11),
#         #legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box = "vertical") +
#   expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
#   ggtitle("All diseases [n=31]") +
#   guides(size = guide_legend(order = 1))

# all diseases
graph = build_graph(pp1)
pg1 = ggraph(graph, layout="circle") +
  geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
  geom_node_point(aes(size=Occurrence), pch=21, fill="white", color="grey75", alpha=1) +
  geom_node_point(aes(size=Occurrence, fill=driver_type), pch=21, color="white", alpha=0.35) +
  geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
  scale_edge_colour_gradientn(colors=colorRampPalette(colors=c("grey70", "grey40"))(200), guide="none") +
  # scale_color_gradientn(colors=pals::ocean.deep(n=30)[10:26], guide="none") +
  # scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[10:26]), guide="none") +
  maptheme +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type", guide="none") + 
  scale_edge_width(range = c(0.5, 3.5), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
  scale_edge_alpha(range=c(0.3, 1), labels=c(1, 5, 10), breaks=c(1, 5, 10)) +
  scale_size(range=c(3, 15), labels=c(0, 1, 5, 10, 15, 20), breaks=c(0, 1, 5, 10, 15, 20)) +
  theme(legend.position="bottom", 
        legend.text = element_text(size=10), 
        legend.title = element_text(size=11),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical") +
  expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
  ggtitle("All diseases [n=31]") +
  guides(size = guide_legend(order = 1))

# zoonotic and vector-borne (n=17 dz)
graph = build_graph(
  pp1 %>% dplyr::filter(human_infection_route %in% c("Wildlife + vector")) %>% dplyr::mutate(Disease = as.vector(Disease))
)
# pg2 = ggraph(graph, layout="circle") +
#   geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
#   geom_node_point(aes(size=Occurrence), pch=21, fill="grey98", color="grey75", alpha=0.95) +
#   geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
#   scale_color_gradientn(colors=pals::ocean.deep(n=30)[10:26], guide="none") +
#   scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[10:26]), guide="none") +
#   maptheme +
#   scale_edge_width(range = c(0.5, 3.5), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
#   scale_edge_alpha(range=c(0.3, 1), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
#   scale_size(range=c(3, 15), labels=c(0, 1, 5, 10), breaks=c(0, 1, 5, 10)) +
#   theme(legend.position="bottom", 
#         legend.text = element_text(size=10), 
#         legend.title = element_blank(),
#         #legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box = "vertical") +
#   expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
#   ggtitle("Zoonotic (vector-borne) [n=16]") +
#   guides(size = guide_legend(order = 1))

pg2 = ggraph(graph, layout="circle") +
  geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
  geom_node_point(aes(size=Occurrence), pch=21, fill="white", color="grey75", alpha=1) +
  geom_node_point(aes(size=Occurrence, fill=driver_type), pch=21, color="white", alpha=0.35) +
  geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
  scale_edge_colour_gradientn(colors=colorRampPalette(colors=c("grey70", "grey40"))(200), guide="none") +
  maptheme +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type", guide="none") + 
  scale_edge_width(range = c(0.5, 3.5), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
  scale_edge_alpha(range=c(0.3, 1), labels=c(2, 4, 6), breaks=c(2, 4, 6)) +
  scale_size(range=c(3, 15), labels=c(0, 1, 5, 10), breaks=c(0, 1, 5, 10)) +
  theme(legend.position="bottom", 
        legend.text = element_text(size=10), 
        legend.title = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical") +
  expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
  ggtitle("Zoonotic (vector-borne) [n=17]") +
  guides(size = guide_legend(order = 1))


# zoonotic
graph = build_graph(
  pp1 %>% dplyr::filter(human_infection_route %in% c("Wildlife")) %>% dplyr::mutate(Disease = as.vector(Disease))
)
# pg3 = ggraph(graph, layout="circle") +
#   geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
#   geom_node_point(aes(size=Occurrence), pch=21, fill="grey98", color="grey75", alpha=0.95) +
#   geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
#   scale_color_gradientn(colors=pals::ocean.deep(n=30)[12:26], guide="none") +
#   scale_edge_colour_gradientn(colors=rev(pals::ocean.deep(n=30)[12:26]), guide="none") +
#   maptheme +
#   scale_edge_width(range = c(0.5, 2), labels=c(1,2), breaks=c(1,2)) +
#   scale_edge_alpha(range=c(0.38, 1), labels=c(1,2), breaks=c(1,2)) +
#   scale_size(range=c(3, 15), labels=c(0, 1, 2, 4), breaks=c(0, 1, 2, 4)) +
#   theme(legend.position="bottom", 
#         legend.text = element_text(size=10), 
#         legend.title = element_blank(),
#         #legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box = "vertical") +
#   expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
#   ggtitle("Zoonotic (direct) [n=9]") +
#   guides(size = guide_legend(order = 1))

pg3 = ggraph(graph, layout="circle") +
  geom_edge_link(aes(alpha=`Co-occurrence`, width=`Co-occurrence`, color=co_occ)) +
  geom_node_point(aes(size=Occurrence), pch=21, fill="white", color="grey75", alpha=1) +
  geom_node_point(aes(size=Occurrence, fill=driver_type), pch=21, color="white", alpha=0.35) +
  geom_node_text(aes(label=abbrev2), size=3.5, color="black") +
  scale_edge_colour_gradientn(colors=colorRampPalette(colors=c("grey70", "grey40"))(200), guide="none") +
  maptheme +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 5, 2, 4, 3)], name="Driver type", guide="none") + 
  scale_edge_width(range = c(0.5, 2), labels=c(1,2,3), breaks=c(1,2,3)) +
  scale_edge_alpha(range=c(0.38, 1), labels=c(1,2,3), breaks=c(1,2,3)) +
  scale_size(range=c(3, 15), labels=c(0, 1, 2, 4), breaks=c(0, 1, 2, 4)) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical") +
  expand_limits(x=c(-1.05, 1.05), y=c(-1.05, 1.05)) +
  ggtitle("Zoonotic (direct) [n=10]") +
  guides(size = guide_legend(order = 1))

pc = gridExtra::grid.arrange(pg1, pg3, pg2, ncol=3)
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           fontface = "plain", size = 28, 
                           x = c(0.01, 0.34, 0.69), y = c(0.96, 0.96, 0.96))

ggsave(pc, file="./output/plots/Figure4_CoOccurence_MajorityRule_v3.jpg", device="jpg", units="in", width=16, height=7.2, dpi=600, scale=0.85)






# =============== observed covariance between drivers worldwide ==================

# background points
bg = read.csv("./output/model_df/global_background/globalbackgroundpoints_df_wp.csv") %>%
  dplyr::select(
    -all_of(c("livestock_all", "precip", "tmean", "health_travel", "evi_coefvar"))
  )

# get region
rr = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
bgx = st_as_sf(bg, coords = c("Longitude", "Latitude"), crs=crs(rr))
ii = sf::st_intersects(bgx, rr)
ii = unlist(lapply(ii, function(x) ifelse(length(x)==0, NA, x)))
bg$iso3 = rr$iso3[ ii ]
bg$region = countrycode::countrycode(bg$iso3, origin="iso3c", destination = "region")

# coastline
coast = rnaturalearth::ne_coastline(scale=50, returnclass = "sf") 
ext = raster::extent(c(xmin=-180, xmax=180, ymin=-60, ymax=83))
coast$feature = 1:nrow(coast)
cent = as.data.frame(st_coordinates(st_centroid(coast)))
coast = coast[ cent$Y > -60, ]

# map data
rr$region2 = countrycode::countrycode(rr$iso3, origin="iso3c", destination = "region")
rr$region2[ rr$iso3 == "MYT" & !is.na(rr$iso3) ] = "Sub-Saharan Africa"
rr$region2[ rr$iso3 == "REU" & !is.na(rr$iso3) ] = "Sub-Saharan Africa"

# transform
robinson = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
rr = st_transform(rr, robinson)
coast = st_transform(coast, robinson)
bgt = st_as_sf(bg, coords = c("Longitude", "Latitude"), crs=crs(rr))

map_plot = rr %>%
  dplyr::filter(!region2 %in% c("Middle East & North Africa", "Europe & Central Asia")) %>%
  ggplot() + 
  geom_sf(aes(fill=region2), color=NA) +
  geom_sf(data=coast, color="grey50", fill=NA, linewidth=0.25) + 
  #geom_sf(data=bgt, color="grey30", size=0.01, alpha=0.1) +
  theme_void() +
  MetBrewer::scale_fill_met_d("Archambault", na.value="grey95") +
  theme(legend.position = "none", 
        panel.grid.major = element_line(color="grey92", linewidth=0.3))
#ggsave(map_plot, filename="./output/plots/nexusmap1.jpg", device="jpg", units="in", width=8, height=4, dpi=600)

# viz functions
createCorrMat = function(x, cols_to_exclude=NA){
  
  x = x %>% dplyr::select(-region, -iso3)
  
  if(!is.na(cols_to_exclude[1])){
    x = x %>% dplyr::select(-all_of(cols_to_exclude))
  }
  
  # subset and scale
  x = x[ , 6:ncol(x)]
  x = apply(x, 2, scale)
  x = as.data.frame(x)
  x = x %>% dplyr::filter(complete.cases(.))
  
  cm = cor(x)
  cm = corrplot::corrplot(cm)
  #print(cm)
  
  return(cm)
}

# create network from corr matrix
createNetwork = function(x, vifs){
  
  # nodes
  nodes = data.frame(
    node = unique(x$xName),
    id = 1:length(unique(x$xName)),
    region = x$region[1]
  ) %>%
    left_join(read.csv("./scripts/04_results/driver_abbrevs.csv"), by=c("node"="driver_code"))
  nodes$abbrev3[ nodes$abbrev3 == "Crop\\nexp." ] = "Crop\nexp."
  nodes$abbrev3[ nodes$abbrev3 == "Veg\\nhet" ] = "Veg.\nhet."
  nodes$abbrev3[ nodes$abbrev3 == "Forest\\nLoss" ] = "Forest\nloss"
  nodes$abbrev3[ nodes$abbrev3 == "Precip\\nchange" ] = "Precip\nchange"
  nodes$abbrev3[ nodes$abbrev3 == "Prot.\\narea" ] = "Prot.\narea"
  nodes$abbrev3[ nodes$abbrev3 == "Temp\\nchange" ] = "Temp\nchange"
  nodes$abbrev3[ nodes$abbrev3 == "Urban\\nexp." ] = "Urban\nexp."
  nodes$abbrev3[ nodes$abbrev3 == "Health\\ntravel" ] = "Health\naccess"
  nodes$abbrev3[ nodes$abbrev3 == "Social\\nvul" ] = "Social\nvul."
  # node_order = c("Urban", "Crop", "Crop\nexp.", "Biodiv", "Health\naccess", "Mining", "Veg.\nhet.", "Livestock", "Social\nvul.", "Prot.\narea", "Forest",  "Forest\nloss", "Precip\nchange", "Temp\nchange")
  # nodes = nodes[ match(node_order, nodes$abbrev3), ]
  
  nodes = left_join(nodes, vifs, by=c("node"="Variables"))
  
  # edges
  edges = x %>%
    dplyr::filter(xName != yName) %>%
    dplyr::select(
      xName, yName, corr
    ) %>%
    dplyr::left_join(
      nodes[ , c("node", "id")], by=c("xName" = "node")
    ) %>%
    dplyr::rename("from"=id) %>%
    dplyr::left_join(
      nodes[ , c("node", "id")], by=c("yName" = "node")
    ) %>%
    dplyr::rename("to"=id) %>%
    dplyr::mutate(ac = abs(corr), ac2 = ac, region=nodes$region[1])
  
  # add corr grp
  edges$corr2 = "Strong negative (<-0.7)"
  edges$corr2[ abs(edges$corr) > -0.7 ] = c("Negative (-0.7 to -0.5)")
  edges$corr2[ abs(edges$corr) > -0.5 ] = c("Weak negative (-0.5 to -0.3)")
  edges$corr2[ abs(edges$corr) > -0.3 ] = c("No correlation (-0.3 to 0.3)")
  edges$corr2[ abs(edges$corr) > 0.3 ] = c("Weak positive (0.3 to 0.5)")
  edges$corr2[ abs(edges$corr) > 0.5 ] = c("Positive (0.5 to 0.7)")
  edges$corr2[ abs(edges$corr) > 0.7 ] = c("Strong positive (>0.7)")
  edges$corr2 = factor(edges$corr2, levels=c("Strong negative (<-0.7)",
                                             "Negative (-0.7 to -0.5)",
                                             "Weak negative (-0.5 to -0.3)",
                                             "No correlation (-0.3 to 0.3)",
                                             "Weak positive (0.3 to 0.5)",
                                             "Positive (0.5 to 0.7)",
                                             "Strong positive (>0.7)"), ordered=TRUE)
  
  # one instance for each edge
  ee = apply(edges[ , c("xName", "yName")], 1, function(x) paste(x[ order(x)], collapse=", "))
  edges = edges[ !duplicated(ee), ]
  
  edges
  
  return(
    tbl_graph(nodes=nodes, edges=edges, directed = FALSE, node_key="node")
  )
}

# 
corrMatPlot = function(region, lims = c(-0.9, 0.9), cols_to_exclude="hunting", plot_vif=FALSE, facet=TRUE){
  
  if(region == "Global"){ 
    bgx = bg 
    bgx[ , 6:21 ] = apply(bgx[ , 6:21 ], 2, scale)
  } else{ 
    bgx = bg
    bgx[ , 6:21 ] = apply(bgx[ , 6:21 ], 2, scale)
    bgx = bgx[ bgx$region == region, ]
  }
  
  m1 = createCorrMat(bgx, cols_to_exclude = cols_to_exclude)$corrPos 
  m1$region = region
  print(sum(m1$corr>0.5 | m1$corr < -0.5) / nrow(m1))
  print(sum(m1$corr>0.7 | m1$corr < -0.7) / nrow(m1))
  
  # vifs
  vifs = data.frame(usdm::vif(bgx %>% dplyr::select(-all_of(c(cols_to_exclude, "record_id", "Longitude", "Latitude", "ADMcode", "presence", "iso3", "region")))))
  
  # mean values
  means = apply(bgx[ , 6:21], 2, function(x) mean(x, na.rm=TRUE))
  means = as.data.frame(means); means$Variables = row.names(means); row.names(means)=c()
  means$avg = ifelse(means$means < 0, "Below", "Above" )
  #means$avg[ which(abs(means$means) < 0.1 ) ] = "Avg"
  if(region == "Global"){ means$means = 0; means$avg = "Avg" }
  vifs = left_join(vifs, means)
  print(means$means)
  
  m1n = createNetwork(m1, vifs = vifs)
  
  # lims = max(abs(m1$corr[ m1$xName != m1$yName ]))
  # lims = c(-lims, lims)
  if(plot_vif == FALSE){
    
    p1 = ggraph(m1n, layout="circle") + 
      geom_edge_link(aes(alpha=ac, width=ac2, color=corr)) + 
      geom_node_point(size=19, color="grey70", fill="grey98", pch=21) + 
      geom_node_text(aes(label=abbrev3), size=3.5) +
      scale_edge_colour_gradientn(colors=pals::brewer.brbg(200), limits=lims, name="Spatial\ndriver\ncorrelation", breaks=c(-0.8, -0.4, 0, 0.4, 0.8)) +
      maptheme +
      scale_edge_width(range = c(0.1, 4), breaks = c(0.2, 0.4, 0.6, 0.8), limits=c(0, 0.9), name="Corr.\n(abs)", guide="none") +
      scale_edge_alpha(range=c(0.5, 0.8), breaks = c(0.2, 0.4, 0.6, 0.8), limits=c(0, 0.9), guide = "none") + 
      facet_grid(.~region) +
      theme(strip.text = element_text(hjust=0.5, size=21),
            legend.title = element_text(size=15), 
            strip.background = element_blank(),
            legend.text = element_text(size=19)) +
      expand_limits(x=c(-1.06, 1.06), y=c(-1.06, 1.06))+ 
      guides(edge_alpha = FALSE)
    
  } else{
    
    # p1 = ggraph(m1n, layout="circle") +
    #   geom_edge_link(aes(alpha=ac, width=ac2, color=corr)) +
    #   geom_node_point(aes(size=means), color="grey70", fill="grey98", pch=21) +
    #   geom_node_text(aes(label=abbrev3), size=3.9) +
    #   scale_edge_colour_gradientn(colors=pals::brewer.brbg(200), limits=lims, name="Spatial\ndriver\ncorrelation", breaks=c(-0.8, -0.4, 0, 0.4, 0.8)) +
    #   maptheme +
    #   scale_edge_width(range = c(0.1, 4), breaks = c(0.2, 0.4, 0.6, 0.8), limits=c(0, 0.9), name="Corr.\n(abs)", guide="none") +
    #   scale_edge_alpha(range=c(0.5, 0.8), breaks = c(0.2, 0.4, 0.6, 0.8), limits=c(0, 0.9), guide = "none") +
    #   scale_size(range=c(4, 28), limits=c(-1.2, 1.2), guide="none") +
    #   facet_grid(.~region) +
    #   theme(strip.text = element_text(hjust=0.5, size=21),
    #         legend.title = element_text(size=15), 
    #         strip.background = element_blank(),
    #         legend.text = element_text(size=19)) +
    #   expand_limits(x=c(-1.06, 1.06), y=c(-1.06, 1.06))+
    #   guides(edge_alpha = FALSE, size=FALSE)
    
    p1 = ggraph(m1n, layout="circle") +
      geom_edge_link(aes(alpha=ac, width=ac2, color=corr)) +
      geom_node_point(aes(size=means), color="white", fill="grey98", pch=21) +
      geom_node_point(aes(size=means, color=avg), fill=NA, pch=21, alpha=0.3) +
      geom_node_text(aes(label=abbrev3, color=avg), size=3.9) +
      scale_edge_colour_gradientn(colors=pals::brewer.brbg(200), limits=lims, name="Spatial\ndriver\ncorrelation", breaks=c(-0.8, -0.4, 0, 0.4, 0.8)) +
      maptheme +
      scale_edge_width(range = c(0.1, 4), breaks = c(0.2, 0.4, 0.6, 0.8), limits=c(0, 0.9), name="Corr.\n(abs)", guide="none") +
      scale_edge_alpha(range=c(0.5, 0.8), breaks = c(0.2, 0.4, 0.6, 0.8), limits=c(0, 0.9), guide = "none") +
      scale_size(range=c(3, 23), limits=c(-1.2, 1.2), guide="none") +
      scale_color_manual(values=c("Avg"="black", "Above"="darkred", "Below"="darkblue"), guide="none") +
      facet_grid(.~region) +
      theme(strip.text = element_text(hjust=0.5, size=17),
            legend.title = element_text(size=15), 
            strip.background = element_blank(),
            legend.text = element_text(size=19)) +
      expand_limits(x=c(-1.06, 1.06), y=c(-1.06, 1.06))+
      guides(edge_alpha = FALSE, size=FALSE)
    
    
    if(facet == FALSE){
      p1 = ggraph(m1n, layout="circle") +
        geom_edge_link(aes(alpha=ac, width=ac2, color=corr)) +
        geom_node_point(aes(size=means), color="grey70", fill="grey98", pch=21) +
        geom_node_text(aes(label=abbrev3), size=3.9) +
        scale_edge_colour_gradientn(colors=pals::brewer.brbg(200), limits=lims, name="Spatial\ndriver\ncorrelation", breaks=c(-0.5, -0.25, 0, 0.25, 0.5)) +
        maptheme +
        scale_edge_width(range = c(0.1, 2.5), breaks = c(0.1, 0.25, 0.5), limits=c(0, 0.55), name="Corr.\n(abs)", guide="none") +
        scale_edge_alpha(range=c(0.5, 0.8), breaks = c(0.1, 0.25, 0.5), limits=c(0, 0.55), guide = "none") +
        scale_size(range=c(3, 23), limits=c(-1.2, 1.2), guide="none") +
        theme(strip.text = element_text(hjust=0.5, size=21),
              legend.title = element_text(size=15), 
              strip.background = element_blank(),
              legend.text = element_text(size=12), legend.position = "left") +
        expand_limits(x=c(-1.06, 1.06), y=c(-1.06, 1.06))+
        guides(edge_alpha = FALSE, size=FALSE)
    }
    
  }
  
  return(p1)
}


# # driver correlations
# 
# p1 = corrMatPlot("Global")
# p2 = corrMatPlot("North America")
# p3 = corrMatPlot("Latin America & Caribbean")
# p4 = corrMatPlot("South Asia")
# p5 = corrMatPlot("Sub-Saharan Africa")
# p6 = corrMatPlot("East Asia & Pacific")
# 
# library(patchwork)
# comb = p1 + p2 + p3 + p4 + p5 + p6
# comb = comb + plot_layout(guides = "collect", ncol=2)
# comb
# 
# ggsave(comb,
#        file="./output/plots/SuppFigure_DriverCorrelations.jpg",
#        device="jpg", units="in", width=12.5, height=17, dpi=300)


# driver correlations with node size by VIF

p1 = corrMatPlot("Global", cols_to_exclude = c("hunting", "urban_expansion"), plot_vif = TRUE)
p2 = corrMatPlot("North America", cols_to_exclude = c("hunting", "urban_expansion"), plot_vif = TRUE)
p3 = corrMatPlot("Latin America & Caribbean", cols_to_exclude = c("hunting", "urban_expansion"), plot_vif = TRUE)
p4 = corrMatPlot("South Asia", cols_to_exclude = c("hunting", "urban_expansion"), plot_vif = TRUE)
p5 = corrMatPlot("Sub-Saharan Africa", cols_to_exclude = c("hunting", "urban_expansion"), plot_vif = TRUE)
p6 = corrMatPlot("East Asia & Pacific", cols_to_exclude = c("hunting", "urban_expansion"), plot_vif = TRUE)

#horizontal
library(patchwork)
comb = p1 + p2 + p3 + p4 + p5 + p6
comb = comb + plot_layout(guides = "collect", ncol=3)
ggsave(comb,
       file="./output/plots/SuppFigure_DriverCorrelations_MeanVals_horiz.jpg",
       device="jpg", units="in", height=12, width=19.3, dpi=300, scale=0.85)












# 
# 
# 
# 
# # ============== do drivers co-occur more often than expected by chance given their prevalence? ==============
# 
# params = list.files("./output/model_outputs/disease_params/", pattern=".csv", full.names=TRUE)
# params = params[ -grep("infocriteria", params)]
# pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
# pp$mean = round(pp$mean, 3)
# pp$lower = round(pp$lower, 3)
# pp$upper = round(pp$upper, 3)
# pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0
# 
# # add abbrevs
# pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
#   dplyr::mutate(Disease = abbrev2)
# 
# pp = pp %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
#   dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
#                 param = replace(param, param == "livestock_log", "Livestock density (log)"),
#                 param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
#                 param = replace(param, param == "hunting", "Hunting pressure"),
#                 param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
#                 param = replace(param, param == "forest_cover", "Forest cover"),
#                 param = replace(param, param == "forest_loss", "Forest loss"),
#                 param = replace(param, param == "crop_expansion", "Cropland expansion"),
#                 param = replace(param, param == "evi_dissimilarity", "Vegetation heterogeneity"),
#                 param = replace(param, param == "evi_coefvar", "EVI variance"),
#                 param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
#                 param = replace(param, param == "urban_expansion", "Urban expansion"),
#                 param = replace(param, param == "urban_cover", "Built-up land"),
#                 param = replace(param, param == "crop_cover", "Cropland cover"),
#                 param = replace(param, param == "mining", "Mining"),
#                 param = replace(param, param == "tmean_change", "Temperature change"),
#                 param = replace(param, param == "precip_change", "Precipitation change"),
#                 param = replace(param, param == "social_vulnerability", "Social vulnerability"),
#                 param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
#                 param = replace(param, param == "protected_areas", "Protected area coverage"),
#                 param = replace(param, param == "popdens_log", "Population density (log)")) 
# 
# # order by n points
# np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
# pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)
# 
# # break down to causal broad
# pp = pp %>%
#   dplyr::filter(type == "Causal (broad)")
# 
# # exclude any with no shared drivers (mining)
# # pp = pp %>% 
# #   dplyr::filter(param != "Mining")
# 
# # driver pairs
# dp = expand.grid(unique(pp$param), unique(pp$param)) %>%
#   as.data.frame() 
# dp = dp[ -which(dp$Var1 == dp$Var2), ]
# row.names(dp) = c()
# 
# # odds ratios dataframe
# result = data.frame()
# 
# # tabulate for each driver pair 
# for(i in 1:nrow(dp)){
#   
#   p_a = as.vector(dp$Var1[i])
#   p_b = as.vector(dp$Var2[i])
#   
#   print(paste(p_a, p_b, sep=" + "))
#   
#   tab_i = data.frame()
#   for(d in unique(pp$Disease)){
# 
#     # entire pop all diseases
#     # pp_d = pp %>% dplyr::filter(Disease == d) %>% dplyr::filter(sig == TRUE)
#     # res_d = data.frame(disease = d, 
#     #                    p_a = p_a, p_b = p_b,
#     #                    a_present = as.numeric(p_a %in% pp_d$param),
#     #                    b_present = as.numeric(p_b %in% pp_d$param))
#     
#     # only both tested
#     pp_d = pp %>% dplyr::filter(Disease == d) 
#     if(!(p_a %in% pp_d$param & p_b %in% pp_d$param)){ 
#       next 
#     } else{
#       pp_d = pp_d %>% dplyr::filter(sig == TRUE)
#       res_d = data.frame(disease = d, 
#                          p_a = p_a, p_b = p_b,
#                          a_present = as.numeric(p_a %in% pp_d$param),
#                          b_present = as.numeric(p_b %in% pp_d$param))
#       tab_i = rbind(tab_i, res_d)
#       }
# 
#   }
#   
#   # contingency table
#   # 
#   # # co occurrence
#   # a1_b1 = sum(tab_i$a_present & tab_i$b_present)
#   # a1_b0 = sum(tab_i$a_present & (1-tab_i$b_present))
#   # a0_b1 = sum((1-tab_i$a_present) & tab_i$b_present)
#   # a0_b0 = sum((1-tab_i$a_present) & (1-tab_i$b_present))
#   # 
#   # # independent occurrence
#   # a1 = sum(tab_i$a_present)
#   # b1 = sum(tab_i$b_present)
#   # a0 = sum((1-tab_i$a_present))
#   # b0 = sum((1-tab_i$b_present))
#   # 
#   # # total n 
#   # tot = n_distinct(pp$Disease)
#   # 
#   # # tests
#   # # if( (a1_b1 + a1_b0 + a0_b1 + a0_b0) == tot ) print("Co-occ correct")
#   # # if( (a1 + a0) == tot ) print("Ind A correct")
#   # # if( (b1 + b0) == tot ) print("Ind B correct")
#   # 
#   # # calculate ORs
#   # # https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1182
#   # 
#   # # or1 - odds of A being present when B is / odds of A being present regardless of B
#   # or_AB = (a1_b1 / a0_b1) / (a1 / a0)
#   # 
#   # # or2 - odds of B being present when A is / odds of B being present regardless of A
#   # or_BA = (a1_b1 / a1_b0) / (b1 / b0)
#   # 
#   # # odds of A being present when B is / odds of A being present when B is absent
#   # or_SYM = (a1_b1 / a0_b1) / (a1_b0 / a0_b0)
#   # 
#   # # save output
#   # result = rbind(result,
#   #                data.frame(p_a = p_a,
#   #                           p_b = p_b,
#   #                           or_AB = or_AB,
#   #                           or_BA = or_BA, 
#   #                           or_SYM = or_SYM)
#   # )
#   
#   # only where at least 3 times tested in combination
#   if(nrow(tab_i) > 5){
#     
#     # glm
#     or_dat = tab_i[ , c("a_present", "b_present")]
#     or_AB = sppairs::or.glm(or_dat, complex=TRUE)
#     res_i = data.frame(
#       p_a = p_a,
#       p_b = p_b,
#       or_AB = or_AB[5],
#       pval = or_AB[4],
#       n = nrow(tab_i), 
#       prop_both_present = sum(or_dat$a_present == 1 & or_dat$b_present == 1) / nrow(or_dat)
#     )
#     row.names(res_i) = c()
#     result = rbind(result, res_i)
#     
#   }
#   
# }
# 
# result %>% 
#   dplyr::arrange(prop_both_present)
# 
# 
# 
# 
# 
# 
# 
# # ============== do drivers co-occur more often than expected by chance given their prevalence? ==============
# 
# params = list.files("./output/model_outputs/disease_params/", pattern=".csv", full.names=TRUE)
# params = params[ -grep("infocriteria", params)]
# pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
# pp$mean = round(pp$mean, 3)
# pp$lower = round(pp$lower, 3)
# pp$upper = round(pp$upper, 3)
# pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0
# 
# # add abbrevs
# pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
#   dplyr::mutate(Disease = abbrev2)
# 
# pp = pp %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
#   dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
#                 param = replace(param, param == "livestock_log", "Livestock density (log)"),
#                 param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
#                 param = replace(param, param == "hunting", "Hunting pressure"),
#                 param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
#                 param = replace(param, param == "forest_cover", "Forest cover"),
#                 param = replace(param, param == "forest_loss", "Forest loss"),
#                 param = replace(param, param == "crop_expansion", "Cropland expansion"),
#                 param = replace(param, param == "evi_dissimilarity", "Vegetation heterogeneity"),
#                 param = replace(param, param == "evi_coefvar", "EVI variance"),
#                 param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
#                 param = replace(param, param == "urban_expansion", "Urban expansion"),
#                 param = replace(param, param == "urban_cover", "Built-up land"),
#                 param = replace(param, param == "crop_cover", "Cropland cover"),
#                 param = replace(param, param == "mining", "Mining"),
#                 param = replace(param, param == "tmean_change", "Temperature change"),
#                 param = replace(param, param == "precip_change", "Precipitation change"),
#                 param = replace(param, param == "social_vulnerability", "Social vulnerability"),
#                 param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
#                 param = replace(param, param == "protected_areas", "Protected area coverage"),
#                 param = replace(param, param == "popdens_log", "Population density (log)")) 
# 
# # order by n points
# np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
# pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)
# 
# # break down to causal broad
# pp = pp %>%
#   dplyr::filter(type == "Causal (broad)")
# 
# # exclude any with no shared drivers (mining)
# # pp = pp %>% 
# #   dplyr::filter(param != "Mining")
# 
# # driver pairs
# dp = expand.grid(unique(pp$param), unique(pp$param)) %>%
#   as.data.frame() 
# dp = dp[ -which(dp$Var1 == dp$Var2), ]
# row.names(dp) = c()
# 
# # odds ratios dataframe
# result = data.frame()
# 
# # tabulate for each driver pair 
# for(i in 1:nrow(dp)){
#   
#   p_a = as.vector(dp$Var1[i])
#   p_b = as.vector(dp$Var2[i])
#   
#   print(paste(p_a, p_b, sep=" + "))
#   
#   tab_i = data.frame()
#   for(d in unique(pp$Disease)){
#     
#     # entire pop all diseases
#     # pp_d = pp %>% dplyr::filter(Disease == d) %>% dplyr::filter(sig == TRUE)
#     # res_d = data.frame(disease = d, 
#     #                    p_a = p_a, p_b = p_b,
#     #                    a_present = as.numeric(p_a %in% pp_d$param),
#     #                    b_present = as.numeric(p_b %in% pp_d$param))
#     
#     # only both tested
#     pp_d = pp %>% dplyr::filter(Disease == d) 
#     if(!(p_a %in% pp_d$param & p_b %in% pp_d$param)){ 
#       next 
#     } else{
#       pp_d = pp_d %>% dplyr::filter(sig == TRUE)
#       res_d = data.frame(disease = d, 
#                          p_a = p_a, p_b = p_b,
#                          a_present = as.numeric(p_a %in% pp_d$param),
#                          b_present = as.numeric(p_b %in% pp_d$param),
#                          a_present_positive = as.numeric(pp_d$mean[ pp_d$param == p_a] > 0),
#                          a_present_negative = as.numeric(pp_d$mean[ pp_d$param == p_a] < 0),
#                          b_present_positive = as.numeric(pp_d$mean[ pp_d$param == p_b] > 0) ,
#                          b_present_negative = as.numeric(pp_d$mean[ pp_d$param == p_b] < 0))
#       tab_i = rbind(tab_i, res_d)
#     }
#     
#   }
#   
#   # contingency table
#   # 
#   # # co occurrence
#   # a1_b1 = sum(tab_i$a_present & tab_i$b_present)
#   # a1_b0 = sum(tab_i$a_present & (1-tab_i$b_present))
#   # a0_b1 = sum((1-tab_i$a_present) & tab_i$b_present)
#   # a0_b0 = sum((1-tab_i$a_present) & (1-tab_i$b_present))
#   # 
#   # # independent occurrence
#   # a1 = sum(tab_i$a_present)
#   # b1 = sum(tab_i$b_present)
#   # a0 = sum((1-tab_i$a_present))
#   # b0 = sum((1-tab_i$b_present))
#   # 
#   # # total n 
#   # tot = n_distinct(pp$Disease)
#   # 
#   # # tests
#   # # if( (a1_b1 + a1_b0 + a0_b1 + a0_b0) == tot ) print("Co-occ correct")
#   # # if( (a1 + a0) == tot ) print("Ind A correct")
#   # # if( (b1 + b0) == tot ) print("Ind B correct")
#   # 
#   # # calculate ORs
#   # # https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1182
#   # 
#   # # or1 - odds of A being present when B is / odds of A being present regardless of B
#   # or_AB = (a1_b1 / a0_b1) / (a1 / a0)
#   # 
#   # # or2 - odds of B being present when A is / odds of B being present regardless of A
#   # or_BA = (a1_b1 / a1_b0) / (b1 / b0)
#   # 
#   # # odds of A being present when B is / odds of A being present when B is absent
#   # or_SYM = (a1_b1 / a0_b1) / (a1_b0 / a0_b0)
#   # 
#   # # save output
#   # result = rbind(result,
#   #                data.frame(p_a = p_a,
#   #                           p_b = p_b,
#   #                           or_AB = or_AB,
#   #                           or_BA = or_BA, 
#   #                           or_SYM = or_SYM)
#   # )
#   
#   # only where at least 3 times tested in combination
#   if(nrow(tab_i) > 5){
#     
#     # glm
#     or_dat = tab_i[ , c("a_present", "b_present")]
#     or_AB = sppairs::or.glm(or_dat, complex=TRUE)
#     res_i = data.frame(
#       p_a = p_a,
#       p_b = p_b,
#       or_AB = or_AB[5],
#       pval = or_AB[4],
#       n = nrow(tab_i)
#     )
#     row.names(res_i) = c()
#     result = rbind(result, res_i)
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
# 
# # 1. nodes
# nodes = sig_efs %>%
#   dplyr::group_by(param) %>%
#   dplyr::summarise(n_models = n_distinct(Disease)) %>%
#   dplyr::arrange(desc(n_models)) %>%
#   dplyr::rename(node=1)
# 
# # add abbrevs for visualisation
# nodes = left_join(nodes, read.csv("./scripts/04_results/driver_abbrevs.csv"), by=c("node"="driver"))
# 
# # edges
# edges = sig_efs %>%
#   dplyr::select(Disease, param) %>%
#   dplyr::full_join(
#     sig_efs %>% dplyr::select(Disease, param) %>% dplyr::rename(param2 = param)
#   )
# edges = edges[ -which(edges$param == edges$param2), ]
# 
# edges = edges %>%
#   dplyr::group_by(param, param2) %>%
#   dplyr::summarise(co_occ = n_distinct(Disease)) %>%
#   dplyr::arrange(desc(co_occ)) %>%
#   dplyr::left_join(
#     nodes %>% dplyr::mutate(from = 1:nrow(nodes)) %>% dplyr::rename(param = node) %>% dplyr::select(param, from)
#   ) %>%
#   dplyr::left_join(
#     nodes %>% dplyr::mutate(to = 1:nrow(nodes)) %>% dplyr::rename(param2 = node) %>% dplyr::select(param2, to)
#   ) %>%
#   dplyr::mutate(`Driver\nco-occurrence\n(num. diseases)`=co_occ)
# 
# e1 = edges %>% 
#   left_join(nodes[ c("node", "n_models", "abbrev2")] %>% dplyr::rename("param"=1, "n_models1"=2, "abb1"=3)) %>%
#   left_join(nodes[ c("node", "n_models", "abbrev2")] %>% dplyr::rename("param2"=1, "n_models2"=2, "abb2"=3)) 
# 
# e1$param_pair = apply(e1[ , c("abb1", "abb2")], 1, function(x) paste(x[ order(x)], collapse=" + "))
# e1$n_models_mean = (e1$n_models1 + e1$n_models2)/2
# e1 = e1[ !duplicated(e1$param_pair), ]
# 
# ggplot(e1) + 
#   geom_smooth(aes(x=n_models_mean, y=co_occ), method="glm", method.args = list(family = "poisson")) +
#   geom_text(aes(x=n_models_mean, y=co_occ, label=param_pair), position=position_jitter(), size=2.5, alpha=.8) +
#   theme_classic() + 
#   xlab("N models with individual significant effect of driver (mean a + b)") + 
#   ylab("N models with co-occurring significant effects")
# 
# 
# 
# 
# 
# 
# 
# params = list.files("./output/model_outputs/disease_params/", pattern=".csv", full.names=TRUE)
# params = params[ -grep("infocriteria", params)]
# pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
# pp$mean = round(pp$mean, 3)
# pp$lower = round(pp$lower, 3)
# pp$upper = round(pp$upper, 3)
# pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0
# 
# # add abbrevs
# pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
#   dplyr::mutate(Disease = abbrev2)
# 
# pp = pp %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
#   dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
#                 param = replace(param, param == "livestock_log", "Livestock density (log)"),
#                 param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
#                 param = replace(param, param == "hunting", "Hunting pressure"),
#                 param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
#                 param = replace(param, param == "forest_cover", "Forest cover"),
#                 param = replace(param, param == "forest_loss", "Forest loss"),
#                 param = replace(param, param == "crop_expansion", "Cropland expansion"),
#                 param = replace(param, param == "evi_dissimilarity", "Vegetation heterogeneity"),
#                 param = replace(param, param == "evi_coefvar", "EVI variance"),
#                 param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
#                 param = replace(param, param == "urban_expansion", "Urban expansion"),
#                 param = replace(param, param == "urban_cover", "Built-up land"),
#                 param = replace(param, param == "crop_cover", "Cropland cover"),
#                 param = replace(param, param == "mining", "Mining"),
#                 param = replace(param, param == "tmean_change", "Temperature change"),
#                 param = replace(param, param == "precip_change", "Precipitation change"),
#                 param = replace(param, param == "social_vulnerability", "Social vulnerability"),
#                 param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
#                 param = replace(param, param == "protected_areas", "Protected area coverage"),
#                 param = replace(param, param == "popdens_log", "Population density (log)")) 
# 
# # order by n points
# np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
# pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)
# 
# # break down to causal broad
# pp = pp %>%
#   dplyr::filter(type == "Causal (broad)")
# 
# # exclude any with no shared drivers (mining)
# pp = pp %>% 
#   dplyr::filter(param != "Mining")
# 
# # significant effects only 
# sig_efs = pp 
# # 
# sig_efs = sig_efs %>%
#   left_join(read.csv("./scripts/04_results/driver_abbrevs.csv") %>% 
#               dplyr::rename("param"=2) %>%
#               dplyr::select(param, abbrev2) %>%
#               dplyr::rename("param2"=abbrev2))
# 
# # params
# px_for_rand = sig_efs$param2
# 
# # generate parameter pairs
# ppairs = c()
# for(i in 1:1000){
#   pi = sample(px_for_rand, 2, replace=FALSE)
#   while(pi[1]==pi[2]){
#     pi = sample(px_for_rand, 2, replace=FALSE)
#   }
#   pi = paste(pi[ order(pi) ], collapse=" + ")
#   ppairs = c(ppairs, pi)
# }
# as.data.frame(table(ppairs)) %>%
#   dplyr::arrange(desc(Freq))
# 
# 
# 
# 
# # # =============== other viz =============
# # 
# # graph2 = build_graph(pp %>% dplyr::filter(!param %in% c("Built-up land", "Healthcare travel time (log)")))
# # 
# # library(ggraph)
# # 
# # p2 = ggraph(graph2, layout="stress") + 
# #   geom_edge_link(aes(alpha=co_occ, width=co_occ, color=co_occ)) + 
# #   geom_node_point(aes(size=n_models, color=n_models)) + 
# #   geom_node_text(aes(label=abbrev2)) +
# #   scale_color_gradientn(colors=pals::ocean.deep(n=30)[9:26]) +
# #   scale_edge_colour_gradientn(colors=pals::ocean.deep(n=30)[9:26]) +
# #   maptheme +
# #   scale_edge_width(range = c(0.25, 3.5)) +
# #   scale_size(labels = c(1, 5, 10, 15, 20), breaks = c(1, 5, 10, 15, 20), range=c(3, 20)) +
# #   theme(legend.position="none") +
# #   scale_alpha(range=c(0.3, 1))
# # 
# # 
# # library(patchwork)
# # 
# # comb = p1 + p2
# # ggsave(comb, file="./output/plots/figure4_test_drivernetwork.jpg", device="jpg", units="in", width=12, height=7, dpi=600)
# # 
# # 
# # 
# # aa = pp %>%
# #   dplyr::group_by(Disease, param) %>%
# #   dplyr::summarise(sig = ifelse(sig == TRUE, 1, 0)) %>%
# #   tidyr::pivot_wider(names_from = c("param"), values_from = "sig")
# # aa = as.data.frame(apply(aa, 2, function(x) replace(x, is.na(x), 0)))
# # 
# # 
# # param_results = pp
# # 
# # 
# # build_graph = function(param_results){
# # 
# #   # significant effects only
# #   sig_efs = param_results %>%
# #     dplyr::filter(sig == TRUE)
# # 
# #   # 1. nodes
# #   nodes = sig_efs %>%
# #     dplyr::group_by(Disease) %>%
# #     dplyr::summarise(n_drivers = n_distinct(param)) %>%
# #     dplyr::arrange(desc(n_drivers)) %>%
# #     dplyr::rename(node=1)
# # 
# #   # add abbrevs for visualisation
# #   #nodes = left_join(nodes, read.csv("./scripts/04_results/driver_abbrevs.csv"), by=c("node"="driver"))
# # 
# #   # edges
# #   edges = sig_efs %>%
# #     dplyr::select(Disease, param) %>%
# #     dplyr::full_join(
# #       sig_efs %>% dplyr::select(Disease, param) %>% dplyr::rename(Disease2 = Disease)
# #     )
# #   edges = edges[ -which(edges$Disease == edges$Disease2), ]
# # 
# #   edges = edges %>%
# #     dplyr::group_by(Disease, Disease2) %>%
# #     dplyr::summarise(co_occ = n_distinct(param)) %>%
# #     dplyr::arrange(desc(co_occ)) %>%
# #     dplyr::left_join(
# #       nodes %>% dplyr::mutate(from = 1:nrow(nodes)) %>% dplyr::rename(Disease = node) %>% dplyr::select(Disease, from)
# #     ) %>%
# #     dplyr::left_join(
# #       nodes %>% dplyr::mutate(to = 1:nrow(nodes)) %>% dplyr::rename(Disease2 = node) %>% dplyr::select(Disease2, to)
# #     )
# # 
# #   # edges %>%
# #   #   ggplot() +
# #   #   geom_tile(aes(x = param, y=param2, fill=co_occ)) +
# #   #   theme_classic() +
# #   #   scale_fill_viridis_c(option="magma")
# # 
# #   # combine into network
# #   graph = tbl_graph(nodes=nodes, edges=edges, directed = FALSE, node_key="node")
# #   return(graph)
# # 
# # }
# # 
# # graph = build_graph(pp)
# # 
# # library(ggraph)
# # 
# # p1 = ggraph(graph, layout="circle") +
# #   geom_edge_link(aes(alpha=co_occ, width=co_occ, color=co_occ)) +
# #   geom_node_point(aes(size=n_drivers, color=n_drivers)) +
# #   geom_node_text(aes(label=node)) +
# #   scale_color_gradientn(colors=pals::ocean.deep(n=30)[9:26]) +
# #   scale_edge_colour_gradientn(colors=pals::ocean.deep(n=30)[9:26]) +
# #   maptheme +
# #   scale_edge_width(range = c(0.25, 3.5)) +
# #   scale_size(labels = c(1, 5, 10, 15, 20), breaks = c(1, 5, 10, 15, 20), range=c(3, 20)) +
# #   theme(legend.position="none") +
# #   scale_alpha(range=c(0.3, 1))
# # #
# # 
# # 
# # 
# # comb# ggraph(graph, layout="linear") + 
# # #   geom_edge_arc(aes(alpha=co_occ, width=co_occ, color=co_occ)) + 
# # #   geom_node_point(aes(size=n_models, color=n_models)) + 
# # #   maptheme +
# # #   scale_edge_width(range = c(0.1, 4)) +
# # #   scale_size(labels = c(1, 5, 10, 15, 20), breaks = c(1, 5, 10, 15, 20), range=c(2, 10))
# # 
# # 
# # library(tidygraph)
# # 
# # graph <- play_erdos_renyi(n = 10, p = 0.2) %>% 
# #   activate(nodes) %>% 
# #   mutate(class = sample(letters[1:4], n(), replace = TRUE)) %>% 
# #   activate(edges) %>% 
# # 
# # 
# # 
# # corr_mat = pp %>%
# #   dplyr::select(Disease, param, mean) %>%
# #   dplyr::full_join(
# #     pp %>% dplyr::select(Disease, param, mean) %>% dplyr::rename(param2 = param, mean2 = mean)
# #   )
# # 
# # 
# # pp %>%
# #   dplyr::group_by(param, param2)
# # 
# # 
# # 
# # a = corr_mat %>% dplyr::filter(param == "Forest cover" & param2 == "Vegetation heterogeneity")
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# dx = rbind(
#   c(1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1),
#   c(1, 0, 0, 1, 0, 0, 0, 0 ,0, 1, 1),
#   c(0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0),
#   c(1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
#   c(0, 0, 1, 0, 0 ,0, 0,0,0,0,0)
#   )
# 1-dist(dx, method="binary")
# 
# 
# 
# # ================= jaccard based dissimilarity ================
# 
# params = list.files("./output/model_outputs/disease_params/", pattern=".csv", full.names=TRUE)
# params = params[ -grep("infocriteria", params)]
# pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
# pp$mean = round(pp$mean, 3)
# pp$lower = round(pp$lower, 3)
# pp$upper = round(pp$upper, 3)
# pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0
# 
# # add abbrevs
# pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
#   dplyr::mutate(Disease = abbrev2)
# 
# pp = pp %>%
#   dplyr::filter(param != "Intercept") %>%
#   dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
#   dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
#                 param = replace(param, param == "livestock_log", "Livestock density (log)"),
#                 param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
#                 param = replace(param, param == "hunting", "Hunting pressure"),
#                 param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
#                 param = replace(param, param == "forest_cover", "Forest cover"),
#                 param = replace(param, param == "forest_loss", "Forest loss"),
#                 param = replace(param, param == "crop_expansion", "Cropland expansion"),
#                 param = replace(param, param == "evi_dissimilarity", "Vegetation heterogeneity"),
#                 param = replace(param, param == "evi_coefvar", "EVI variance"),
#                 param = replace(param, param == "precip_anomaly", "Precipitation anomaly"),
#                 param = replace(param, param == "urban_expansion", "Urban expansion"),
#                 param = replace(param, param == "urban_cover", "Built-up land"),
#                 param = replace(param, param == "crop_cover", "Cropland cover"),
#                 param = replace(param, param == "mining", "Mining"),
#                 param = replace(param, param == "tmean_change", "Temperature change"),
#                 param = replace(param, param == "precip_change", "Precipitation change"),
#                 param = replace(param, param == "social_vulnerability", "Social vulnerability"),
#                 param = replace(param, param == "biodiv_intact", "Biodiversity Intactness"),
#                 param = replace(param, param == "protected_areas", "Protected area coverage"),
#                 param = replace(param, param == "popdens_log", "Population density (log)")) 
# 
# # subset
# # pp = pp %>%
# #   dplyr::filter(human_infection_route %in% c("Vector", "Wildlife + vector"))
# 
# # order by n points
# np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
# pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)
# 
# # break down to causal broad
# pp = pp %>%
#   dplyr::filter(type == "Causal (broad)") %>%
#   dplyr::filter(!param %in% c("Built-up land", "Healthcare travel time (log)")) #%>%
#   #dplyr::filter(Num_observations > 70)
# 
# # create matrix of all-diseases all-drivers
# all_dat = expand.grid(unique(pp$Disease), unique(pp$param)) %>%
#   dplyr::rename("Disease"=1, "param"=2)
# sigs = pp %>%
#   dplyr::select(Disease, param, sig) %>%
#   distinct()
# all_dat = left_join(all_dat, sigs) %>%
#   dplyr::mutate(sig = ifelse(sig == TRUE, 1, 0),
#                 sig = replace(sig, is.na(sig), 0))
# all_dat = all_dat %>%
#   tidyr::pivot_wider(names_from = param, values_from = sig)
# 
# # calcuate jaccard distances
# jac = dist(
#   as.matrix(all_dat[ , 2:ncol(all_dat)]),
#   method = "binary"
# )
# 
# 
# # https://stats.stackexchange.com/questions/195446/choosing-the-right-linkage-method-for-hierarchical-clustering
# cor(jac,cophenetic(hclust(as.dist(jac), method="complete")))
# cor(jac,cophenetic(hclust(as.dist(jac), method="single")))
# cor(jac,cophenetic(hclust(as.dist(jac), method="ward.D")))
# cor(jac,cophenetic(hclust(as.dist(jac), method="average")))
# 
# # hierarchical clust
# hc = hclust(as.dist(jac), method="complete")
# plot(hc, labels = as.vector(all_dat$Disease))
# 
# rect.hclust(hc,
#             k = 25, # k is used to specify the number of clusters
#             border = "blue"
# )
# 
# # elbow point
# barplot(hc$height,
#         names.arg = (nrow(jac) - 1):1 # show the number of cluster below each bars
# )
# 
# ggdendrogram(hc, rotate = TRUE, size = 2)
# 
