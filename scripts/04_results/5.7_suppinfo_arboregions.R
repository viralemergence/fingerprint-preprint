

# =============== Region-specific models for dengue and yellow fever ===================

# setup objects
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")


# ---------- files --------------

ff = list.files("./output/model_outputs/disease_pooled_arboregions/", pattern="csv", full.names=TRUE)
ff = ff[ grep("summary", ff)]
fx = do.call(rbind.data.frame, lapply(ff, read.csv))

# ---------- rename parameters ------------

fx = fx %>% 
  dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
              param = replace(param, param == "livestock_log", "Livestock density (log)"),
              param = replace(param, param == "tmean_anomaly", "Tmean anomaly"),
              param = replace(param, param == "forest_cover", "Forest cover"),
              param = replace(param, param == "forest_loss", "Forest loss"),
              param = replace(param, param == "crop_expansion", "Cropland expansion"),
              param = replace(param, param == "evi_dissimilarity", "Landscape heterogeneity"),
              param = replace(param, param == "precip_anomaly", "Precip anomaly"),
              param = replace(param, param == "urban_expansion", "Urban expansion"),
              param = replace(param, param == "urban_cover", "Urban cover"),
              param = replace(param, param == "hunting", "Hunting pressure"),
              param = replace(param, param == "crop_cover", "Cropland cover"),
              param = replace(param, param == "mining", "Mining cover"),
              param = replace(param, param == "tmean_change", "Temperature change"),
              param = replace(param, param == "precip_change", "Precipitation change"),
              param = replace(param, param == "social_vulnerability", "Social vulnerability"),
              param = replace(param, param == "biodiv_intact", "Biodiversity Intactness Index"),
              param = replace(param, param == "protected_areas", "Protected area cover"),
              param = replace(param, param == "popdens_log", "Population density (log)")) %>%
  dplyr::filter(param != "Intercept")

#
fac_ord = fx %>%
  dplyr::select(param, median) %>%
  dplyr::group_by(param) %>%
  dplyr::summarise(max = max(median))

# variable type
fac_ord$vartype = "Land use intensity"
fac_ord$vartype[ fac_ord$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
fac_ord$vartype[ fac_ord$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
fac_ord$vartype[ fac_ord$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
fac_ord$vartype[ fac_ord$param %in% c("Forest cover", "Cropland cover", "Landscape heterogeneity", "Biodiversity Intactness Index") ] = "Ecosystem structure"
fac_ord$vartype = factor(fac_ord$vartype, levels=c("Detection", "Ecosystem structure", "Socioeconomic",  "Climate change", "Land use intensity"))

# order by max effect and vartype
fac_ord = fac_ord %>%
  dplyr::arrange(vartype, desc(max))

fx = fx %>%
  dplyr::mutate(param = factor(param, levels = rev(fac_ord$param), ordered=TRUE)) %>%
  dplyr::mutate(param_cd = as.numeric(param))

fac_ord = fac_ord %>%
  dplyr::mutate(param = factor(param, levels = rev(fac_ord$param), ordered=TRUE)) %>%
  dplyr::mutate(param_cd = as.numeric(param))

fx = fx %>%
  dplyr::left_join(
    data.frame(
      region = c("AFR", "ASIA", "LAC"),
      region2 = c("Africa", "Asia & Pacific", "Latin America & Caribbean")
    )
  )
fx$region2 = factor(fx$region2, levels=rev(c("Latin America & Caribbean", "Africa", "Asia & Pacific")), ordered=TRUE)

p1 = ggplot() + 
  geom_rect(data = fac_ord, aes(xmin=param_cd - 0.5, xmax=param_cd + 0.5, ymin=-Inf, ymax=Inf, fill=vartype), alpha=0.4) +
  theme_classic() +
  geom_vline(xintercept=seq(0.5, nrow(fac_ord)+0.5), linewidth=0.3, color="white") +
  geom_hline(yintercept=0, lty=2, color="grey40", linewidth=0.3) +
  geom_linerange(data = fx, aes(param_cd, ymin=lower, ymax=upper, group=region2), position=position_dodge(width=0.9), color="grey20", size=0.3, alpha=1) +
  geom_point(data = fx, aes(param_cd, median, group=region2, pch=region2), position=position_dodge(width=0.9), color="grey15", size=2) +
  theme(
    legend.position="right",
    plot.title = element_text(size=16, hjust=0.5),
    axis.text.y = element_text(size=13),
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=13),
    legend.title = element_text(size=13),
    legend.text = element_text(size=11.2),
    #panel.grid.major.x = element_line(color="grey70", linewidth=0.5),
    panel.border = element_rect(fill=NA, color=NA),
    axis.line.y = element_blank(),
    axis.line.x = element_line(color="grey50"),
    axis.ticks.y = element_line(color="grey50"),
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=14.5)) +
  #scale_color_manual(values = cols[1:2], guide=guide_legend(reverse=TRUE)) +
  scale_x_continuous(breaks=fac_ord$param_cd, labels=fac_ord$param) +
  #scale_y_continuous(breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
  ylab("Z-score scaled effect on outbreak event risk (log odds)") + 
  xlab("Driver") +
  coord_flip() +
  facet_wrap(~Disease, nrow=1, scales="free_x") +
  scale_shape_manual(values=c(19, 17, 15), name="Region") +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6) ),
         pch = guide_legend(reverse=TRUE))
  
# # save
# ggsave(p1, 
#        file="./output/plots/ExtendedData_ArbovirusesByRegion.jpg", 
#        device="jpg", units="in", width=17, height=9, dpi=600, scale=0.6)




# -------------- create map to accompany ----------------

# yellow fever
dz1 = read.csv("./output/spillovers_processed/spillovers_yfv.csv") %>%
  distinct() %>%
  dplyr::mutate(record_id = 1:length(Disease)) %>%
  dplyr::filter(!is.na(Longitude)) %>%
  dplyr::filter(Year >= 1985)

# dengue
dz2 = read.csv("./output/spillovers_processed/spillovers_dengue.csv") %>%
  distinct() %>%
  dplyr::mutate(record_id = 1:length(Disease)) %>%
  dplyr::filter(!is.na(Longitude)) %>%
  dplyr::filter(Year >= 1985)

# combine
dz = rbind(
  dz1 %>% dplyr::select(Disease, Longitude, Latitude),
  dz2 %>% dplyr::select(Disease, Longitude, Latitude)
)

# create study area polygon with country/coastal boundaries
study_area = createStudyAreaPolygon(lat=dz$Latitude, lon=dz$Longitude, buffer=FALSE, extent_border = 8)

# adding a map
# ne = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
# ne = ne[ ne$name != "Antarctica", ]
robinson = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#robinson = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

ne2 = st_transform(study_area, robinson)
dz2 = dz %>% dplyr::filter(!is.na(Longitude))
coordinates(dz2) = ~Longitude + Latitude
dz2 = st_as_sf(dz2)
st_crs(dz2) = st_crs(ne)
dz2 = st_transform(dz2, robinson)

dz2$Disease = factor(dz2$Disease, levels=c("Yellow fever", "Dengue"), ordered=TRUE)

spills = ggplot() + 
  geom_sf(data=ne2, fill="grey94", color="grey75", size=0.2) + 
  maptheme + 
  geom_sf(data=dz2, aes(col=Disease), size=0.45, alpha=0.4) +
  guides(colour = guide_legend(override.aes = list(size=5, alpha=1), ncol=2, reverse=TRUE)) + 
  theme(legend.position="bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size=13)) +
  scale_color_manual(values = pals::brewer.brbg(n=9)[c(2, 7)]) 


pc = gridExtra::grid.arrange(p1, spills, ncol=1, heights=c(1, 0.6))
ggsave(pc, 
       file="./output/plots/SuppFigure_ArbovirusesByRegion.jpg", 
       device="jpg", units="in", width=17, height=14.5, dpi=600, scale=0.6)
