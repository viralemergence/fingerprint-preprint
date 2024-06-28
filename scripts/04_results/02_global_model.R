

# =============== Results from global modelling exercise (all diseases) ===================

# Models fitted in the 03_global_model.R script in the modelling directory

# setup objects
setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")




# ----------- 1. Global model including all diseases (n=31) -------------

# load three submodels
load(file="./output/model_outputs/global_models/global_m1_nonspatial.R")
load(file="./output/model_outputs/global_models/global_m2_spatial.R")
load(file="./output/model_outputs/global_models/global_m3_spatial.R")



# ----------- plot fixed effects ---------

# colours
cols = pals::brewer.brbg(200)[ c(1, 50, 170, 200)]
cols = pals::ocean.haline(200)[ c(35, 90, 140) ]

# effects
all_models = do.call(
  rbind.data.frame,
  list(m_nonspatial1$fx, m_spatial1$fx, m_spatial2$fx)
) %>%
  dplyr::mutate(
    model_type = ifelse(grepl("Nonspatial", type), "Non-spatial", "Spatial"),
    detection = ifelse(grepl("nondetection", type), "Non-detection", "Detection"),
    type2 = paste(model_type, detection, sep=" + "),
    type3 = "No adjustment",
    type3 = replace(type3, type2 == "Spatial + Non-detection", "Geospatial"),
    type3 = replace(type3, type2 == "Spatial + Detection", "Geospatial + Detection"),
    type2 = factor(type2, levels = rev(c("Non-spatial + Non-detection", "Non-spatial + Detection", "Spatial + Detection")), ordered=TRUE),
    type3 = factor(type3, levels = c("No adjustment", "Geospatial", "Geospatial + Detection"), ordered=TRUE)
  )

all_models$param = replace(all_models$param, all_models$param == "Travel time to health facility (log)", "Healthcare travel time (log)")
all_models$param = replace(all_models$param, all_models$param == "EVI dissimilarity", "Vegetation heterogeneity")
all_models$param = replace(all_models$param, all_models$param == "Biodiversity Intactness", "Biodiversity Intactness index")
all_models$param = replace(all_models$param, all_models$param == "Tmean change", "Temperature change")
all_models$param = replace(all_models$param, all_models$param == "Precip change", "Precipitation change")
all_models$param = replace(all_models$param, all_models$param == "Mining coverage", "Mining cover")
all_models$param = replace(all_models$param, all_models$param == "Protected area coverage", "Protected area cover")

# order by size of effect in full spatial model
param_order = all_models %>% dplyr::filter(type3 == "Geospatial + Detection" & param != "Intercept") %>% dplyr::arrange(desc(abs(mean)))

# plot
# p1 = all_models %>% 
#   dplyr::filter(paramname != "Intercept") %>%
#   dplyr::mutate(param = factor(param, levels=rev(param_order$param), ordered=TRUE)) %>%
#   ggplot() + 
#   geom_hline(yintercept=0, lty=2) +
#   geom_point(aes(param, mean, group=type3, col=type3), position=position_dodge(width=0.5), size=2) + 
#   geom_linerange(aes(param, ymin=lower, ymax=upper, group=type3, col=type3), position=position_dodge(width=0.5), size=0.4, alpha=0.7) +
#   theme_classic() + 
#   ylab("Scaled covariate effect on outbreak risk (log odds)") + 
#   xlab("Driver") + 
#   theme(legend.position = c(0.78, 0.12), 
#         plot.title = element_text(size=16, hjust=0.5),
#         axis.text = element_text(size=12), 
#         axis.title = element_text(size=13),
#         legend.title = element_blank(),
#         legend.text = element_text(size=10), 
#         #panel.grid.major.y = element_line(color="grey90"),
#         #legend.background = element_rect(fill=NA, color="grey70"),
#         legend.background = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(size=14)) +
#   scale_color_manual(values = cols[3:1], guide=guide_legend(reverse=TRUE)) +
#   coord_flip()

# p1 = all_models %>%
#   dplyr::filter(paramname != "Intercept") %>%
#   dplyr::mutate(param = factor(param, levels=rev(param_order$param), ordered=TRUE)) %>%
#   ggplot() +
#   geom_hline(yintercept=0, lty=2) +
#   geom_point(aes(param, mean, group=type3, col=type3), position=position_dodge(width=0.5), size=1.5) +
#   geom_linerange(aes(param, ymin=lower, ymax=upper, group=type3, col=type3), position=position_dodge(width=0.5), linewidth=0.4, alpha=1) +
#   theme_classic() +
#   ylab("Scaled covariate effect on outbreak risk (log odds)") +
#   xlab("Driver") +
#   scale_y_continuous(labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), limits=c(-1.6, 1.6)) +
#   theme(legend.position = "none",
#         plot.title = element_text(size=16, hjust=0.5),
#         axis.text.y = element_text(size=13),
#         axis.text.x = element_text(size=10),
#         axis.title = element_text(size=14),
#         legend.title = element_blank(),
#         legend.text = element_text(size=10),
#         panel.grid.major.y = element_line(color="grey90", linewidth=0.2),
#         panel.grid.major.x = element_line(color="grey95", linewidth=0.1),
#         #legend.background = element_rect(fill=NA, color="grey70"),
#         legend.background = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(size=14)) +
#   scale_color_manual(values = cols[3:1], guide=guide_legend(reverse=TRUE)) +
#   coord_flip() +
#   facet_wrap(~type3)

p1 = all_models %>%
  dplyr::filter(paramname != "Intercept") %>%
  dplyr::mutate(param = factor(param, levels=rev(param_order$param), ordered=TRUE)) %>%
  ggplot() +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type3), color="grey10", position=position_dodge(width=0.5), size=1.5) +
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type3), color="grey10", position=position_dodge(width=0.5), linewidth=0.4, alpha=1) +
  theme_classic() +
  ylab("Scaled covariate effect on outbreak risk (log odds)") +
  xlab("Driver") +
  scale_y_continuous(labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), limits=c(-1.6, 1.6)) +
  theme(legend.position = "none",
        plot.title = element_text(size=16, hjust=0.5),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        panel.grid.major.y = element_line(color="grey90", linewidth=0.2),
        panel.grid.major.x = element_line(color="grey95", linewidth=0.1),
        #legend.background = element_rect(fill=NA, color="grey70"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=14)) +
  scale_color_manual(values = cols[3:1], guide=guide_legend(reverse=TRUE)) +
  coord_flip() +
  facet_wrap(~type3)


# ----------- plot spatial field including outbreak locs ---------

ll = list.files("./output/model_df/", pattern=".csv", full.names=TRUE)
ll = ll[ -grep("brazil|usa|filo", ll)]
rr = read.csv(ll[1])
nn = names(rr)
nn = nn[ - which(nn %in% c("livestock_all", "tmean", "tmean_anomaly", "precip", "precip_anomaly", "evi_coefvar", "health_travel")) ]
dd = do.call(rbind.data.frame, lapply(ll, function(x) read.csv(x) %>% dplyr::select(all_of(nn)) %>% dplyr::mutate(file=x)))
dd$Disease = unlist(
  lapply(
    strsplit(
      unlist(lapply(strsplit(dd$file, "/"), "[", 4)), 
      "_"
    ),
    "[", 1
  )
)
dd = cbind(dd[ , (ncol(dd)-1) : ncol(dd), drop=FALSE], dd[ , 1:(ncol(dd)-2)])
dd = dd %>% 
  dplyr::select(-ADMcode, -record_id) %>%
  dplyr::filter(presence == 1) %>% 
  dplyr::filter(Longitude > -130) %>% 
  dplyr::filter(Latitude < 50)

ne = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
ne = ne[ ne$name != "Antarctica", ]
robinson = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
ne2 = st_transform(ne, robinson)

dd2 = dd %>% sf::st_as_sf(coords=c("Longitude", "Latitude")) %>%
  sf::st_set_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84") %>%
  sf::st_transform(robinson)
d1 = sf::st_coordinates(dd2)

spde_ras = m_spatial2$spde
crs(spde_ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
spde_ras = raster::projectRaster(spde_ras, crs=robinson)
spde_ras = spde_ras %>% as.data.frame(xy=TRUE) 

study_area = m_nonspatial1$shp
study_area = st_transform(study_area, robinson)
sf::sf_use_s2(FALSE)

lims = c(-max(abs(spde_ras$field), na.rm=TRUE), max(abs(spde_ras$field), na.rm=TRUE))

nec = rnaturalearth::ne_coastline(scale = 50, returnclass = c("sf"))
x = as.data.frame(sf::st_coordinates(sf::st_centroid(nec)))
nec = nec[ -which(x$Y < -60), ]
nec = st_transform(nec, robinson)

p2 = spde_ras %>%
  ggplot() + 
  geom_sf(data=ne2, fill="grey85", color=NA) +
  geom_raster(aes(x, y, fill=field)) + 
  geom_point(data=as.data.frame(d1), aes(X, Y), size=0.005, color="grey20", alpha=0.025) +
  scale_fill_gradientn(colors=rev(pals::brewer.brbg(200)), na.value="white", limits=lims, name="Spatial\nsampling\nbias") +
  maptheme +
  geom_sf(data=nec, fill=NA, col="grey50", size=0.05) + 
  theme(legend.position = c(0.16, 0.3), 
        plot.title = element_text(hjust=0.5, size=15), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.background = element_rect(fill="white", color="white"))

# combine
comb = gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(1, 1.2))

# add annotations
comb = ggpubr::as_ggplot(comb)  +
  cowplot::draw_plot_label(label = c("a", "b", "c", "d"), 
                           fontface = "plain", size = 28, 
                           x = c(0.25, 0.50, 0.745, 0.1), y = c(0.97, 0.97, 0.97, 0.52))

ggsave(comb, 
       file="./output/plots/Figure2_GlobalModel.jpg", 
       device="jpg", units="in", width=18, height=15, dpi=600, scale=0.6)







# ------------- 2. Global models of transmission subsets ----------------

load(file="./output/model_outputs/global_models/global_m3_spatial_zzvb.R")
m_zoo = mz_spatial2
load(file="./output/model_outputs/global_models/global_m3_spatial_vb.R")
m_vb = mz_spatial2
#load(file="./output/model_outputs/global_models/global_m3_spatial_vbo.R")

all_models = do.call(
  rbind.data.frame,
  list(
    m_zoo$fx %>% dplyr::mutate(type = "Zoonotic"), 
    m_vb$fx %>% dplyr::mutate(type = "Vector-borne")
  )
) %>%
  dplyr::filter(param != "Intercept")

all_models$param = replace(all_models$param, all_models$param == "Travel time to health facility (log)", "Healthcare travel time (log)")
all_models$param = replace(all_models$param, all_models$param == "EVI dissimilarity", "Vegetation heterogeneity")
all_models$param = replace(all_models$param, all_models$param == "Biodiversity Intactness", "Biodiversity Intactness index")
all_models$param = replace(all_models$param, all_models$param == "Tmean change", "Temperature change")
all_models$param = replace(all_models$param, all_models$param == "Precip change", "Precipitation change")
all_models$param = replace(all_models$param, all_models$param == "Mining coverage", "Mining cover")
all_models$param = replace(all_models$param, all_models$param == "Protected area coverage", "Protected area cover")

# order by size of parameter estimate in full global model
load(file="./output/model_outputs/global_models/global_m3_spatial.R")
px = m_spatial2$fx %>% 
  dplyr::filter(param != "Intercept") %>%
  dplyr::arrange(desc(abs(mean))) %>%
  dplyr::select(param)

px$param = replace(px$param, px$param == "Travel time to health facility (log)", "Healthcare travel time (log)")
px$param = replace(px$param, px$param == "EVI dissimilarity", "Vegetation heterogeneity")
px$param = replace(px$param, px$param == "Biodiversity Intactness", "Biodiversity Intactness index")
px$param = replace(px$param, px$param == "Tmean change", "Temperature change")
px$param = replace(px$param, px$param == "Precip change", "Precipitation change")
px$param = replace(px$param, px$param == "Mining coverage", "Mining cover")
px$param = replace(px$param, px$param == "Protected area coverage", "Protected area cover")

all_models$param = factor(all_models$param, levels=rev(px$param), ordered=TRUE)

#cols = pals::brewer.brbg(200)[ c(1, 50, 170, 200)]
cols = pals::ocean.haline(200)[ c(35, 140) ]

p1 = all_models %>%
  ggplot() +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(param, mean, group=type, col=type), position=position_dodge(width=0.5), size=1.5) +
  geom_linerange(aes(param, ymin=lower, ymax=upper, group=type, col=type), position=position_dodge(width=0.5), size=0.4, alpha=1) +
  theme_classic() +
  ylab("Scaled covariate effect on outbreak risk (log odds)") +
  xlab("Driver") +
  #scale_y_continuous(labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), limits=c(-1.6, 1.6)) +
  theme(legend.position = c(0.8, 0.15),
        plot.title = element_text(size=16, hjust=0.5),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=13),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        #panel.grid.major.y = element_line(color="grey90"),
        panel.border = element_rect(fill=NA, color="grey20"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=14)) +
  scale_color_manual(values = cols[1:2], guide=guide_legend(reverse=TRUE)) +
  coord_flip()

# spatial field 1

spde_ras = m_zoo$spde %>%
  as.data.frame(xy=TRUE) 

study_area = m_zoo$shp
sf::sf_use_s2(FALSE)

lims = c(-max(abs(spde_ras$field), na.rm=TRUE), max(abs(spde_ras$field), na.rm=TRUE))

p2 = spde_ras %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=field)) + 
  #geom_point(data=dd, aes(Longitude, Latitude), size=0.005, color="grey20", alpha=0.025) +
  scale_fill_gradientn(colors=rev(pals::brewer.brbg(200)), na.value="white", limits=lims, name="Spatial\nrandom\neffect") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black", size=0.4) + 
  theme(legend.position = "right", plot.title = element_text(hjust=0.5, size=13)) +
  ggtitle("Zoonotic (any route) [n=26]")

# spatial field 2

spde_ras = m_vb$spde %>%
  as.data.frame(xy=TRUE) 

study_area = m_vb$shp
sf::sf_use_s2(FALSE)

lims = c(-max(abs(spde_ras$field), na.rm=TRUE), max(abs(spde_ras$field), na.rm=TRUE))

p3 = spde_ras %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=field)) + 
  #geom_point(data=dd, aes(Longitude, Latitude), size=0.005, color="grey20", alpha=0.025) +
  scale_fill_gradientn(colors=rev(pals::brewer.brbg(200)), na.value="white", limits=lims, name="Spatial\nrandom\neffect") +
  maptheme +
  geom_sf(data=study_area, fill=NA, col="black", size=0.4) + 
  theme(legend.position = "right", plot.title = element_text(hjust=0.5, size=13)) +
  ggtitle("Vector-borne (any host) [n=20]")

pc = gridExtra::grid.arrange(p1,
                        gridExtra::grid.arrange(p2, p3, ncol=1),
                        ncol=2,
                        widths=c(1, 1.3))
ggsave(pc, file="./output/plots/ExtendedData5_GlobalModels_TransmissionTypes.jpg", 
       device="jpg", units="in", width=16, height=6, dpi=600)

