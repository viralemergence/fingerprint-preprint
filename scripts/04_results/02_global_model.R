

# =============== Results from global modelling exercise (all diseases) ===================

# Script reads overall global model, plus subsets by transmission type and region
# Combines posterior samples across all submodels then generates overall results figures

# setup objects
PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")




# ----------- 1. Global model including all diseases (n=31) -------------

# get fixed effects estimates across bootstrap ensemble for each model type

getBSparams = function(loc, model_type){
  
  f1 = list.files(loc, pattern=model_type, full.names=TRUE)
  f1 = f1[ grep(".csv", f1) ]
  fx = do.call(rbind.data.frame, lapply(f1, read.csv))
  fxd = data.frame(
    param = names(fx[ , !names(fx) %in% c("type", "id")]),
    lower95 = apply(fx[ , !names(fx) %in% c("type", "id")], 2, quantile, c(0.025)),
    lower67 = apply(fx[ , !names(fx) %in% c("type", "id")], 2, quantile, c(0.1666)),
    median = apply(fx[ , !names(fx) %in% c("type", "id")], 2, quantile, c(0.5)),
    upper67 = apply(fx[ , !names(fx) %in% c("type", "id")], 2, quantile, c(0.8333)),
    upper95 = apply(fx[ , !names(fx) %in% c("type", "id")], 2, quantile, c(0.975))
  ) %>%
    dplyr::mutate(param = sapply(strsplit(param, "[.]"), "[", 1),
                  type = model_type)
  row.names(fxd) = c()
  return(fxd)
  
}

locx = "./output/model_outputs/global_model_bs/outputs/all_diseases/"

fx = do.call(
  rbind.data.frame,
  list(getBSparams(locx, "Nonspatial_nondetection"), 
       getBSparams(locx, "Spatial_nondetection"), 
       getBSparams(locx, "Spatial_detection")
  )
)

# get SPDE rasters
f = list.files("./output/model_outputs/global_model_bs/outputs/all_diseases/", pattern=".R", full.names=TRUE)
f = f[ grep("Spatial_detection", f) ]
ss = raster::stack()
for(i in 1:length(f)){
  load(f[i])
  ss = raster::stack(ss, spde)
}


# ----------- plot fixed effects ---------

# colours
cols = pals::brewer.brbg(200)[ c(1, 50, 170, 200)]
cols = pals::ocean.haline(200)[ c(35, 90, 140) ]

# effects
all_models = fx %>%
  dplyr::mutate(
    model_type = ifelse(grepl("Nonspatial", type), "Non-spatial", "Spatial"),
    detection = ifelse(grepl("nondetection", type), "Non-detection", "Detection"),
    type2 = paste(model_type, detection, sep=" + "),
    type3 = "No adjustment",
    type3 = replace(type3, type2 == "Spatial + Non-detection", "Geospatial"),
    type3 = replace(type3, type2 == "Spatial + Detection", "Geospatial + Detection"),
    type2 = factor(type2, levels = rev(c("Non-spatial + Non-detection", "Non-spatial + Detection", "Spatial + Detection")), ordered=TRUE),
    type3 = factor(type3, levels = c("No adjustment", "Geospatial", "Geospatial + Detection"), ordered=TRUE)
  ) %>%
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
                param = replace(param, param == "popdens_log", "Population density (log)"))

# add detection covars to df for non detecton models
dcs = all_models %>% dplyr::filter(param %in% c("Healthcare travel time (log)", "Urban cover"))
dcs = rbind(
  dcs %>% dplyr::mutate(type3 = "No adjustment"),
  dcs %>% dplyr::mutate(type3 = "Geospatial")
) %>% 
  dplyr::mutate(
  lower95 = NA, lower67=NA, median=NA, upper67=NA, upper95=NA
)
all_models = all_models %>%
  rbind(dcs)

# variable type
all_models$vartype = "Land use intensity"
all_models$vartype[ all_models$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
all_models$vartype[ all_models$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
all_models$vartype[ all_models$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
all_models$vartype[ all_models$param %in% c("Forest cover", "Cropland cover", "Landscape heterogeneity", "Biodiversity Intactness Index") ] = "Ecosystem structure"
all_models$vartype = factor(all_models$vartype, levels=c("Detection", "Ecosystem structure", "Socioeconomic",  "Climate change", "Land use intensity"))

# order by vartype and then by size of effect in full spatial model
param_order = all_models %>% 
  dplyr::filter(type3 == "Geospatial + Detection" & param != "Intercept") %>% 
  dplyr::arrange(vartype, desc(median))
all_models = all_models %>% 
  dplyr::mutate(param = factor(param, levels=rev(param_order$param), ordered=TRUE)) %>%
  dplyr::mutate(param_cd = as.numeric(param)) # create a numeric version for plotting background rectangles

# order for labels
ords = all_models %>%
  dplyr::filter(type3 == "Geospatial + Detection") %>%
  dplyr::arrange(param_cd)

p1 = all_models %>%
  dplyr::mutate(
    type3 = factor(type3, levels = c("No adjustment", "Geospatial", "Geospatial + Detection"), ordered=TRUE)
  ) %>%  
  ggplot() +
  geom_rect(aes(xmin=param_cd - 0.5, xmax=param_cd + 0.5, ymin=-Inf, ymax=Inf, fill=vartype), alpha=0.4) +
  geom_vline(xintercept=seq(0.5, nrow(ords)+0.5), linewidth=0.3, color="white") +
  geom_hline(yintercept=0, lty=2, color="grey40", linewidth=0.3) +
  geom_linerange(aes(param_cd, ymin=lower95, ymax=upper95), position=position_dodge(width=0.5), color="grey15", size=0.4, alpha=1) +
  geom_linerange(aes(param_cd, ymin=lower67, ymax=upper67), position=position_dodge(width=0.5), color="grey15", size=1.3, alpha=1) +
  geom_point(aes(param_cd, median), position=position_dodge(width=0.5), color="grey10", size=1.8) +
  theme_classic() +
  theme(
    legend.position="bottom",
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
  scale_x_continuous(breaks=ords$param_cd, labels=ords$param) +
  #scale_y_continuous(breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
  ylab("Z-score scaled effect on outbreak event risk (log odds)") + 
  xlab("Driver") +
  coord_flip() +
  facet_wrap(~type3, nrow=1) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6) ) )


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

study_area = createStudyAreaPolygon(lat=dd$Latitude, lon=dd$Longitude, buffer=FALSE, extent_border = 1)

ne = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
ne = ne[ ne$name != "Antarctica", ]
robinson = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
ne2 = st_transform(ne, robinson)

dd2 = dd %>% sf::st_as_sf(coords=c("Longitude", "Latitude")) %>%
  sf::st_set_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84") %>%
  sf::st_transform(robinson)
d1 = sf::st_coordinates(dd2)

spde_ras = mean(ss)
crs(spde_ras) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
spde_ras = raster::projectRaster(spde_ras, crs=robinson)
spde_ras = spde_ras %>% as.data.frame(xy=TRUE) 
names(spde_ras)[3] = "field"

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
  scale_fill_gradientn(colors=rev(pals::brewer.brbg(200)), na.value="white", limits=lims, name="Residual\n(unexplained)\noutbreak risk") +
  maptheme +
  geom_sf(data=nec, fill=NA, col="grey50", size=0.05) + 
  theme(legend.position = c(0.16, 0.3), 
        plot.title = element_text(hjust=0.5, size=15), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.background = element_rect(fill="white", color="white"))

# combine
comb = gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(1.1, 1))

# add annotations
comb = ggpubr::as_ggplot(comb)  +
  cowplot::draw_plot_label(label = c("a", "b", "c", "d"), 
                           fontface = "plain", size = 28, 
                           x = c(0.24, 0.485, 0.735, 0.1), y = c(1, 1, 1, 0.505))

# save
ggsave(comb, 
       file="./output/plots/Figure2_GlobalModel_BS_CS.jpg", 
       device="jpg", units="in", width=18.5, height=15, dpi=600, scale=0.6)





# ------------- 2. Global models of transmission subsets ----------------

# 1. zoonotic (any route)

locx = "./output/model_outputs/global_model_bs/outputs/zoonotic"
fx_zoo = getBSparams(locx, "Spatial_detection")

# 2. vector-borne (any host)

locx = "./output/model_outputs/global_model_bs/outputs/vectorborne/"
fx_vec = getBSparams(locx, "Spatial_detection")

# combine and plot effects
models = fx_vec %>% 
  dplyr::mutate(model = "Vector-borne") %>%
  rbind(
    fx_zoo %>% dplyr::mutate(model = "Zoonotic")
  ) %>%
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
                param = replace(param, param == "popdens_log", "Population density (log)"))

models = models %>% 
  dplyr::mutate(param = factor(param, levels=rev(param_order$param), ordered=TRUE)) %>%
  dplyr::mutate(param_cd = as.numeric(param)) # create a numeric version for plotting background rectangles

# variable type
models$vartype = "Land use intensity"
models$vartype[ models$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
models$vartype[ models$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
models$vartype[ models$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
models$vartype[ models$param %in% c("Forest cover", "Cropland cover", "Landscape heterogeneity", "Biodiversity Intactness Index") ] = "Ecosystem structure"
models$vartype = factor(models$vartype, levels=c("Detection", "Ecosystem structure", "Socioeconomic",  "Climate change", "Land use intensity"))

# order for labels
ords = models %>%
  dplyr::filter(model == "Vector-borne") %>%
  dplyr::arrange(param_cd)

# add names
models$route = ifelse(models$model == "Vector-borne", "Vector-borne (any host) [n=20]", "Zoonotic (any route) [n=26]")

p_trans = models %>%
  dplyr::mutate(
    route = factor(route, levels = c("Zoonotic (any route) [n=26]", "Vector-borne (any host) [n=20]"), ordered=TRUE)
  ) %>%  
  ggplot() +
  geom_rect(aes(xmin=param_cd - 0.5, xmax=param_cd + 0.5, ymin=-Inf, ymax=Inf, fill=vartype), alpha=0.25) +
  geom_vline(xintercept=seq(0.5, nrow(ords)+0.5), linewidth=0.3, color="white") +
  geom_hline(yintercept=0, lty=2, color="grey40", linewidth=0.3) +
  geom_linerange(aes(param_cd, ymin=lower95, ymax=upper95, group=route), position=position_dodge(width=0.75), color="grey25", size=0.4, alpha=1) +
  geom_linerange(aes(param_cd, ymin=lower67, ymax=upper67, group=route), position=position_dodge(width=0.75), color="grey25", size=1, alpha=1) +
  geom_point(aes(param_cd, median, pch=route), position=position_dodge(width=0.75), color="grey10", size=2, fill="white") +
  theme_classic() +
  theme(
    legend.position="right",
    plot.title = element_text(size=15, hjust=0.5, face="bold"),
    axis.text.y = element_text(size=11),
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=11),
    legend.title = element_text(size=12),
    legend.text = element_text(size=11),
    legend.box.background = element_rect(color="grey20"),
    #panel.grid.major.x = element_line(color="grey70", linewidth=0.5),
    panel.border = element_rect(fill=NA, color=NA),
    axis.line.y = element_blank(),
    axis.line.x = element_line(color="grey50"),
    axis.ticks.y = element_line(color="grey50"),
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=12)) +
  #scale_color_manual(values = cols[1:2], guide=guide_legend(reverse=TRUE)) +
  scale_x_continuous(breaks=ords$param_cd, labels=ords$param) +
  scale_shape_manual(values=c(15, 17), name="Transmission") +
  #scale_y_continuous(breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
  ylab("Z-score scaled effect on outbreak event risk (log odds)") + 
  xlab("Driver") +
  coord_flip() +
  guides(shape = guide_legend(reverse=TRUE)) +
  ggtitle("Global model by pathogen transmission route") + 
  guides(fill = guide_legend(override.aes = list(alpha = 0.6) ) )



# ------------ Region-specific models including all diseases per region ---------------

locx = "output/model_outputs/global_model_bs/outputs/regions/"

getBSparams = function(loc, model_type){
  
  f1 = list.files(loc, pattern=model_type, full.names=TRUE)
  f1 = f1[ grep(".csv", f1) ]
  fx = do.call(rbind.data.frame, lapply(f1, read.csv))
  fxd = data.frame(
    param = names(fx[ , !names(fx) %in% c("region", "id")]),
    lower95 = apply(fx[ , !names(fx) %in% c("region", "id")], 2, quantile, c(0.025)),
    lower67 = apply(fx[ , !names(fx) %in% c("region", "id")], 2, quantile, c(0.1666)),
    median = apply(fx[ , !names(fx) %in% c("region", "id")], 2, quantile, c(0.5)),
    upper67 = apply(fx[ , !names(fx) %in% c("region", "id")], 2, quantile, c(0.8333)),
    upper95 = apply(fx[ , !names(fx) %in% c("region", "id")], 2, quantile, c(0.975))
  ) %>%
    dplyr::mutate(param = sapply(strsplit(param, "[.]"), "[", 1),
                  region = fx$region[1])
  row.names(fxd) = c()
  return(fxd)
  
}

fx_rg = do.call(
  rbind.data.frame,
  list(
    getBSparams(locx, model_type = "EAP"),
    getBSparams(locx, model_type = "SSA"),
    getBSparams(locx, model_type = "NAM"),
    getBSparams(locx, model_type = "LAC"),
    getBSparams(locx, model_type = "SAS")
  )
)

# add "NA" for social vul for missing regions
sv = fx_rg %>% dplyr::filter(param == "social_vulnerability" & region == "Sub-Saharan Africa")
sv = do.call(
  rbind.data.frame,
  list(sv, sv, sv)
) %>%
  dplyr::mutate(
    region = c("North America", "East Asia & Pacific", "South Asia"),
    lower95 = NA, lower67=NA, median=NA, upper67=NA, upper95=NA
  )

ht= fx_rg %>% dplyr::filter(param == "health_travel_log" & region == "Sub-Saharan Africa") %>%
  dplyr::mutate(
    region = c("North America"),
    lower95 = NA, lower67=NA, median=NA, upper67=NA, upper95=NA
  )


fp = fx_rg %>%
  rbind(sv) %>%
  rbind(ht) %>%
  #dplyr::filter(region == "Latin America & Caribbean") %>%
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
  dplyr::mutate(param = factor(param, levels=rev(param_order$param), ordered=TRUE)) %>%
  dplyr::mutate(param_cd = as.numeric(param))

fp$vartype = "Land use change/intensity"
fp$vartype[ fp$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
fp$vartype[ fp$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
fp$vartype[ fp$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
fp$vartype[ fp$param %in% c("Forest cover", "Cropland cover", "Landscape heterogeneity", "Biodiversity Intactness Index") ] = "Ecosystem structure"
fp$vartype = factor(fp$vartype, levels=c("Detection", "Ecosystem structure", "Socioeconomic",  "Climate change", "Land use change/intensity"))

ords = fp %>%
  dplyr::filter(region == "Latin America & Caribbean") %>%
  dplyr::arrange(param_cd)

fp$region = factor(fp$region,
                   levels=c("North America", "South Asia", "East Asia & Pacific", "Latin America & Caribbean", "Sub-Saharan Africa"))

p_reg = fp %>%
  ggplot() +
  geom_rect(aes(xmin=param_cd - 0.5, xmax=param_cd + 0.5, ymin=-Inf, ymax=Inf, fill=vartype), alpha=0.4) +
  geom_vline(xintercept=seq(0.5, nrow(ords)+0.5), linewidth=0.3, color="white") +
  geom_hline(yintercept=0, lty=2, color="grey40", linewidth=0.3) +
  geom_linerange(aes(param_cd, ymin=lower95, ymax=upper95), position=position_dodge(width=0.5), color="grey15", size=0.6, alpha=1) +
  geom_linerange(aes(param_cd, ymin=lower67, ymax=upper67), position=position_dodge(width=0.5), color="grey15", size=1.2, alpha=1) +
  geom_point(aes(param_cd, median), position=position_dodge(width=0.5), color="grey10", size=1.7) +
  theme_classic() +
  theme(#legend.position.inside = c(0.8, 0.15),
    legend.position="none",
    plot.title = element_text(size=15, hjust=0.5, face="bold"),
    axis.text.y = element_text(size=11),
    axis.text.x = element_text(size=10),
    axis.title = element_text(size=11),
    legend.title = element_blank(),
    legend.text = element_text(size=11),
    #panel.grid.major.y = element_line(color="grey90"),
    panel.border = element_rect(fill=NA, color=NA),
    axis.line.y = element_blank(),
    axis.line.x = element_line(color="grey50"),
    axis.ticks.y = element_line(color="grey50"),
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=12)) +
  #scale_color_manual(values = cols[1:2], guide=guide_legend(reverse=TRUE)) +
  scale_x_continuous(breaks=ords$param_cd, labels=ords$param) +
  #scale_y_continuous(breaks=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5), labels=c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)) +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
  ylab("Z-score scaled effect on outbreak event risk (log odds)") + 
  xlab("Driver") +
  coord_flip() +
  facet_wrap(~region, nrow=2) +
  ggtitle("Region-specific models")


# ---------- combine transmission route and region into a single ED plot -----------

# combine
comb = gridExtra::grid.arrange(p_trans, p_reg, nrow=2, heights=c(1, 1.55))

# add annotations
comb = ggpubr::as_ggplot(comb)  +
  cowplot::draw_plot_label(label = c("a", "b"), 
                           fontface = "plain", size = 28, 
                           x = c(0.05, 0.05), y = c(1, 0.6))

# save
ggsave(comb, 
       file="./output/plots/SuppFigure_GlobalModel_Subsets.jpg", 
       device="jpg", units="in", width=15, height=18, dpi=600, scale=0.6)

  
