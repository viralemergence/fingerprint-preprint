
# ====================== SUPP FIGURE: VISUALISE BACKGROUND POINTS AND GEOSPATIAL EFFECT =========================

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")

# result location
loc = "./output/model_outputs/disease_models/"
ll = list.files(loc, full.names=TRUE)

# load model
ll = ll[ grep("chagas", ll)]
load(file=ll)

pts = model_obj$data
shp = model_obj$study_area
field = model_obj$spde; names(field) = "spde"

# read population
pop = raster::raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_10km_agg.tif")
pop = crop(pop, shp)
names(pop) = "pop"

# mask
pop = raster::mask(x=pop, mask=shp, field=id)

# plot
p1 = pop %>%
  as.data.frame(xy=TRUE) %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=log(pop+1)), alpha=0.95) + 
  geom_sf(data=shp, fill=NA, col="grey10") + 
  scale_fill_gradientn(colors=rev(viridisLite::mako(200))[1:160], na.value="white", name="Population\n(log)") + 
  maptheme +
  geom_point(data=pts[ pts$presence==0, ], aes(Longitude, Latitude), color="black", size=0.6, alpha=0.6) + 
  geom_point(data=pts[ pts$presence==1, ], aes(Longitude, Latitude), color="red", size=0.9, alpha=0.9) + 
  theme(legend.position=c(0.2, 0.2),
        plot.title = element_text(hjust=0.5, size=18),
        legend.title = element_text(size=15)) + 
  ggtitle("Outbreak presence and background points")

field2 = field %>% as.data.frame(xy=TRUE)
lims = c(-max(abs(field2$spde), na.rm=TRUE), max(abs(field2$spde), na.rm=TRUE))

p2 = field2 %>%
  ggplot() + 
  geom_raster(aes(x, y, fill=spde), alpha=0.95) + 
  geom_sf(data=shp, fill=NA, col="grey10") + 
  #scale_fill_viridis_c(na.value="white", name="Outbreak\nrisk\n(log odds)") + 
  scale_fill_gradientn(colors=rev(pals::brewer.brbg(200)), na.value="white", limits=lims, name="Effect\n(log odds)") +
  maptheme +
  geom_point(data=pts[ pts$presence==1, ], aes(Longitude, Latitude), color="grey20", size=0.4, alpha=0.4) +
  theme(legend.position=c(0.2, 0.2),
        plot.title = element_text(hjust=0.5, size=18),
        legend.title = element_text(size=15)) + 
  ggtitle("Spatial random effect (Gauss-Markov random field)")

library(patchwork)
comb = p1 + p2  
#comb

ggsave(comb, 
       file="./output/plots/SuppFigure_BackgroundPoints_Chagas.jpg", 
       device="jpg", 
       units="in", 
       width=16, height=8, dpi=600, scale=0.85)












# subset to 12
set.seed(8204)
ll = ll[ sample(1:length(ll), 12, replace=FALSE)]
res = vector("list", length=12)

for(i in 1:length(ll)){
  
  load(file=ll[i])
  
  pts = model_obj$data
  ras = model_obj$spde; names(ras) = "field"
  shp = model_obj$study_area
  
  px = ras %>%
    as.data.frame(xy=TRUE) %>%
    ggplot() + 
    geom_raster(aes(x, y, fill=field)) +
    scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
    maptheme +
    geom_sf(data=shp, fill=NA, col="black") + 
    geom_point(data=pts[ pts$presence == 1, ], aes(Longitude, Latitude), col="red", size=0.25, alpha=0.5) + 
    theme(legend.position = "bottom", plot.title = element_text(hjust=0.5, size=16),
          legend.title = element_blank()) +
    ggtitle(pts$Disease[1])
  
  res[[i]] = px
  
}

pcomb = gridExtra::grid.arrange(grobs = res, 
                        nrow=3)
ggsave(pcomb, file="./output/plots/SuppFigure_SPDEs.jpg", 
       device="jpg", units="in", width=15, height=11, dpi=600)
