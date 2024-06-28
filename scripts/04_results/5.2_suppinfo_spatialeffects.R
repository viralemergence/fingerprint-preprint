
# ====================== SUPP FIGURE: VISUALISE GEOSPATIAL EFFECTS =========================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")

# result location
loc = "./output/model_outputs/disease_models/"
ll = list.files(loc, full.names=TRUE)

# subset to 12 (randomly generated, specified here for reproducibility)
ll = ll[ grep("lassa|junin|rickett|jev|plague|melioi|oropouche|mayaro|hanta|sle|lyme|mpx", ll)]
#ll = ll[ sample(1:length(ll), 12, replace=FALSE)]
res = vector("list", length=12)

for(i in 1:length(ll)){
  
  load(file=ll[i])
  
  pts = model_obj$data
  ras = model_obj$spde; names(ras) = "field"
  shp = model_obj$study_area
  
  if(pts$Disease[1] == "Mayaro virus disease"){
    pts$Disease = "Mayaro fever"
  }
  if(pts$Disease[1] == "Monkeypox"){
    pts$Disease = "Mpox"
  }
  
  ras2 = ras %>% as.data.frame(xy=TRUE)
  lims = c(-max(abs(ras2$field), na.rm=TRUE), max(abs(ras2$field), na.rm=TRUE))
  
  px = ras %>%
    as.data.frame(xy=TRUE) %>%
    ggplot() + 
    geom_raster(aes(x, y, fill=field)) +
    #scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white") +
    scale_fill_gradientn(colors=rev(pals::brewer.brbg(200)), na.value="white", limits=lims) +
    theme_classic() +
    geom_sf(data=shp, fill=NA, col="black") + 
    geom_point(data=pts[ pts$presence == 1, ], aes(Longitude, Latitude), col="grey20", size=0.1, alpha=0.2) + 
    theme(legend.position = "bottom", plot.title = element_text(hjust=0.5, size=16),
          legend.title = element_blank(), axis.title = element_blank(), 
          axis.text = element_text(size=9)) +
    ggtitle(pts$Disease[1])
  
  res[[i]] = px
  
}

pcomb = gridExtra::grid.arrange(grobs = res, 
                        nrow=3)
ggsave(pcomb, file="./output/plots/SuppFigure_GMRFs.jpg", 
       device="jpg", units="in", width=15, height=11, dpi=600)
