
# ================== Map and description of entire fingerprint database ==================

setwd("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")

library(patchwork); library(ggplot2)


# combine spillover data

read_spillovers = function(x){
  cc = read.csv(x) %>%
    dplyr::mutate(IfPolygon_AdminLevel = as.character(IfPolygon_AdminLevel))
  if("ADM2code" %in% names(cc)){ 
    cc$ADM2code = as.character(cc$ADM2code)
  }
  if("NumCases" %in% names(cc)){ 
    cc$NumCases = as.integer(cc$NumCases)
  }
  if("ADM1code" %in% names(cc)){ 
    cc$ADM1code = as.character(cc$ADM1code)
  }
  if("ADMcode" %in% names(cc)){ 
    cc$ADMcode = as.character(cc$ADMcode)
  }
  if("Year" %in% names(cc)){ 
    cc$Year = as.integer(cc$Year)
  }
  cc
}

sp_files = list.files("./output/spillovers_processed/", pattern=".csv", full.names=TRUE)
sp = do.call(
  dplyr::bind_rows,
  lapply(sp_files, read_spillovers)
)

# ones that do not enter the analyses
sp = sp[ -which(sp$Disease == "Chikungunya" & sp$Source == "Brazil DATASUS system"), ]
sp = sp[ -which(sp$Disease == "Crimean-Congo haemorrhagic fever" & sp$Source == "Ak2020"), ]
sp$Disease[ sp$Disease == "Mayaro fever" ] = "Mayaro virus disease"

# arbonet
arbo = read.csv("C:/Users/roryj/Documents/Research/projects_current/202011_fingerprint/arbonet/fingerprint_formatted/spillovers_arbonet.csv") %>%
  dplyr::filter(Presentation == "All cases") %>%
  dplyr::filter(!Disease %in% c("Zika", "Chikungunya"))
arbo$Longitude = arbo$Longitude + sample(seq(-0.01, 0.01, length.out=500), nrow(arbo), replace=TRUE)
arbo$Latitude = arbo$Latitude + sample(seq(-0.01, 0.01, length.out=500), nrow(arbo), replace=TRUE)

# combined
sp = sp %>%
  dplyr::select(Disease, Year, Longitude, Latitude, NumCases, Source) %>%
  rbind(
    arbo %>% 
      dplyr::select(Disease, Year, Longitude, Latitude, NumCases, Source) 
  )

sp = sp %>%
  left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) 

sp$abbrev = sp$abbrev6


# ========== description ============

desc = sp %>% 
  dplyr::group_by(Disease) %>%
  dplyr::summarise(
    Years = paste(c(min(Year, na.rm=TRUE), max(Year, na.rm=TRUE)), collapse="-"),
    Num_outbreaks = length(NumCases)
  )


# =========== visualise ============

# records per disease
p1 = as.data.frame(table(sp$abbrev )) %>%
  dplyr::rename(abbrev=1) %>%
  dplyr::arrange(desc(Freq)) %>%
  dplyr::left_join(
    sp %>% dplyr::select(abbrev,  dz_class) %>% distinct()
  ) %>%
  dplyr::mutate(logFreq = log(Freq), 
                abbrev = factor(abbrev, levels=rev(abbrev), ordered=TRUE)) %>%
  ggplot() + 
  geom_segment(aes(x = 0, xend=Freq, y=abbrev, yend=abbrev, color= dz_class), show.legend = FALSE) + 
  geom_point(aes(Freq, abbrev, fill=dz_class), size=3, pch=23, color="grey20", stroke=0.35) +
  theme_classic() + 
  scale_x_log10() +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(4), name="Transmission") + 
  scale_color_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(4), name="Transmission") + 
  xlab("Total outbreak events") + ylab("Disease") +
  theme(legend.position = "none",
        axis.title = element_text(size=15),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))

# records over time
time = as.data.frame(table(sp$Year, sp$dz_class)) %>%
  dplyr::rename(Year=1, dz_class = 2) %>%
  dplyr::mutate(Year = as.integer(as.vector(Year)))

time = left_join(
  data.frame(Year=rep(1950:2020, n_distinct(sp$dz_class)),
             dz_class = rep(unique(sp$dz_class), each=length(1950:2020))),
  time
) %>%
  dplyr::mutate(Freq = replace(Freq, is.na(Freq), 0),
                logFreq = log(Freq+1)) %>%
  dplyr::filter(!is.na(dz_class))

p2 = time %>%
  ggplot() + 
  geom_area(aes(Year, Freq, fill=dz_class)) +
  geom_area(aes(Year, Freq, group=dz_class), color="grey20", fill=NA, size=0.15, alpha=0.6) +
  #geom_vline(xintercept=1985, lty=2, size=0.6, alpha=0.7, color="grey20") +
  scale_x_continuous(labels=seq(1950, 2020, by=10), breaks=seq(1950, 2020, by=10)) +
  scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(4), name="Transmission") + 
  theme_classic() + 
  theme(legend.position=c(0.3, 0.8),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14)) +
  ylab("Number of outbreak events") + xlab("Year") 


# combine and save
comb = p2 + p1 
comb
ggsave(comb, file="./output/plots/Figure1_histo.jpg",
       device="jpg",
       units="in",
       dpi=900,
       width=12, height=5.5, scale=0.9)


# --------- map ------------

# adding a map
ne = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
ne = ne[ ne$name != "Antarctica", ]
robinson = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#robinson = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

ne2 = st_transform(ne, robinson)
sp2 = sp %>% dplyr::filter(!is.na(Longitude))
coordinates(sp2) = ~Longitude + Latitude
sp2 = st_as_sf(sp2)
st_crs(sp2) = st_crs(ne)
sp2 = st_transform(sp2, robinson)

spills = ggplot() + 
  geom_sf(data=ne2, fill="grey94", color="grey75", size=0.2) + 
  maptheme + 
  geom_sf(data=sp2, aes(col=abbrev), size=0.065, alpha=0.4) + 
  guides(colour = guide_legend(override.aes = list(size=5), ncol=2, alpha=1)) + 
  theme(legend.position="right", legend.title = element_blank(), 
        legend.text = element_text(size=11)) +
  scale_color_manual(values = rev(pals::tol.rainbow(n=n_distinct(sp2$Disease))))

#  save
ggsave(spills, file="./output/plots/Figure1_map.jpg", device="jpg", units="in", dpi=900, width=12, height=5)

# map and time series figure are combined separately


spills2 = ggplot() + 
  geom_sf(data=ne2, fill="grey94", color="grey75", size=0.2) + 
  maptheme + 
  geom_sf(data=sp2, aes(col=abbrev), size=0.065, alpha=0.4) + 
  guides(colour = guide_legend(override.aes = list(size=5), ncol=8, alpha=1)) + 
  theme(legend.position=c(0.55, -0.12), legend.title = element_blank(), 
        legend.text = element_text(size=10.5)) +
  scale_color_manual(values = rev(pals::tol.rainbow(n=n_distinct(sp2$Disease))))
ggsave(spills2, file="./output/plots/Figure1_map_vert.jpg", device="jpg", units="in", dpi=900, width=12, height=8)


