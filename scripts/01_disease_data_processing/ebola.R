
# ======================= Ebola ===========================

# read in dataset from Pigott et al, and supplement of spillover index cases since 2015
# only at polygon level for many spillovers; combine polygons from pigott with health zone polygons from DRC
# extract lat-lon centroids for spatial locs, and save polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")



# -------- 1. Pigott et al 2016 eLife ----------------

ebo1 = read.csv("./data/spillovers/Ebola/pigott2016/elife_advances_human_index.csv")

# points
ebo1p = ebo1 %>%
  dplyr::filter(shape == "point") %>%
  dplyr::select(Year.Start, Virus, Country, Apparent.Origin, Long, Lat, shape, Outbreak_ID) %>%
  dplyr::rename("Longitude"=Long, "Latitude"=Lat)

# polyons
poly = ebo1 %>%
  dplyr::filter(shape == "polygon") %>%
  dplyr::select(Year.Start, Virus, Country, Apparent.Origin, shape, Outbreak_ID)
shp = sf::st_read("./data/spillovers/Ebola/pigott2016/index_human_case_polygon_2015.shp")
shp = shp %>%
  cbind(
    st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
  ) %>%
  dplyr::rename("ADMcode"=OUTBREAK)
poly = poly %>%
  left_join(
    shp %>% st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude), by=c("Outbreak_ID"="ADMcode")
  )

# combine
ebo1 = rbind(ebo1p, poly) %>%
  dplyr::mutate(
    Disease = "Ebola virus disease",
    IfPolygon_AdminLevel = NA, 
    Source = "Pigott2016",
    Admin_name = NA,
    ADMSource = "Pigott2016"
  ) %>%
  dplyr::rename(
    "Year"=Year.Start,
    "LocationName" = Apparent.Origin,
    "DataType"=shape,
    "ADMcode"=Outbreak_ID
  )



# -------- 2. Updates since 2015 -----------

# updated ebola data
ebo2 = read.csv("./data/spillovers/Ebola/gibb_2022/gibb_updates_2022.csv")
ebo2$ADM1code[ ebo2$ADM1code == "" ] = ebo2$ADM2code[ ebo2$ADM1code == "" ]
ebo2 = ebo2 %>% dplyr::select(-ADM2code) %>% dplyr::rename("ADMcode"=ADM1code)

# shapefile: areas and zones
shpx = sf::st_read("./data/spillovers/Ebola/gibb_2022/drc_shp/RDC_Zones de santé.shp") %>%
  dplyr::filter(Pcode %in% c(ebo2$ADMcode)) %>%
  dplyr::select(Pcode, Nom, PROVINCE) %>%
  dplyr::rename("ADMcode" = 1, "Name"=2, "Province"=3)
shpy = sf::st_read("./data/spillovers/Ebola/gibb_2022/drc_shp/RDC_Aires de santé.shp") %>%
  dplyr::filter(PCODE %in% c(ebo2$ADMcode)) %>%
  dplyr::select(PCODE, AS_, Province) %>%
  dplyr::rename("ADMcode" = 1, "Name"=2, "Province"=3)
shp2 = rbind(shpx, shpy)

# get district centroids
shp2 = shp2 %>%
  cbind(
    st_coordinates(st_centroid(shp2)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
    )

# add into ebola data
foo = ebo2 %>%
  dplyr::filter( DataType == "polygon" ) %>%
  dplyr::select(-Longitude, -Latitude) %>%
  left_join(
    shp2 %>% st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)
  ) 
ebo2 = rbind(ebo2[ ebo2$DataType == "point", ], foo)

# harmonise names with other ebola data
ebo2 = ebo2 %>%
  dplyr::rename("Year"=1) %>%
  dplyr::mutate(Disease = "Ebola virus disease",
                Virus = "Ebola Zaire",
                IfPolygon_AdminLevel=NA) %>%
  dplyr::select(names(ebo1))



# -------- combine both datasets and combine shapefiles ------------

ebo = rbind(ebo1, ebo2) %>%
  dplyr::mutate(
    ADMcode = replace(ADMcode, DataType=="point", NA)
  ) %>%
  dplyr::arrange(Year)

ebo = ebo %>%
  dplyr::rename("Pathogen"=Virus) %>%
  dplyr::mutate(
    Pathogen = replace(Pathogen, Pathogen == "Ebola Zaire", "Zaire ebolavirus"),
    Pathogen = replace(Pathogen, Pathogen == "Ebola Sudan", "Sudan ebolavirus"),
    Pathogen = replace(Pathogen, Pathogen == "Ebola Ivory Coast", "Tai Forest ebolavirus"),
    Pathogen = replace(Pathogen, Pathogen == "Ebola Bundibugyo", "Bundibugyo ebolavirus")
  ) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, Admin_name, LocationName, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

shp = rbind(
  shp, 
  shp2 %>% dplyr::select(ADMcode, Longitude, Latitude)
)
shp$ADMarea = as.vector(sf::st_area(shp)/10^6)

# save
write.csv(ebo, "./output/spillovers_processed/spillovers_ebola.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_ebola_shp.R")





# # viz
# library(maptools)
# data("wrld_simpl")
# ws = sf::st_as_sf(wrld_simpl)
# ggplot() + 
#   geom_sf(data=ws, fill="grey80", color=NA) + 
#   maptheme + 
#   geom_point(data=ebo, aes(x=Longitude, y=Latitude), col="black", size=0.25)

