
# ================ P knowlesi ======================

# GAUL polygons accessed from Google Earth Engine

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")

# # install miniconda / setup python dependencies
# # numpy and ee (earth engine api for python)
#rgee::ee_install()
#rgee::ee_install_upgrade()
#ee_check()

# # 1. Initialize the Python Environment  
ee_Initialize()




# --------------- Pk occurrences from Shearer et al Plos NTDs ----------------

pk = read.csv("./data/spillovers/Pknowlesi/knowlesei_xx.csv") %>% 
  dplyr::filter(Host == "human") %>%
  dplyr::filter(Presence == 1) %>%
  dplyr::filter(Unconfirmed..1. == 0 & Training.data == 1) 

# centroid lat/lon
pk$Latitude[ is.na(pk$Latitude) ] = pk$Centroid.Latitude[ is.na(pk$Latitude) ]
pk$Longitude[ is.na(pk$Longitude) ] = pk$Centroid.Longitude[ is.na(pk$Longitude) ]

# select cols
pk = pk %>%
  dplyr::select(Year, Country, Site_name, Subnational_area, Geometry_type, Longitude, Latitude,
                Admin_level, Gaul_code, Polygon_code) %>%
  dplyr::rename("LocationName"=3, "Admin_name"=4, "DataType"=5, "IfPolygon_AdminLevel"=8, "ADMcode"=9, "Polycode"=10) %>%
  dplyr::mutate(NumCases = NA, Source="Shearer2016", Disease="Plasmodium knowlesi malaria", Pathogen="Plasmodium knowlesi", ADMSource="GAUL") %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, LocationName, Admin_name, IfPolygon_AdminLevel, 
                ADMcode, Polycode, ADMSource, Source)

# get polygons via GEE
# admin1 and admin2 polygons
adm1_polys = unique(pk$ADMcode[ pk$IfPolygon_AdminLevel == "Admin1"])
adm2_polys = unique(pk$ADMcode[ pk$IfPolygon_AdminLevel == "Admin2"])

# adm1 polys from GEE
gee1 = ee$FeatureCollection("FAO/GAUL/2015/level1")$
  filter( 
    ee$Filter$Or(
      ee$Filter$inList('ADM1_CODE', adm1_polys)
    )
  )
gee1 = ee_as_sf(gee1)

# adm2 polys from GEE
gee2 = ee$FeatureCollection("FAO/GAUL/2015/level2")$
  filter( 
    ee$Filter$Or(
      ee$Filter$inList('ADM2_CODE', adm2_polys)
    )
  )
gee2 = ee_as_sf(gee2)

# combine into full polygons
p1 = gee1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = gee2 %>% 
  dplyr::mutate(ADMcode = paste("GAUL2", ADM2_CODE, sep="_"),
                Admin_name = ADM2_NAME,
                Admin_level = 2) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
shp = rbind(p1, p2)

# harmonise admcode with polygons
pk$ADMcode[ pk$IfPolygon_AdminLevel == "Admin1" ] = paste("GAUL1", pk$ADMcode[ pk$IfPolygon_AdminLevel == "Admin1" ], sep="_")
pk$ADMcode[ pk$IfPolygon_AdminLevel == "Admin2" ] = paste("GAUL2", pk$ADMcode[ pk$IfPolygon_AdminLevel == "Admin2" ], sep="_")

# other non GAUL polygon identifiers: all Malaysia
foo = pk %>%
  dplyr::filter(grepl("para_ras", IfPolygon_AdminLevel)) %>%
  dplyr::select(LocationName, Admin_name) %>%
  distinct()
polys = sf::st_as_sf(raster::getData(name="GADM", country="MYS", level=2))
foo2 = read.csv("./data/spillovers/Pknowlesi/gibb_gadm_lookup.csv") %>% dplyr::select(-X) # manual crossref
polys = polys %>% 
  dplyr::filter(GID_2 %in% foo2$Polycode) %>%
  dplyr::select(GID_2, NAME_2) %>%
  dplyr::rename("ADMcode"=1, 
                "Admin_name"=2) %>%
  dplyr::mutate(ADMcode = paste("GADM2_", ADMcode, sep=""))
foo2 = foo2 %>%
  dplyr::mutate(ADMcode = paste("GADM2_", Polycode, sep="")) %>%
  dplyr::select(-ADMlevel, -Polycode) 
foo2$ADMcode[ foo2$ADMcode == "GADM2_" ] = ""

# combine
pk = pk %>%
  dplyr::left_join(
    foo2 %>% dplyr::rename("ADMcodeX" = ADMcode)
  )
pk$ADMcode[ !is.na(pk$ADMcodeX) ] = pk$ADMcodeX[ !is.na(pk$ADMcodeX) ]
pk = pk %>% dplyr::select(-ADMcodeX, -Polycode)

# combine polygons
shp = rbind(shp, polys %>% dplyr::mutate(Admin_level = 2))

# save
write.csv(pk, "./output/spillovers_processed/spillovers_pknowlesi.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_pknowlesi_shp.R")



# viz
library(maptools)
data("wrld_simpl")
ws = sf::st_as_sf(wrld_simpl)
ggplot() + 
  geom_sf(data=ws, fill="grey80", color=NA) + 
  maptheme + 
  geom_point(data=pk, aes(x=Longitude, y=Latitude), col="black", size=0.25)


