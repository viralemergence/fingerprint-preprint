
# ================== Marburg virus ========================

# Data to 2014 from Pigott et al 2015
# Data from 2014 onwards geolocated for this study (RG)

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")



# -------------- Pigott data and shapefiles ------------------

mar = read.csv("./data/spillovers/Marburg/pigott2015/occ_data.csv") %>%
  dplyr::rename("LocationName"=Apparent.Origin,
                "Year"=Year.Start,
                "Latitude"=Lat, 
                "Longitude"=Long,
                "DataType"=shape,
                "ADMcode"=Outbreak_ID) %>%
  dplyr::mutate(Disease = "Marburg virus disease",
                Pathogen = "Marburg marburgvirus",
                NumCases = NA,
                CasesMetric = "Index case locations",
                DiagnosticMetric = "Not specified",
                ADMSource = "Pigott2015",
                Source = "Pigott2015",
                SourceType = "Compiled from literature and surveillance reports", 
                IfPolygon_AdminLevel = "custom",
                Admin_name = "Not specified (see LocationName)") %>%
  dplyr::select(Disease, Pathogen, Country, Year, Longitude, Latitude, Admin_name, LocationName, 
                DataType, NumCases, CasesMetric, DiagnosticMetric, 
                ADMcode, ADMSource, IfPolygon_AdminLevel,
                Source, SourceType)

# shapefiles
shp = sf::st_read("./data/spillovers/Marburg/pigott2015/occ_polygons.shp") %>%
  dplyr::select(outbreakid, Lat_Long_N) %>%
  dplyr::rename(ADMcode = outbreakid,
                Locality_information = Lat_Long_N)




# ------------- Outbreaks post-2014 ---------------

# not including Equatorial Guinea as geographical information not precise and situation evolving
mar2 = read.csv("./data/spillovers/Marburg/gibb2023/marv_locations_2017_2023.csv")

# polygons

# level 4 uganda
gx = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level4")
gx = gx[ gx$ID_4 == "UGA.20.2.3.2_1", ]
gx = gx %>%
  dplyr::mutate(
    ADMcode = ID_4, 
    Admin_name = NAME_4,
    Admin_level = 4
  ) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level)
sf::st_geometry(gx) = "geometry"
gx1 = gx

# level 2 ghana
gx = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level2")
gx = gx[ gx$ID_2 == "GHA.1.1_1", ]
gx = gx %>%
  dplyr::mutate(
    ADMcode = ID_2, 
    Admin_name = NAME_2,
    Admin_level = 2
  ) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level)
sf::st_geometry(gx) = "geometry"
gx2 = gx

# level 2 ghana gaul
gx = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")
gx = gx[ gx$ADM2_CODE == 190580, ]
gx = gx %>%
  dplyr::mutate(
    ADMcode = ADM2_CODE, 
    Admin_name = ADM2_NAME,
    Admin_level = 2
  ) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level)
gx3 = gx

# level 2 tanzania gaul
gx = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")
gx = gx[ gx$ADM2_CODE == 48404, ]
gx = gx %>%
  dplyr::mutate(
    ADMcode = ADM2_CODE, 
    Admin_name = ADM2_NAME,
    Admin_level = 2
  ) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level)
gx4 = gx

# shapes combine
shp2 = 
  rbind(gx1, gx2) %>%
  rbind(gx3) %>%
  rbind(gx4) %>%
  dplyr::mutate(Locality_information = NA)



# ============== combine everything and save ===============

shp = shp %>% 
  dplyr::mutate(Admin_name = NA, Admin_level = NA) %>%
  rbind(shp2)

mar = rbind(mar, mar2)

write.csv(mar, "./output/spillovers_processed/spillovers_marv.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_marv_shp.R")

