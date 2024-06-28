

# ==================== Mayaro virus =======================

# matches records to GAUL polygons
# outstanding issue: need to properly harmonise admin3 polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")

# shapefile for ref
shp_ref = sf::st_read("./data/spillovers/Chagas/brazil_moh/Brazil_shp_harm_2022.shp")


# ======================= Data wrangling ====================

may = read.csv("./data/spillovers/Mayaro/celone2023/S1_MAYV_compendium_georefs_29JAN23.csv") %>%
  dplyr::filter(Host_Type == "Human") %>%
  dplyr::filter(Year_MAYV_End >= 1980) %>%
  dplyr::select(Location_ID, Year_MAYV_End, Adm1, Adm2, GAUL_code, Finer_Res, Location_Type, Admin_Level, X_Coord, Y_Coord,
               Uncertainty_km, Diagnostic_Test, Positive_n) %>%
  dplyr::filter(Admin_Level != 0) %>%
  dplyr::rename(ID = 1, Year = 2, Admin_name = 3, LocationName= 4, ADMcode = 5, IfPolygon_AdminLevel = 8, 
                Longitude = 9, Latitude = 10, Buffer_radius = 11, DiagnosticMethod = 12, NumCases = 13) %>%
  dplyr::mutate(DataType = ifelse(Location_Type == "Point", "point", "polygon")) %>%
  dplyr::mutate(Disease = "Mayaro fever",
                Pathogen = "Mayaro alphavirus (MAYV)", 
                Source = "Celone et al 2023",
                SourceType = "Compiled from literature and surveillance reports",
                CasesMetric = "Outbreak/case locations",
                ADMSource = "GAUL",
                ADM1code = replace(ADMcode, IfPolygon_AdminLevel != 1, NA),
                ADM2code = replace(ADMcode, IfPolygon_AdminLevel != 2, NA),
                ADM3code = replace(ADMcode, IfPolygon_AdminLevel != 3, NA)) %>%
  dplyr::mutate(Longitude = as.numeric(Longitude),
                Latitude = as.numeric(Latitude))
  


# --------- 1. points -------------

# points defined as those with uncertainty <5km 
m1 = may %>%
  dplyr::filter(DataType == "point") %>%
  dplyr::mutate(LocationName = paste(LocationName, Finer_Res, sep=", ")) %>%
  dplyr::mutate(IfPolygon_AdminLevel = "", 
                ADMcode = "")



# ------- 2. admin polygons level 1 --------

# GAUL polys (ADM1)
adm1_polys = may %>% dplyr::filter(IfPolygon_AdminLevel == 1) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")

# extract (12/12)
gps = unique(adm1_polys$ADM1code)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

adm1_polys = adm1_polys %>%
  dplyr::mutate(ADMcode = paste("GAUL1", ADMcode, sep="_"))

# ------- 3. admin polygons level 2 --------

# GAUL polys (ADM2)
adm2_polys = may %>% dplyr::filter(IfPolygon_AdminLevel == 2) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")

# extract (77/77)
gps = unique(adm2_polys$ADM2code)
gaul2 = g1 %>% dplyr::filter(ADM2_CODE %in% gps)
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))

adm2_polys = adm2_polys %>%
  dplyr::mutate(ADMcode = paste("GAUL2", ADMcode, sep="_"))

# ------- 4. admin polygons level 3 --------

# GADM polys for level 3 
# cross-ref via lat-lon column
# adm3_polys = may %>% dplyr::filter(IfPolygon_AdminLevel == 3) 
# g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer = "level3")
# 
# # missing with lat-lons - extract polygon directly
# missing = adm3_polys
# missingx = missing; coordinates(missingx) = ~Longitude+Latitude
# missingx = st_as_sf(missingx)
# st_crs(missingx) = st_crs(g1)
# sf::sf_use_s2(FALSE)
# ii = st_intersects(missingx, g1)
# g1_missing = g1[ as.numeric(ii), ]
# adm3_polys$ADM3code = g1_missing$ID_3
# adm3_polys$ADMSource = "GADM"
# 
# # 
# gps = unique(adm3_polys$ADM3code)
# gaul3 = g1 %>% dplyr::filter(ID_3 %in% gps)
# print(paste("Found", nrow(gaul3), "of", length(gps), "GADM3 polygons", sep=" "))

# set adm3 polys to point for now
adm3_polys = may %>% 
  dplyr::filter(IfPolygon_AdminLevel == 3) %>%
  dplyr::mutate(DataType = "point")


# ------- 5. custom polygons -------

# polygons
cps = sf::st_read("./data/spillovers/Mayaro/celone2023/S4_MAYV_compendium_custom_polys_27Oct22/MAYV_custom_polys_27SEP22.shp")

custom_polys = may %>%
  dplyr::filter(ID %in% cps$Location_i) %>%
  dplyr::mutate(ADMcode = ID) %>% 
  dplyr::mutate(IfPolygon_AdminLevel = "custom")



# ----- combine data ------

# data
may1 = m1 %>%
  rbind(adm1_polys) %>%
  rbind(adm2_polys) %>%
  rbind(adm3_polys) %>%
  rbind(custom_polys) %>%
  dplyr::select(-ID, -Finer_Res, -Location_Type, -ADM1code, -ADM2code, -ADM3code)


# shapefiles 
p1 = gaul1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = gaul2 %>% 
  dplyr::mutate(ADMcode = paste("GAUL2", ADM2_CODE, sep="_"),
                Admin_name = ADM2_NAME,
                Admin_level = 2) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 

# custom
cps = cps %>%
  dplyr::filter(Location_i %in% custom_polys$ADMcode) %>%
  dplyr::mutate(ADMcode = Location_i, 
                Admin_name = Name, 
                Admin_level = "custom") %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) %>%
  sf::st_transform(st_crs(p2)) %>%
  sf::st_zm(drop=TRUE)

# combine
shp = rbind(p1, p2) %>%
  rbind(cps)

# save data and shapefile
data.table::fwrite(may1, "./output/spillovers_processed/spillovers_mayaro.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_mayaro_shp.R")


# # viz
# library(maptools)
# data("wrld_simpl")
# ws = st_as_sf(wrld_simpl)
# ggplot() +
#   geom_sf(data=ws, fill="grey80", color=NA) +
#   maptheme +
#   geom_point(data=je1, aes(x=Longitude, y=Latitude), col="black", size=0.25)

