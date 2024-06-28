

# ==================== Zika =======================

# From Messina et al (2016) mapping and global compendium of Zika virus occurrence

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")


# ======================= Data wrangling ====================

# https://figshare.com/articles/dataset/Global_compendium_of_human_Zika_virus_occurrence/2573629

# 1. Points
zik = read.csv("./data/spillovers/Zika/Occ_Standard_JM_10032016.csv") %>% 
  dplyr::filter(Admin != 0) %>%
  dplyr::rename("IfPolygon_AdminLevel"=Admin,
                "ADMarea" = Area) %>%
  dplyr::mutate(Disease = "Zika", 
                Pathogen = "Zika virus",
                DataType = ifelse(IfPolygon_AdminLevel == -999, "point", "polygon"), 
                IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, DataType == "point", NA), 
                Country = "",
                ADMcode = replace(UniqueID, DataType == "point", NA),
                ADM1code = replace(UniqueID, IfPolygon_AdminLevel==2, NA),
                ADM2code = replace(UniqueID, IfPolygon_AdminLevel==1, NA), 
                Source = "") %>%
  dplyr::select(-UniqueID) %>%
  dplyr::mutate(ADMarea = replace(ADMarea, ADMarea == -999, NA))




# GAUL polys (ADM1)
adm1_polys = zik %>% dplyr::filter(IfPolygon_AdminLevel == 1) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")

# missing with lat-lons - extract polygon directly
missing = adm1_polys %>% dplyr::filter(!ADM1code %in% g1$ADM1_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADM1code = g1_missing$ADM1_CODE

# combine
adm1_polys = rbind(adm1_polys %>% dplyr::filter(ADM1code %in% g1$ADM1_CODE),
                   missing)

# extract (54/54)
gps = unique(adm1_polys$ADM1code)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))


# GAUL polys (ADM2)
adm2_polys = zik %>% dplyr::filter(IfPolygon_AdminLevel == 2) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")

# missing with lat-lons - extract polygon directly
missing = adm2_polys %>% dplyr::filter(!ADM2code %in% g1$ADM2_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADM2code = g1_missing$ADM2_CODE

# combine
adm2_polys = rbind(adm2_polys %>% dplyr::filter(ADM2code %in% g1$ADM2_CODE),
                   missing)

# extract (210/210)
gps = unique(adm2_polys$ADM2code)
gaul2 = g1 %>% dplyr::filter(ADM2_CODE %in% gps)
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))


# combine all data
zik1 = adm1_polys %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1code, sep="_")) %>%
  rbind(adm2_polys %>% dplyr::mutate(ADMcode = paste("GAUL2", ADM2code, sep="_"))) %>%
  rbind(zik %>% dplyr::filter(DataType == "point") %>% dplyr::mutate(ADMcode = "")) %>%
  dplyr::select(-ADM1code, -ADM2code)

# all retained? TRUE
nrow(zik1) == nrow(zik)


# combine and harmonise shapefiles 
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

shp1 = rbind(p1, p2)
all(shp1$ADMcode %in% zik1$ADMcode)


# add extra metadata and save
zik1 = zik1 %>%
  dplyr::mutate(NumCases = NA,
                CasesMetric = "Outbreak locations", 
                DiagnosticMetric = "Not specified",
                ADMSource = "GAUL",
                Source = "Messina2016",
                SourceType = "Geolocated outbreak locations from literature and surveillance") %>%
  dplyr::select(
    Disease, Pathogen, DataType, Longitude, Latitude, Year, Country, NumCases, 
    CasesMetric, DiagnosticMetric, IfPolygon_AdminLevel,
    ADMcode, ADMSource, Source, SourceType
  )


# save data and shapefile
write.csv(zik1, "./output/spillovers_processed/spillovers_zika.csv", row.names=FALSE)
save(shp1, file="./output/spillovers_processed/spillovers_zika_shp.R")


# # viz
# library(maptools)
# data("wrld_simpl")
# ws = st_as_sf(wrld_simpl)
# ggplot() + 
#   geom_sf(data=ws, fill="grey80", color=NA) + 
#   maptheme + 
#   geom_point(data=zik1, aes(x=Longitude, y=Latitude), col="black", size=0.25)

