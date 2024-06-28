

# ==================== Japanese encephalitis =======================

# matches records to GAUL polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")


# ======================= Data wrangling ====================

je = read.csv("./data/spillovers/JEV/occurrence_clean_updated_diagnostics_checked.csv") %>%
  dplyr::select(longitude, latitude, year, country, geometry_type, diagnostic_method,
                admin_level, GAUL) %>%
  dplyr::rename("Longitude"=longitude, "Latitude"=latitude, 
                "Year"=year, "Country"=country, "DataType"=geometry_type,
                "IfPolygon_AdminLevel"=admin_level,
                "ADMcode"=GAUL,
                "DiagnosticMetric"=diagnostic_method) %>%
  dplyr::mutate(IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, IfPolygon_AdminLevel== -999, NA),
                DiagnosticMetric = replace(DiagnosticMetric, DiagnosticMetric %in% c("serological", "Serological", "Seroogical"), "Serology"),
                DiagnosticMetric = replace(DiagnosticMetric, DiagnosticMetric %in% c("genotyping", "Genotyping", "PCR"), "Sequencing"),
                DiagnosticMetric = replace(DiagnosticMetric, DiagnosticMetric %in% c("reported", "Reported", "Reported "), "Reported"),
                DiagnosticMetric = replace(DiagnosticMetric, !DiagnosticMetric %in% c("Serology", "Sequencing", "Reported"), "Unspecified/unknown")) %>%
  dplyr::mutate(Disease = "Japanese encephalitis",
                Pathogen = "Japanese encephalitis virus (JEV)",
                Source = "Shearer/Longbottom",
                NumCases = NA,
                CasesMetric = "Outbreak location",
                SourceType = "Compiled from literature and surveillance reports (IHME)", 
                ADMSource = "GAUL",
                ADM1code = replace(ADMcode, IfPolygon_AdminLevel != 1, NA),
                ADM2code = replace(ADMcode, IfPolygon_AdminLevel != 2, NA),
                ADM3code = replace(ADMcode, IfPolygon_AdminLevel != 3, NA))

# exclude country level records (too imprecise)
je = je[ -which(je$IfPolygon_AdminLevel == 0), ]



# ------ polygons --------

# GAUL polys (ADM1)
adm1_polys = je %>% dplyr::filter(IfPolygon_AdminLevel == 1) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")

# missing with lat-lons - extract polygon directly
missing = adm1_polys %>% dplyr::filter(!ADM1code %in% g1$ADM1_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
#sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADM1code = g1_missing$ADM1_CODE

# missing without latlons
missing2 = adm1_polys %>% dplyr::filter(!ADM1code %in% g1$ADM1_CODE) %>% dplyr::filter(is.na(Longitude))

# combine
adm1_polys = rbind(adm1_polys %>% dplyr::filter(ADM1code %in% g1$ADM1_CODE),
                   missing, 
                   missing2)

# extract (194/195)
gps = unique(adm1_polys$ADM1code)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

# fix country name for Taiwan
adm1_polys$Country[ adm1_polys$ADMcode == 925 ] = "Taiwan"




# GAUL polys (ADM2)
adm2_polys = je %>% dplyr::filter(IfPolygon_AdminLevel == 2) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")

# missing with lat-lons - extract polygon directly
missing = adm2_polys %>% dplyr::filter(!ADM2code %in% g1$ADM2_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
#sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADM2code = g1_missing$ADM2_CODE

# missing without latlons
missing2 = adm2_polys %>% dplyr::filter(!ADM2code %in% g1$ADM2_CODE) %>% dplyr::filter(is.na(Longitude))

# combine
adm2_polys = rbind(adm2_polys %>% dplyr::filter(ADM2code %in% g1$ADM2_CODE),
                   missing, 
                   missing2)

# extract (210/210)
gps = unique(adm2_polys$ADM2code)
gaul2 = g1 %>% dplyr::filter(ADM2_CODE %in% gps)
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))



# GADM polys for level 3 
# cross-ref via lat-lon column
adm3_polys = je %>% dplyr::filter(IfPolygon_AdminLevel == 3) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer = "level3")

# missing with lat-lons - extract polygon directly
missing = adm3_polys
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
# sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
adm3_polys$ADM3code = g1_missing$ID_3
adm3_polys$ADMSource = "GADM"

# 
gps = unique(adm3_polys$ADM3code)
gaul3 = g1 %>% dplyr::filter(ID_3 %in% gps)
print(paste("Found", nrow(gaul3), "of", length(gps), "GADM3 polygons", sep=" "))



# combine all data
je1 = adm1_polys %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1code, sep="_")) %>%
  rbind(adm2_polys %>% dplyr::mutate(ADMcode = paste("GAUL2", ADM2code, sep="_"))) %>%
  rbind(je %>% dplyr::filter(DataType == "point") %>% dplyr::mutate(ADMcode = "")) %>%
  rbind(adm3_polys %>% dplyr::mutate(ADMcode = paste("GADM3", ADM3code, sep="_"))) %>%
  dplyr::select(-ADM1code, -ADM2code, -ADM3code)

# all records retained? TRUE
nrow(je1) == nrow(je)


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
p3 = gaul3 %>%
  dplyr::mutate(ADMcode = paste("GADM3", ID_3, sep="_"),
                Admin_name = NAME_3,
                Admin_level = 3) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level)
sf::st_geometry(p3) = "geometry"

shp1 = rbind(p1, p2) %>%
  rbind(p3)
all(shp1$ADMcode %in% je1$ADMcode)

# finally order cols and save
je1 = je1 %>%
  dplyr::select(
    Disease, Pathogen, DataType, Longitude, Latitude, Year, Country, NumCases, 
    CasesMetric, DiagnosticMetric, IfPolygon_AdminLevel,
    ADMcode, ADMSource, Source, SourceType
  )

# save data and shapefile
data.table::fwrite(je1, "./output/spillovers_processed/spillovers_jev.csv", row.names=FALSE)
save(shp1, file="./output/spillovers_processed/spillovers_jev_shp.R")


# # viz
# library(maptools)
# data("wrld_simpl")
# ws = st_as_sf(wrld_simpl)
# ggplot() +
#   geom_sf(data=ws, fill="grey80", color=NA) +
#   maptheme +
#   geom_point(data=je1, aes(x=Longitude, y=Latitude), col="black", size=0.25)

