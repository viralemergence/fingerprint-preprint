

# ================ Rift Valley fever ======================

# IHME dataset for RVF: process to link observations to polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")

# # install miniconda / setup python dependencies
# # numpy and ee (earth engine api for python)
#rgee::ee_install()
#rgee::ee_install_upgrade()
#ee_check()

# # 1. Initialize the Python Environment  
#ee_Initialize()




# -------------- IHME dataset of RVF occurrences -------------------

# acute human index cases (i.e. spillovers)
rvf = read.csv("./data/spillovers/RVF/IHME_RVF_OCCURRENCE_Y2020M02D28.CSV") %>%
  dplyr::filter(organism_type == "Human" & patient_type != "import") %>% # autochthonous human cases
  dplyr::filter(is.na(serosurvey) | serosurvey == "diagnostic") %>% # remove exploratory serosurveys
  dplyr::select(pathogen, transmission_route, diagnostic, country, origin, problem_geography, lat, long, latlong_source,
                loc_confidence, shape_type, poly_type, poly_reference, buffer_radius, poly_id, poly_field, year_start) 

# remove country level records
rvf = rvf[ -which(rvf$poly_field == "ADM0_CODE"), ]



# ------------- ADMIN1 POLYGONS ---------------

# recode GADM polys
rvf$poly_reference[ rvf$poly_id == 48360 & rvf$origin == "Iringa Region" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 48360 & rvf$origin == "Iringa Region" ] = "TZA.5_1"
rvf$poly_reference[ rvf$poly_id == 115001 & rvf$origin == "Arusha Region" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 115001 & rvf$origin == "Arusha Region" ] = "TZA.1_1"
rvf$poly_reference[ rvf$poly_id == 48370 & rvf$origin == "Mwanza Region" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 48370 & rvf$origin == "Mwanza Region" ] = "TZA.16_1"
rvf$poly_reference[ rvf$poly_id == 12226 & rvf$origin == "Mahdia" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 12226 & rvf$origin == "Mahdia" ] = "TUN.12_1"
rvf$poly_reference[ rvf$poly_id == 11153 & rvf$origin == "Tagant" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 11153 & rvf$origin == "Tagant" ] = "MRT.11_1"
rvf$poly_reference[ rvf$poly_id == 12231 & rvf$origin == "Kabale District" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 12231 & rvf$origin == "Kabale District" ] = "UGA.12_1"

# select
adm1_polys = rvf %>%
  dplyr::filter(poly_field == "ADM1_CODE") %>%
  dplyr::select(country, origin, poly_field, poly_id, poly_reference) %>% 
  distinct()

# adm1 polys from GAUL
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")
gps = unique(adm1_polys$poly_id[ grep("GAUL", adm1_polys$poly_reference) ])
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

# adm1 polys from GADM
gadm = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level1")
gas = unique(adm1_polys$poly_id[ grep("GADM", adm1_polys$poly_reference) ])
gadm1 = gadm %>% dplyr::filter(ID_1 %in% gas)
print(paste("Found", nrow(gadm1), "of", length(gas), "GADM1 polygons", sep=" "))



# ------------- ADMIN2 POLYGONS ---------------

# recode GADM polys
rvf$poly_reference[ rvf$poly_id == 18543 & rvf$origin == "Garissa District" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 18543 & rvf$origin == "Garissa District" ] = "KEN.7.4_1"
rvf$poly_reference[ rvf$poly_id == 48457 & rvf$origin == "Magu" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 48457 & rvf$origin == "Magu" ] = "TZA.16.4_1"
rvf$poly_reference[ rvf$poly_id == 48490 & rvf$origin == "Iramba" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 48490 & rvf$origin == "Iramba" ] = "TZA.25.2_1"
rvf$poly_reference[ rvf$poly_id == 48393 & rvf$origin == "Kondoa" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 48393 & rvf$origin == "Kondoa" ] = "TZA.3.5_1"
rvf$poly_reference[ rvf$poly_id == 48444 & rvf$origin == "Kilosa" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 48444 & rvf$origin == "Kilosa" ] = "TZA.14.3_1"
rvf$poly_reference[ rvf$poly_id == 3008153 & rvf$origin == "Tamchekett/Hodh El Gharbi" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 3008153 & rvf$origin == "Tamchekett/Hodh El Gharbi" ] = "MRT.8.3_1"
rvf$poly_reference[ rvf$poly_id == 5003153 & rvf$origin == "Maghta-Lahjar/Brakna" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 5003153 & rvf$origin == "Maghta-Lahjar/Brakna" ] = "MRT.3.5_1"
rvf$poly_reference[ rvf$poly_id == 3013153 & rvf$origin == "Mederdra/Trarza" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 3013153 & rvf$origin == "Mederdra/Trarza" ] = "MRT.13.3_1"
rvf$poly_reference[ rvf$poly_id == 3007153 & rvf$origin == "Djigueni/Hodh Charghi" ] = "GADM"
rvf$poly_id[ rvf$poly_id == 3007153 & rvf$origin == "Djigueni/Hodh Charghi" ] = "MRT.7.3_1"

# select
adm2_polys = rvf %>%
  dplyr::filter(poly_field == "ADM2_CODE") %>%
  dplyr::select(country, origin, poly_field, poly_id, poly_reference) %>% 
  distinct()

# adm2 polys from GAUL
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")
gps = unique(adm2_polys$poly_id[ grep("GAUL", adm2_polys$poly_reference) ])
gaul2 = g1 %>% dplyr::filter(ADM2_CODE %in% as.numeric(gps))
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))

# adm2 polys from GADM
gadm = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level2")
gas = unique(adm2_polys$poly_id[ grep("GADM", adm2_polys$poly_reference) ])
gadm2 = gadm %>% dplyr::filter(ID_2 %in% gas)
print(paste("Found", nrow(gadm2), "of", length(gas), "GADM1 polygons", sep=" "))



# ---------- adm3 polygons: map to GADM ------------

rvf$poly_reference[ rvf$poly_id == 51609 & rvf$origin == "Fafi-Jarajila" ] = "GADM3"
rvf$poly_id[ rvf$poly_id == 51609 & rvf$origin == "Fafi-Jarajila" ] = "KEN.7.3.4_1"
rvf$poly_reference[ rvf$poly_id == 51446 & rvf$origin == "Bahari" ] = "GADM3"
rvf$poly_id[ rvf$poly_id == 51446 & rvf$origin == "Bahari" ] = "KEN.21.2.1_1"
rvf$poly_reference[ rvf$poly_id == 51723 & rvf$origin == "Keserian, Makutani, Baringo" ] = "GADM3"
rvf$poly_id[ rvf$poly_id == 51723 & rvf$origin == "Keserian, Makutani, Baringo" ] = "KEN.1.4.4_1"

adm3_polys = rvf %>%
  dplyr::filter(poly_field == "ADM3_CODE") %>%
  dplyr::select(country, origin, poly_field, poly_id, poly_reference) %>% 
  distinct()

# adm3 polys from GADM
gadm = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level3")
gas = unique(adm3_polys$poly_id[ grep("GADM", adm3_polys$poly_reference) ])
gadm3 = gadm %>% dplyr::filter(ID_3 %in% gas)
print(paste("Found", nrow(gadm3), "of", length(gas), "GADM3 polygons", sep=" "))



# ----------- harmonise shapfiles and merge to match database -------------

# add gaul names
rvf = rvf %>%
  dplyr::mutate(
    poly_reference = replace(poly_reference, grepl("GAUL", poly_reference) & poly_field == "ADM1_CODE", "GAUL1"),
    poly_reference = replace(poly_reference, grepl("GAUL", poly_reference) & poly_field == "ADM2_CODE", "GAUL2"),
    poly_reference = replace(poly_reference, grepl("GADM", poly_reference) & poly_field == "ADM1_CODE", "GADM1"),
    poly_reference = replace(poly_reference, grepl("GADM", poly_reference) & poly_field == "ADM2_CODE", "GADM2"),
    poly_reference = replace(poly_reference, grepl("GAUL_2015", poly_reference) & poly_field == "ADM3_CODE", "GADM3")
  )
rvf$poly_id[ grep("GAUL", rvf$poly_reference) ] = paste(rvf$poly_reference[ grep("GAUL", rvf$poly_reference) ], rvf$poly_id[ grep("GAUL", rvf$poly_reference) ], sep="_")

# combine and harmonise shapefiles 

# gaul
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

# gadm
p1 = gadm1 %>% 
  dplyr::mutate(ADMcode = ID_1,
                Admin_name = NAME_1,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = gadm2 %>% 
  dplyr::mutate(ADMcode = ID_2,
                Admin_name = NAME_2,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p3 = gadm3 %>% 
  dplyr::mutate(ADMcode = ID_3,
                Admin_name = NAME_3,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
shp2 = rbind(p1, p2) %>%
  rbind(p3)
st_geometry(shp2) = "geometry"

# combine
shp = rbind(shp1, shp2)



# -------- add custom shapefiles ------

rr = unique(rvf$poly_reference)
cps = rr[ !is.na(rr) & !grepl("GADM|GAUL", rr)]

cust = rvf %>% dplyr::filter(poly_reference %in% cps) %>%
  dplyr::mutate(poly_id = poly_reference)

# read in shapefiles one at a time
cust_shps = list.files("./data/spillovers/RVF/custom_polygons/", pattern=".shp", full.names = TRUE)
shp3 = data.frame()
for(i in 1:nrow(cust)){
  cs1 = sf::st_read( cust_shps[ grep(cust$poly_reference[i], cust_shps)] )
  st_crs(cs1) = st_crs(shp)
  cs1 = cs1 %>%
    dplyr::mutate(ADMcode = cust$poly_id[i],
                  Admin_name = cust$poly_id[i],
                  Admin_level = "custom") %>%
    dplyr::select(ADMcode, Admin_name, Admin_level)
  shp3 = rbind(shp3, cs1)
}

# combine tweaked results with rvf data
rvf = rbind(
  rvf %>% dplyr::filter(!poly_reference %in% cps),
  cust
)

# combine shapefiles
shp = rbind(shp, shp3)



# --------- final processing and save ---------

# add lat-lons
sf::sf_use_s2(FALSE)
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
)
sf::sf_use_s2(TRUE)

# add in lat lons of polygon centroids
rvf = left_join(rvf, shp %>% st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude), by=c("poly_id"="ADMcode"))
rvf$lat[ is.na(rvf$lat) ] = rvf$Latitude[ is.na(rvf$lat) ]
rvf$long[ is.na(rvf$long) ] = rvf$Longitude[ is.na(rvf$long) ]
rvf = rvf %>% dplyr::select(-Latitude, -Longitude)

# rename "shape type" for buffer points to point
rvf$shape_type[ rvf$shape_type == "polygon" & rvf$poly_type == "buffer" ] = "point"

#rename and rearrange cols 
rvf = rvf %>%
  dplyr::mutate(Disease = "Rift Valley Fever",
                Source = "IHME",
                pathogen = "Rift Valley Fever virus") %>%
  dplyr::rename("Pathogen"=pathogen,
                "Country"=country,
                "LocationName" = origin,
                "Latitude"=lat,
                "Longitude"=long,
                "DataType" = shape_type, 
                "PolyType"=poly_type,
                "BufferRadius"=buffer_radius,
                "ADMSource"=poly_reference,
                "ADMcode"=poly_id,
                "IfPolygon_AdminLevel"=poly_field,
                "Year"=year_start) %>%
  dplyr::mutate("CasesMetric" = ifelse(diagnostic == "reported", "Case report", "Confirmed (see DiagnosticMetric column)"),
                "DiagnosticMetric" = Hmisc::capitalize(diagnostic)) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, CasesMetric, DiagnosticMetric, Longitude, Latitude, LocationName,
                IfPolygon_AdminLevel, ADMcode, ADMSource, PolyType, BufferRadius, Source)
rvf$IfPolygon_AdminLevel[ rvf$IfPolygon_AdminLevel=="ADM1_CODE" ] = 1
rvf$IfPolygon_AdminLevel[ rvf$IfPolygon_AdminLevel=="ADM2_CODE" ] = 2
rvf$IfPolygon_AdminLevel[ rvf$IfPolygon_AdminLevel=="ADM3_CODE" ] = 3

# save
write.csv(rvf, "./output/spillovers_processed/spillovers_rvf.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_rvf_shp.R")
