

# ==================== Hendra and Nipah viruses =======================

# matches records to GAUL polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")

# shapefile for ref
shp_ref = sf::st_read("./data/shapefiles/africa_shp.shp")




# ============= Initial data wrangling =============

hen = read.csv("./data/spillovers/Henipaviruses/experiment1_all_henipaviruses_2019-07-15.csv") %>%
  dplyr::filter(organism_type == "human") %>%
  dplyr::filter(poly_field != "ADM0_CODE") %>% # too imprecise
  dplyr::filter(pathogen %in% c("hendra virus", "nipah virus")) %>%
  dplyr::filter(transmission_route == "zoonotic") %>%
  dplyr::mutate(
    Disease = ifelse(pathogen == "hendra virus", "Hendra virus disease", "Nipah virus disease"),
    Pathogen = Hmisc::capitalize(pathogen),
    Year = year_start,
    LocationName = origin,
    Longitude = long,
    Latitude = lat,
    Country = countrycode::countrycode(country, origin = "iso3c", destination = "country.name"),
    DataType = shape_type,
    ADMcode = poly_id,
    CasesMetric = "Geolocated spillover/outbreak locations",
    DiagnosticMetric = diagnostic,
    IfPolygon_AdminLevel = "",
    IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, poly_type == "buffer", "buffer"),
    IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, poly_field == "ADM1_CODE", "1"),
    IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, poly_field == "ADM2_CODE", "2"),
    Source = "Pigott_unpublished",
    SourceType = "Outbreak locations from literature and surveillance reports"
  ) %>%
  dplyr::select(Disease, Pathogen, Year, Country, LocationName, Longitude, Latitude, DataType, ADMcode,
                IfPolygon_AdminLevel, CasesMetric, DiagnosticMetric, buffer_radius, Source, SourceType) %>%
  distinct()



# ================ sort geolocation ================


# ----- 1. points - no further processing ----------

hen1 = hen %>% dplyr::filter(DataType == "point")



# ---- 2. point + uncertainty buffer -----------

hen2 = hen %>%
  dplyr::filter(DataType == "polygon" & IfPolygon_AdminLevel == "buffer") 
hen2$ADMcode = paste("HENbuffer", 1:nrow(hen2), sep="_")
for_buf = hen2 %>% dplyr::select(ADMcode, Longitude, Latitude, buffer_radius)
for_buf = st_as_sf(for_buf, coords = c("Longitude", "Latitude")) 
st_crs(for_buf) = crs(shp_ref)
for_buf$radius_deg = for_buf$buffer_radius / 111 # arc degree conversion (approxs)
createBuffer = function(x){
  return(sf::st_buffer( for_buf[x, ], dist=for_buf$radius_deg[x] ))
}
buf = data.frame()
for(i in 1:nrow(for_buf)){
  bx = createBuffer(i)
  buf = rbind(buf, bx)
}


# ---- 3. polygons admin1 ----

# GAUL polys (ADM1)
adm1_polys = hen %>% dplyr::filter(IfPolygon_AdminLevel == "1") 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")

# recode to match and extract
adm1_polys$ADMcode[ adm1_polys$LocationName == "Perak" ] = 1898
adm1_polys$ADMcode[ adm1_polys$LocationName == "Negeri Sembilan" ] = 1896
adm1_polys$ADMcode[ adm1_polys$LocationName == "Selangor" ] = 1903
  
# extract (4/4)
gps = unique(adm1_polys$ADMcode)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

# add lat-lons
ll = as.data.frame(st_coordinates(sf::st_centroid(gaul1)))
ll$ADMcode = gaul1$ADM1_CODE
names(ll)[1:2] = c("Longitude", "Latitude")

adm1_polys = adm1_polys %>%
  dplyr::select(-Longitude, -Latitude) %>%
  left_join(ll) %>%
  dplyr::mutate(ADMcode = paste("GAUL1", ADMcode, sep="_"))




# --- 4. polygons admin2 ----

# GAUL polys (ADM2)
adm2_polys = hen %>% dplyr::filter(IfPolygon_AdminLevel == "2") 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")

# recode to match and extract
adm2_polys$ADMcode[ adm2_polys$LocationName == "Kushtia District" ] = 5800
adm2_polys$ADMcode[ adm2_polys$LocationName == "Port Dickson" ] = 37374

# extract (4/4)
gps = unique(adm2_polys$ADMcode)
gaul2 = g1 %>% dplyr::filter(ADM2_CODE %in% gps)
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))

# add lat-lons
ll = as.data.frame(st_coordinates(sf::st_centroid(gaul2)))
ll$ADMcode = gaul2$ADM2_CODE
names(ll)[1:2] = c("Longitude", "Latitude")

adm2_polys = adm2_polys %>%
  dplyr::select(-Longitude, -Latitude) %>%
  left_join(ll) %>%
  dplyr::mutate(ADMcode = paste("GAUL2", ADMcode, sep="_"))




# ========= combine and harmonise ============

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

# add buffers
shp1 = rbind(shp1, 
             buf %>% 
               dplyr::mutate(Admin_name = NA,
                             Admin_level = "buffer") %>%
               dplyr::select(ADMcode, Admin_name, Admin_level))

# combine together processed data
hen_final = hen1 %>%
  rbind(hen2) %>%
  rbind(adm1_polys) %>%
  rbind(adm2_polys) %>%
  dplyr::select(-buffer_radius)

# check everything retained = TRUE
nrow(hen) == nrow(hen_final)
all(shp1$ADMcode %in% hen_final$ADMcode)

# save data and shapefile
data.table::fwrite(hen_final, "./output/spillovers_processed/spillovers_henipa.csv", row.names=FALSE)
save(shp1, file="./output/spillovers_processed/spillovers_henipa_shp.R")


# # # viz
# library(maptools)
# data("wrld_simpl")
# ws = st_as_sf(wrld_simpl)
# hen_final2 = st_as_sf(hen_final, coords=c("Longitude", "Latitude")) %>%
#   dplyr::filter(Disease == "Hendra virus disease")
# st_crs(hen_final2) = st_crs(ws)
# ggplot() +
#   geom_sf(data=st_crop(ws, extent(hen_final2)+0.5), fill="grey80", color=NA) +
#   maptheme +
#   geom_sf(data=hen_final2, col="black", size=0.25)

