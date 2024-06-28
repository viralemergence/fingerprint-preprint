

# ================ MERS-CoV ======================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")

# shapefile for crs ref
shp_ref = sf::st_read("./data/shapefiles/africa_shp.shp")



# -------------- IHME dataset of MERS CoV occurrences -------------------

# acute human index cases (i.e. spillovers)
mers = read.csv("./data/spillovers/MERS-CoV/ihme/IHME_MERS_COV_DATABASE_1983_2017_Y2019M07D23.CSV") %>%
  dplyr::filter(organism_type == "human" & patient_type == "index" & serosurvey %in% c("", "diagnostic") & problem_geography == "") %>%
  dplyr::filter(pathogen == "MERS-CoV") %>%
  dplyr::select(pathogen, transmission_route, diagnostic, country, origin, problem_geography, lat, long, latlong_source,
                loc_confidence, shape_type, poly_type, buffer_radius, gaul_year_or_custom_shapefile, poly_id, poly_field, year_start,
                clinical, diagnostic, diagnostic_note) %>%
  dplyr::filter(gaul_year_or_custom_shapefile != "ARE_Western_Region.shp") #exclude v large custom shapefile
mers = mers[ -which(mers$poly_field == "ADM0_CODE"), ] # remove everything at country level


# 1. points
mers1 = mers %>% dplyr::filter(shape_type == "point") %>%
  dplyr::mutate(ADMcode = NA)

# 2. polygons wuth uncertainty buffer (create buffer)
mers2 = mers %>% dplyr::filter(shape_type == "polygon" & poly_type == "buffer")

mers2$ADMcode = paste("MERSbuffer", 1:nrow(mers2), sep="_")
for_buf = mers2 %>% dplyr::select(ADMcode, long, lat, buffer_radius)
coordinates(for_buf) = ~long+lat
for_buf = st_as_sf(for_buf) 
st_crs(for_buf) = crs(shp_ref)
for_buf$radius_m = for_buf$buffer_radius * 1000 # metres
createBuffer = function(x){
  return(sf::st_buffer( for_buf[x, ], dist=for_buf$radius_m[x] ))
}
buf = data.frame()
for(i in 1:nrow(for_buf)){
  bx = createBuffer(i)
  buf = rbind(buf, bx)
}


# 3. adm1 polys
mers3 = mers %>%
  dplyr::filter(shape_type == "polygon" & poly_field %in% c("ADM1_CODE")) 

# adm1 polys from GAUL (found all)
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")
gps = unique(mers3$poly_id)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

p1 = gaul1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 

mers3$ADMcode = paste("GAUL1", mers3$poly_id, sep="_")

# add long lat from polys
ct = st_coordinates(st_centroid(p1)) %>%
  as.data.frame() %>%
  dplyr::rename("long"=X, "lat"=Y)
ct$ADMcode = p1$ADMcode
mers3 = mers3 %>% dplyr::select(-long, -lat) %>% left_join(ct)


# combine all mers data
mers = mers1 %>%
  rbind(mers2) %>%
  rbind(mers3)

# combine shapefile
shp = rbind(
  p1,
  buf %>% dplyr::select(ADMcode) %>% dplyr::mutate(Admin_name = NA, Admin_level = "buffer")
)
all(shp$ADMcode %in% mers$ADMcode)


# finally fix fields
mers = mers %>%
  dplyr::mutate(Disease = "Middle East respiratory syndrome (MERS)",
                Source = "IHME",
                ADMSource = ifelse(poly_field == "ADM1_CODE", "GAUL", ""),
                IfPolygon_AdminLevel = ifelse(poly_field == "ADM1_CODE", 1, NA),
                NumCases = NA) %>%
  dplyr::rename("Pathogen"=pathogen,
                "Country"=country,
                "LocationName" = origin,
                "Latitude"=lat,
                "Longitude"=long,
                "DataType" = shape_type, 
                "Year"=year_start,
                "CasesMetric"=clinical,
                "DiagnosticMetric"=diagnostic) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, LocationName,
                IfPolygon_AdminLevel, ADMcode, ADMSource, NumCases, CasesMetric, DiagnosticMetric, Source)

# save
write.csv(mers, "./output/spillovers_processed/spillovers_mers.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_mers_shp.R")

