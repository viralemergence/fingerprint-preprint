

# ================ Monkeypox (Frame/Pigott) ======================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")

# shapefile for ref
shp_ref = sf::st_read("./data/spillovers/Chagas/brazil_moh/Brazil_shp_harm_2022.shp")



# ============= Monkeypox data =================

# known or possible MPX spillovers without imported cases
mpx = read.csv("./data/spillovers/Monkeypox/frame2022/SciData_monkeypox_newIDs.csv") %>%
  dplyr::filter(organism_type == "human" & transmission_route %in% c("zoonotic", "unspecified") & !patient_type %in% c("import", "secondary")) %>%
  dplyr::filter(country != "USA") %>%
  dplyr::filter(!year %in% c("1972", "1976", "1977"))



# ----------- 1. Point data --------------

# --- 1. points with no buffer - can use as-is
m1 = mpx %>%
  dplyr::filter(shape_type == "point" & poly_type == "") %>%
  dplyr::mutate(ADMcode = NA,
                IfPolygon_AdminLevel=NA)

m1 = m1 %>%
  dplyr::rename("Longitude"=long, "Latitude"=lat)


# --- 2. points with uncertainty buffer: create buffer polygons (1 dec degree == 111km at equator)
m2 = mpx %>%
  dplyr::filter(shape_type == "polygon" & poly_type == "buffer") %>%
  dplyr::mutate(IfPolygon_AdminLevel="buffer")
m2$ADMcode = paste("MPXbuffer", 1:nrow(m2), sep="_")
for_buf = m2 %>% dplyr::select(ADMcode, long, lat, buffer_radius)
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

m2 = m2 %>%
  dplyr::rename("Longitude"=long, "Latitude"=lat) %>%
  dplyr::mutate(poly_reference = "buffer")



# ---------- 2. Polygon data ----------

# ------- a. ADMIN1 (GAUL and GADM) -----------

adm1_polys = mpx %>%
  dplyr::filter(shape_type == "polygon" & poly_field %in% c("ADM1_CODE", "GID_1")) 

# # recode missing GAULs to GADMs
# adm1_polys$poly_reference[ adm1_polys$poly_id == 2210 & adm1_polys$origin == "Aba" ] = "GADM"
# adm1_polys$poly_id[ adm1_polys$poly_id == 2210 & adm1_polys$origin == "Aba" ] = "NGA.1_1"
# adm1_polys$poly_reference[ adm1_polys$poly_id == 14920 & adm1_polys$origin == "Aben Gourrou Department" ] = "GADM"
# adm1_polys$poly_id[ adm1_polys$poly_id == 14920 & adm1_polys$origin == "Aben Gourrou Department" ] = ""

# adm1 polys from GAUL (none)
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")
gps = unique(adm1_polys$poly_id[ grep("GAUL", adm1_polys$poly_reference) ])
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

# adm1 polys from GADM
gadm = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level1")
gas = unique(adm1_polys$poly_id[ grep("GADM", adm1_polys$poly_reference) ])
gadm1 = gadm %>% dplyr::filter(ID_1 %in% gas)
print(paste("Found", nrow(gadm1), "of", length(gas), "GADM1 polygons", sep=" "))

# gadm ref
adm1_polys$poly_reference = "GADM1"
adm1_polys$IfPolygon_AdminLevel = 1
adm1_polys$ADMcode = adm1_polys$poly_id
  

# ------ b. ADMIN2 (GAUL and GADM) --------

adm2_polys = mpx %>%
  dplyr::filter(shape_type == "polygon" & poly_field %in% c("ADM2_CODE", "GID_2")) 

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

adm2_polys = adm2_polys %>% dplyr::mutate(
  poly_reference = replace(poly_reference, grepl("GAUL", poly_reference), "GAUL2"),
  poly_reference = replace(poly_reference, grepl("GADM", poly_reference), "GADM2"))
adm2_polys$poly_id[ grep("GAUL", adm2_polys$poly_reference) ] = paste(adm2_polys$poly_reference[ grep("GAUL", adm2_polys$poly_reference) ], adm2_polys$poly_id[ grep("GAUL", adm2_polys$poly_reference) ], sep="_")
adm2_polys$IfPolygon_AdminLevel = 2
adm2_polys$ADMcode = adm2_polys$poly_id


# # ------ b. ADMIN3 (GAUL recode to GADM) --------

adm3_polys = mpx %>%
  dplyr::filter(shape_type == "polygon" & poly_field %in% c("ADM3_CODE"))

# centroids
cc = read.csv("./data/spillovers/Monkeypox/frame2022/ADM3_centroids_occ_ID.csv")
adm3_polys = left_join(adm3_polys, cc)
adm3_polys = adm3_polys %>%
  dplyr::filter(!is.na(x)) %>%
  dplyr::rename("Longitude"=x, "Latitude"=y) %>%
  dplyr::mutate(shape_type = "point",
                IfPolygon_AdminLevel = 3,
                ADMcode = paste("ADM3_point", 1:length(Longitude))) %>%
  dplyr::select(-long, -lat)



# ----------- harmonise shapfiles and merge to match database -------------

# combine all polys
mpx_polys = adm1_polys %>%
  rbind(adm2_polys) 

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
shp1 = p2

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
shp2 = rbind(p1, p2)
st_geometry(shp2) = "geometry"

# combine
shp = rbind(shp1, shp2)

# everything cross-references ok
#all(shp$ADMcode %in% mpx_polys$poly_id)

# add centroid lat-lon for polygons
ct = st_coordinates(st_centroid(shp)) %>%
  as.data.frame() %>%
  dplyr::rename("Longitude"=X, "Latitude"=Y)
shp = cbind(shp, ct)
mpx_polys = mpx_polys %>%
  dplyr::select(-long, -lat) %>%
  left_join(
    shp %>% st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)
  )

# combine with adm3 points
mpx_polys = rbind(mpx_polys, adm3_polys)


# -------- add custom shapefiles ------

m5 = mpx %>% 
  dplyr::filter(poly_type == "custom") %>%
  dplyr::mutate(IfPolygon_AdminLevel = NA)

# custom_polys
cp = unique(m5$poly_reference)

# get polys, combine and align IDs
shp_custom = data.frame()
res = data.frame()
for(i in cp){
  
  res_x = m5 %>% dplyr::filter(poly_reference == i)
  poly_x = list.files("./data/spillovers/Monkeypox/frame2022/custom_polygons/", pattern=".shp", full.names=TRUE)
  poly_x = poly_x[ -grep(".xml", poly_x) ]
  poly_x = poly_x[ grep(i, poly_x)]
  poly_x = sf::st_read(poly_x)
  
  if(is.na(sf::st_crs(poly_x))){
    st_crs(poly_x) = st_crs(shp)
  }
  poly_x = sf::st_transform(poly_x, st_crs(shp))
  
  if(! c("GAUL_CODE") %in% names(poly_x)){ poly_x$GAUL_CODE = unique(res_x$poly_id) }
  poly_x$ADMcode = paste(i, poly_x$GAUL_CODE, sep="_")
  res_x$ADMcode = paste(i, res_x$poly_id, sep="_")
  print(all(res_x$ADMcode %in% poly_x$ADMcode))
  poly_x = poly_x %>%
    dplyr::select(ADMcode) %>%
    cbind(
      st_coordinates(st_centroid(poly_x)) %>%
        as.data.frame() %>%
        dplyr::rename("Longitude"=X, "Latitude"=Y)
    )
  res_x = res_x %>%
    dplyr::select(-long, -lat) %>%
    left_join(
      poly_x %>% st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)
    )
  shp_custom = rbind(shp_custom, poly_x)
  res = rbind(res, res_x)
  
}
m5 = res




# --------- combine everything ----------

# combine
mpx_pts = m1 %>% 
  dplyr::mutate(DataType = "point")
mpx_polys = m2 %>% 
  rbind(mpx_polys) %>%
  rbind(m5) %>%
  dplyr::mutate(DataType = shape_type)
mpx = rbind(mpx_pts, mpx_polys)

# full dataframe
mpx = mpx %>%
  dplyr::mutate(Disease = "Monkeypox", Pathogen = Hmisc::capitalize(pathogen), Year = year_start,
                Country = countrycode::countrycode(country, origin="iso3c", destination="country.name"),
               LocationName = origin,
                ADMSource = poly_reference, Source="Frame2022",
               CasesMetric = "Geolocated outbreak or case locations",
               DiagnosticMetric = diagnostic, 
               NumCases = NA) %>%
  dplyr::select(Disease, Pathogen, DataType, Year, Country, LocationName, Longitude, Latitude, NumCases, CasesMetric, DiagnosticMetric, IfPolygon_AdminLevel, ADMcode,
                ADMSource, Source) %>%
  dplyr::mutate(Pathogen = paste(Pathogen, "virus", sep=" "))  

# combined shapefile
shp_processed = buf %>% 
  sf::st_transform(st_crs(shp)) %>%
  dplyr::left_join(mpx[ , c("ADMcode", "Longitude", "Latitude")]) %>% 
  dplyr::mutate(Admin_level = NA) %>% 
  dplyr::select(ADMcode, Longitude, Latitude, Admin_level) %>%
  rbind(
    shp %>% dplyr::select(ADMcode, Longitude, Latitude, Admin_level) 
  ) %>%
  rbind(
    shp_custom %>% dplyr::mutate(Admin_level = NA)
  )
all(shp_processed$ADMcode %in% mpx$ADMcode) # TRUE
shp = shp_processed

# save
write.csv(mpx, "./output/spillovers_processed/spillovers_mpx.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_mpx_shp.R")


