library(magrittr); library(dplyr)

# data from Kucharshi & Edmunds 2015
ff = read.csv("./data/spillovers/Influenza/kucharski2015/H5N1_list_2006_2014.csv")

# read in and string trim
influenza = ff %>% dplyr::filter(!Country %in% c("Egypt", "Egypt ")) %>%
  dplyr::mutate(
    District = stringr::str_trim(District),
    Province = stringr::str_trim(Province),
    District = replace(District, is.na(District), "")
  ) 

# geo locations data - cross referenced to GADM
geo = read.csv("./data/spillovers/Influenza/kucharski2015/h5n1_list_locs.csv")

# combine
influenza = left_join(influenza, geo) %>%
  dplyr::filter(data_type %in% c("point", "polygon")) %>% 
  dplyr::select(-X, -X.1, -X.2, -X.3) %>%
  dplyr::mutate(
    data_type = replace(data_type, GADM_level==3, "point")
  )


# ------- polygons admin1 --------

p1.1 = sf::st_read("./data/spillovers/Influenza/kucharski2015/shapefiles/gadm41_KHM_shp/gadm41_KHM_1.shp") %>%
  dplyr::select(COUNTRY, NAME_1, GID_1) %>%
  dplyr::filter(GID_1 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_1, "ADMcode"=GID_1) %>% 
  dplyr::mutate(Admin_level = 1)

p1.2 = sf::st_read("./data/spillovers/Influenza/kucharski2015/shapefiles/gadm41_CHN_shp/gadm41_CHN_1.shp") %>%
  dplyr::select(COUNTRY, NAME_1, GID_1) %>%
  dplyr::filter(GID_1 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_1, "ADMcode"=GID_1) %>% 
  dplyr::mutate(Admin_level = 1)

p1.3 = sf::st_read("./data/spillovers/Influenza/kucharski2015/shapefiles/gadm41_VNM_shp/gadm41_VNM_1.shp") %>%
  dplyr::select(COUNTRY, NAME_1, GID_1) %>%
  dplyr::filter(GID_1 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_1, "ADMcode"=GID_1) %>% 
  dplyr::mutate(Admin_level = 1)

p1.4 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level1") %>%
  dplyr::select(COUNTRY, NAME_1, ID_1) %>%
  dplyr::filter(ID_1 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_1, "ADMcode"=ID_1) %>% 
  dplyr::mutate(Admin_level = 1) %>%
  dplyr::filter(COUNTRY %in% c("Indonesia", "Thailand", "Laos", "Myanmar"))
sf::st_geometry(p1.4) = "geometry"

# combine
p1 = p1.1 %>%
  rbind(p1.2) %>%
  rbind(p1.3) %>%
  rbind(p1.4) 

# check all matched
all( unique(influenza$GADM_code[ influenza$GADM_level == 1 & !is.na(influenza$GADM_level)]) %in% p1$ADMcode )


# -------- polygons admin2 -----------

p2.1 = sf::st_read("./data/spillovers/Influenza/kucharski2015/shapefiles/gadm41_KHM_shp/gadm41_KHM_2.shp") %>%
  dplyr::select(COUNTRY, NAME_2, GID_2) %>%
  dplyr::filter(GID_2 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_2, "ADMcode"=GID_2) %>% 
  dplyr::mutate(Admin_level = 2)

p2.2 = sf::st_read("./data/spillovers/Influenza/kucharski2015/shapefiles/gadm41_CHN_shp/gadm41_CHN_2.shp") %>%
  dplyr::select(COUNTRY, NAME_2, GID_2) %>%
  dplyr::filter(GID_2 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_2, "ADMcode"=GID_2) %>% 
  dplyr::mutate(Admin_level = 2)

p2.3 = sf::st_read("./data/spillovers/Influenza/kucharski2015/shapefiles/gadm41_VNM_shp/gadm41_VNM_2.shp") %>%
  dplyr::select(COUNTRY, NAME_2, GID_2) %>%
  dplyr::filter(GID_2 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_2, "ADMcode"=GID_2) %>% 
  dplyr::mutate(Admin_level = 2)

p2.4 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level2") %>%
  dplyr::select(COUNTRY, NAME_2, ID_2) %>%
  dplyr::filter(ID_2 %in% influenza$GADM_code) %>%
  dplyr::rename("NAME"=NAME_2, "ADMcode"=ID_2) %>% 
  dplyr::mutate(Admin_level = 2) %>%
  dplyr::filter(COUNTRY %in% c("Indonesia", "Thailand", "Laos", "Myanmar"))
sf::st_geometry(p2.4) = "geometry"

# combine
p2 = p2.1 %>%
  rbind(p2.2) %>%
  rbind(p2.3) %>%
  rbind(p2.4) 

# check all matched
all( unique(influenza$GADM_code[ influenza$GADM_level == 2 & !is.na(influenza$GADM_level)]) %in% p2$ADMcode )

# save polys
shp = rbind(p1, p2) %>%
  dplyr::rename(Admin_name = NAME)

# library(ggplot2)
# p2 %>%
#   dplyr::filter(COUNTRY == "Cambodia") %>%
#   ggplot() + 
#   geom_sf(color="black", fill=NA) +
#   geom_point(data = influenza[ influenza$Country == "Cambodia" & influenza$GADM_level == 2, ], aes(Long, Lat), col="red")




# ================= format influenza data ==========================

# format for fingerprint
flu = influenza %>%
  dplyr::select(Year, Country, District, Province, District, Long, Lat, Cluster, Source, GADM_code, GADM_level, data_type) %>%
  distinct() %>%
  dplyr::rename("LocationName" = District,
                "Admin_name" = Province,
                "NumCases" = Cluster,
                "Latitude_orig" = Lat,
                "Longitude_orig" = Long) %>%
  dplyr::mutate(Disease = "Influenza H5N1", 
                Pathogen = "Influenza A H5N1",
                NumCases = replace(NumCases, is.na(NumCases), 1),
                DataType = data_type,
                IfPolygon_AdminLevel = GADM_level,
                ADMcode = GADM_code,
                ADMSource = "GADM",
                DiagnosticMetric = "Not specified",
                CasesMetric = "Confirmed cases",
                Source = paste(Source, "(Kucharski2015)", sep=" "),
                SourceType = "Geolocated case/outbreak locations, some uncertainty in geolocation method") %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Latitude_orig, Longitude_orig, Admin_name, LocationName, 
                NumCases, DiagnosticMetric, CasesMetric, IfPolygon_AdminLevel,
                ADMcode, ADMSource, Source, SourceType)

# add latlons from shapefile
shp = cbind(
  shp, sf::st_coordinates(sf::st_centroid(shp))
)
shp = shp %>%
  dplyr::rename("Longitude"=X, "Latitude"=Y)

flu = flu %>%
  dplyr::left_join(
    shp %>% dplyr::select(ADMcode, "Longitude", "Latitude") %>% sf::st_drop_geometry()
  )
flu$Longitude[ is.na(flu$Longitude) ] = flu$Longitude_orig[ is.na(flu$Longitude) ]
flu$Latitude[ is.na(flu$Latitude) ] = flu$Latitude_orig[ is.na(flu$Latitude) ]


# save
write.csv(flu, "./output/spillovers_processed/spillovers_h5n1.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_h5n1_shp.R")
