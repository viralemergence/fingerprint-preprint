

# ==================== Chikungunya =======================


setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")



# ----------------- 1. Nsoesie et al 2017 ------------------

# chikungunya dataset
chk = read.csv("./data/spillovers/Chikungunya/nsoesie_2016/chikv_11_05_2015.csv") %>%
  dplyr::filter(precision != "COUNTRY") # remove country level records

# 1. process database
chk1 = chk %>% 
  dplyr::select(occurrence.date, longitude, latitude, location.name, precision, GAUL_code) %>%
  dplyr::rename("Year" = occurrence.date,
                "Longitude"=longitude, 
                "Latitude"=latitude,
                "LocationName"=location.name,
                "ADMcode"=GAUL_code) %>%
  dplyr::mutate(Disease = "Chikungunya", 
                Pathogen = "Chikungunya virus",
                DataType = ifelse(precision == "PRECISE", "point", "polygon"),
                Year = substr(Year,7, 10),
                IfPolygon_AdminLevel = NA,
                Country = sapply(strsplit(LocationName, ", "), function(x, ...) x[ length(x)]),
                Source = "Nsoesie2016",
                NumCases = NA,
                CasesMetric = "Outbreak location",
                DiagnosticMetric = "Not specified",
                SourceType = "Compiled from literature and surveillance reports (IHME)")
chk1$IfPolygon_AdminLevel[ chk1$precision == "ADMIN1" ] = 1
chk1$IfPolygon_AdminLevel[ chk1$precision == "ADMIN2" ] = 2
chk1$IfPolygon_AdminLevel[ chk1$precision == "ADMIN3" ] = 3
chk1$IfPolygon_AdminLevel[ chk1$precision == "ADMIN4" ] = 4


# get polygons and correct IDs for admin1 and admin2 level reports

# GAUL polys (ADM1)
adm1_polys = chk1 %>% dplyr::filter(IfPolygon_AdminLevel == 1) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")

# missing with lat-lons - extract polygon directly
missing = adm1_polys %>% dplyr::filter(!ADMcode %in% g1$ADM1_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADMcode = g1_missing$ADM1_CODE

# combine
ad1 = rbind(
  adm1_polys %>% dplyr::filter(ADMcode %in% g1$ADM1_CODE),
  missing)
shp1 = g1 %>% dplyr::filter(ADM1_CODE %in% ad1$ADMcode)
print(paste("Found", nrow(shp1), "of", n_distinct(ad1$ADMcode), "GAUL1 polygons", sep=" "))


# GAUL polys (ADM2)
adm1_polys = chk1 %>% dplyr::filter(IfPolygon_AdminLevel == 2) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")

# missing with lat-lons - extract polygon directly
missing = adm1_polys %>% dplyr::filter(!ADMcode %in% g1$ADM2_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADMcode = g1_missing$ADM2_CODE

# combine
ad2 = rbind(
  adm1_polys %>% dplyr::filter(ADMcode %in% g1$ADM2_CODE),
  missing)
shp2 = g1 %>% dplyr::filter(ADM2_CODE %in% ad2$ADMcode)
print(paste("Found", nrow(shp2), "of", n_distinct(ad2$ADMcode), "GAUL2 polygons", sep=" "))



#  combine and harmonise polys and datasets 

# recombine chik dataset
ad1$ADMcode = paste("GAUL1", ad1$ADMcode, sep="_")
ad2$ADMcode = paste("GAUL2", ad2$ADMcode, sep="_")
chk1 = do.call(
  rbind.data.frame,
  list(chk1 %>% dplyr::filter(DataType == "point" | IfPolygon_AdminLevel %in% 3:4), ad1, ad2)
) %>%
  dplyr::select(-precision) %>%
  dplyr::mutate(Admin_name = NA, 
                ADMSource = "GAUL")

# combine polys
p1 = shp1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = shp2 %>% 
  dplyr::mutate(ADMcode = paste("GAUL2", ADM2_CODE, sep="_"),
                Admin_name = ADM2_NAME,
                Admin_level = 2) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
shp1 = rbind(p1, p2)
all(shp1$ADMcode %in% chk1$ADMcode)





# ------------ Brazilian surveillance data (2017 onward) -------------

# n.b. not using these in models 
chk2 = read.csv("./data/spillovers/Chikungunya/brazil_moh/chikungunya_monthlyTS_2017_2020.csv") %>%
  dplyr::filter(Disease == "Chikungunya" & Location_type == "Notification_municipality") %>%
  dplyr::group_by(IBGE6, Year) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE),
                   MunipName = head(MunipName, 1)) %>%
  dplyr::mutate(Disease = "Chikungunya", 
                Pathogen = "Chikungunya virus",
                DataType = "polygon",
                Country = "Brazil",
                ADMSource = "Brazil_shapefile",
                IfPolygon_AdminLevel = 2,
                Source = "Brazil DATASUS system") %>%
  dplyr::rename("Admin_name" = MunipName,
                "ADMcode" = IBGE6) %>%
  dplyr::mutate(LocationName = Admin_name) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Admin_name, LocationName, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

shp = sf::st_read("./data/spillovers/Chikungunya/brazil_moh/Brazil_shp_harm_2022.shp") %>%
  dplyr::select(IBGE6, name_mn) %>%
  dplyr::rename("ADMcode"=IBGE6, "Admin_name"=name_mn) %>%
  sf::st_transform(crs = 4326)
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) 

# combine with lat-lon info
chk2 = chk2 %>%
  dplyr::left_join(
    shp %>% st_drop_geometry() %>% dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% dplyr::mutate(ADMcode = as.numeric(ADMcode)) %>% dplyr::select(-Admin_name)
  ) %>%
  dplyr::arrange(Year, ADMcode)

Encoding(chk2$Admin_name) = "latin1"
chk2$CasesMetric = "All case counts (municipality of notification)"
chk2$DiagnosticMetric = "Not specified"
chk2$SourceType = "National case surveillance system"



# ------------ combine datasets and save ------------

# combine both datasets
chk = do.call(
  rbind.data.frame, 
  list(chk1, chk2)
)

# combine shapefiles
shp = rbind(
  shp1, 
  shp %>% dplyr::select(-Longitude, -Latitude) %>% dplyr::mutate(Admin_level = 2) %>% dplyr::filter(ADMcode %in% chk$ADMcode)
)

# save data and shapefile
data.table::fwrite(chk, "./output/spillovers_processed/spillovers_chikv.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_chikv_shp.R")















# save data and shapefile
write.csv(den, "./output/spillovers_processed/spillovers_dengue.csv", row.names=FALSE)
#save(shp, file="./output/spillovers_processed/spillovers_dengue_shp.R")


# viz
library(maptools)
data("wrld_simpl")
ws = st_as_sf(wrld_simpl)
ggplot() + 
  geom_sf(data=ws, fill="grey80", color=NA) + 
  maptheme + 
  geom_point(data=cchf, aes(x=Longitude, y=Latitude), col="black", size=0.25)

