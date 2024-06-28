

# ==================== CCHF =======================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")


# ======================= Data wrangling ====================

# ------------- 1. CCHF Messina global compendium: yearly occurrences of CCHF at point and polygon level globally
# match to GAUL database for polygon boundaries

cchf1 = read.csv("./data/spillovers/CCHF/messina_globalcompendium/CCHF_1953_2012_Messina.csv") %>%
  dplyr::select(2, 3, 4, 5, 7, 8, 9, 10) %>%
  dplyr::rename(
    "DataType" = LOCATION_TYPE,
    "IfPolygon_AdminLevel" = ADMIN_LEVEL,
    "Year"=YEAR,
    "Latitude"=LATITUDE,
    "Longitude"=LONGITUDE,
    "Country"=COUNTRY,
  ) %>%
  dplyr::mutate(ADMcode = GAUL_AD2,
                Disease = "Crimean-Congo haemorrhagic fever",
                ADMSource = "GAUL", 
                NumCases = NA, 
                LocationName = NA,
                Pathogen = "Crimean-Congo haemorrhagic fever virus")

# fix ADMcodes; none for points, correct level for polygons
cchf1$ADMcode[ cchf1$IfPolygon_AdminLevel == 1 ] = cchf1$GAUL_AD1[ cchf1$IfPolygon_AdminLevel == 1 ] 
cchf1$ADMcode[ cchf1$DataType == "point" ] = NA
cchf1 = cchf1 %>%
  dplyr::select(-GAUL_AD1, -GAUL_AD2) %>%
  dplyr::mutate(Source = "Messina2015",
                IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, IfPolygon_AdminLevel == -999, NA))

# get admin1 polygons
adm1_polys = cchf1 %>% dplyr::filter(DataType  == "polygon" & IfPolygon_AdminLevel == 1) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")
gps = unique(adm1_polys$ADMcode)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))

# adm2 polys from GAUL
adm2_polys = cchf1 %>% dplyr::filter(DataType  == "polygon" & IfPolygon_AdminLevel == 2) 
g2 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")
gps = unique(adm2_polys$ADMcode)
gaul2 = g2 %>% dplyr::filter(ADM2_CODE %in% as.numeric(gps))
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))

# use lat-lons to look up missing polys
missing = adm2_polys[ !adm2_polys$ADMcode %in% gaul2$ADM2_CODE, ]
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g2)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g2)
g2_missing = g2[ as.numeric(ii), ]
missing$ADMcode = g2_missing$ADM2_CODE

# combine together
adm2_polys = 
  rbind(
    adm2_polys[ adm2_polys$ADMcode %in% gaul2$ADM2_CODE, ],
    missing
  )
gaul2 = rbind(
  gaul2,
  g2_missing
)
all(gaul2$ADM2_CODE %in% adm2_polys$ADMcode)

# create combined cchf data
cchf1 = rbind(
  cchf1[ cchf1$DataType == "point", ],
  adm1_polys %>% dplyr::mutate(ADMcode = paste("GAUL1", ADMcode, sep="_")),
  adm2_polys %>% dplyr::mutate(ADMcode = paste("GAUL2", ADMcode, sep="_"))
)

# create combined shapefile
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



# ---------- 2. Ak et al, yearly occurrences at province-level in Turkey

cchf2 = read.csv("./data/spillovers/CCHF/ak_cchfsurveillance_turkey/Ak2020_cchf_cases.csv") %>%
  dplyr::mutate(Disease = "Crimean-Congo haemorrhagic fever", 
                Pathogen = "Crimean-Congo haemorrhagic fever virus",
                DataType = "polygon",
                Country = "Turkey", 
                IfPolygon_AdminLevel = 1,
                ADMSource = "GADM", 
                Source = "Ak2020") %>%
  dplyr::rename("ADMcode"=Code_GADM,
                "NumCases"=Cases, 
                "LocationName"=Name_GADM) %>%
  dplyr::select(-Province)

# get lat-lon coordinates
shp2 = sf::st_read("./data/spillovers/CCHF/ak_cchfsurveillance_turkey/shapefile/gadm40_TUR_1.shp")
coords = cbind(
  shp2, 
  st_coordinates(st_centroid(shp2)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
  ) %>%
  dplyr::select(ID_1, Longitude, Latitude) %>%
  st_drop_geometry()
cchf2 = dplyr::left_join(cchf2, coords, by=c("ADMcode"="ID_1"))

# align names with shp1
shp2 = shp2 %>%
  dplyr::mutate(ADMcode = ID_1,
                Admin_name = NAME_1,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level)



# ------------ 3. Lule et al - African CCHF occurrences 2012 to present

cchf3 = read.csv("./data/spillovers/CCHF/gibb_uganda/RoryGibb_CCHFAfrica_2010_2020.csv") %>%
  dplyr::filter(Type == "Human" & Sero==FALSE) %>%
  dplyr::select(Year, Country, Locale, Longitude, Latitude, AcuteCases) %>%
  dplyr::rename("LocationName"=Locale, 
                "NumCases"=AcuteCases) %>%
  dplyr::mutate(Source = "Gibb2021", ADMcode=NA, DataType="point", IfPolygon_AdminLevel=NA, ADMSource=NA,
                Disease = "Crimean-Congo haemorrhagic fever", 
                Pathogen = "Crimean-Congo haemorrhagic fever virus") %>%
  dplyr::mutate(Year = replace(Year, Year == "2015-16", "2015"))
  
cchf3.1 = read.csv("./data/spillovers/CCHF/lule_africa/cchf_location_datasources_georef.csv", stringsAsFactors = FALSE) %>%
  dplyr::filter(Country != "Uganda") %>%
  dplyr::filter(Sero == FALSE | Sero_Type == "Diagnostic") %>%
  dplyr::filter(Type == "Human") %>%
  dplyr::filter(!is.na(Longitude)) %>%
  dplyr::select(Year, Country, Locale, Longitude, Latitude, AcuteCases) %>%
  dplyr::mutate(Year = replace(Year, Year == "2015-16", "2015")) %>%
  dplyr::rename("LocationName"=Locale, 
                "NumCases"=AcuteCases)
cchf3.2 = read.csv("./data/spillovers/CCHF/lule_africa/balinandietal_2022.csv") %>%
  dplyr::select(-Source) %>%
  dplyr::rename("Year"=Years) %>%
  dplyr::mutate(LocationName = "", NumCases = 1, Year = 2015)  
cchf3 = rbind(cchf3.1, cchf3.2) %>%
  dplyr::mutate(Source = "Lule2022", ADMcode=NA, DataType="point", IfPolygon_AdminLevel=NA, ADMSource=NA,
                Disease = "Crimean-Congo haemorrhagic fever", 
                Pathogen = "Crimean-Congo haemorrhagic fever virus") 


# combine and save
cchf = do.call(
  rbind.data.frame,
  list(cchf1, cchf2, cchf3)
)

shp = rbind(
  shp1, shp2
)

# save data and shapefile
write.csv(cchf, "./output/spillovers_processed/spillovers_cchf.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_cchf_shp.R")


# viz
library(maptools)
data("wrld_simpl")
ws = st_as_sf(wrld_simpl)
ggplot() + 
  geom_sf(data=ws, fill="grey80", color=NA) + 
  maptheme + 
  geom_point(data=cchf, aes(x=Longitude, y=Latitude), col="black", size=0.25)

