
# ====================== Lassa fever ===========================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")



# ------------ 1. spillovers database from Gibb et al 2017 --------------

las1 = read.csv("./data/spillovers/Lassa/gibb2017/Lassafever-outbreaks-seroprev-DATABASE_updatedJul2019.csv") %>%
  dplyr::mutate(Year_start = as.numeric(Year_start)) %>%
  dplyr::filter(Year_start >= 1970) %>%
  dplyr::filter(!is.na(Latitude)) %>%
  dplyr::filter(Rodent_or_human == "Human") %>%
  #dplyr::filter(Acutedisease_or_seroprevalence == "Acute") %>%
  dplyr::mutate(Number_acutecases = replace(Number_acutecases, Acutedisease_or_seroprevalence == "Seroprevalence", "Sero")) %>%
  dplyr::filter(Number_acutecases != "0") %>%
  dplyr::select(Country, Admin_region, Town, Latitude, Longitude, Year_start, Number_acutecases) %>%
  dplyr::mutate(Disease = "Lassa fever", Pathogen = "Lassa arenavirus", DataType = "point", IfPolygon_AdminLevel=NA, ADMcode=NA, ADMSource=NA, Source="Gibb2017") %>%
  dplyr::rename("Admin_name"=Admin_region,
                "Year"=Year_start,
                "NumCases"=Number_acutecases,
                "LocationName"=Town) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, NumCases, LocationName, Admin_name, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

# fix num cases
las1$NumCases[ las1$NumCases == ">6" ] = "6"
las1 = las1 %>%
  dplyr::mutate(NumCases = as.numeric(NumCases)) %>%
  dplyr::filter(!is.na(NumCases))



# -------------- 2. Nigeria surveillance data 2012-2019 ----------------

# incidence database from Nigeria 2012-2019, from Redding et al 2021 (supp data)
las2 <- read.csv("./data/spillovers/Lassa/redding2021/Lassa_LGA_2012-2019_WER-SitRep_Mar2021.csv") %>%
  dplyr::filter(variable == "Confirmed" & LGA_known == 1) %>%
  dplyr::select(LGA, Year, Source, admin1Name, value) 

# filter to data sources following Redding 2021
las2 = las2[ which(las2$Source == "WER" & !las2$Year %in% 2017:2018 | las2$Source == "SitRep_1" & las2$Year == 2017 | las2$Source == "SitRep_2" & !is.na(las2$value)),  ]
las2$Year[ las2$Year == 2020 ] = 2019

# combine, summarise and filter to presence-only
las2 = las2 %>%
  dplyr::group_by(LGA, Year) %>%
  dplyr::summarise(LocationName = paste(head(admin1Name, 1), "state", sep=" "),
                   NumCases = sum(value, na.rm=TRUE)) %>%
  dplyr::filter(NumCases > 0)

# read shapefile and get latlons
shp = sf::st_read("./data/spillovers/Lassa/redding2021/LF_nigeria_shapefile.shp")
locs = shp %>%
  cbind(
    st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
  ) %>%
  dplyr::select(LGA, Longitude, Latitude, admin2Pcod) %>%
  st_drop_geometry()
las2 = dplyr::left_join(las2, locs)

# process and rename accordingly
las2 = las2 %>%
  dplyr::mutate(Disease = "Lassa fever", Pathogen = "Lassa arenavirus", DataType = "polygon", Country="Nigeria", IfPolygon_AdminLevel=2, ADMSource="Redding2021_shapefile", Source="Redding2021") %>%
  dplyr::rename("Admin_name"=LGA, "ADMcode"=admin2Pcod) %>%
  dplyr::select(names(las1))


# ------------ 3. Liberia surveillance data 2019-20 ----------------

las3 = read.csv("./data/spillovers/Lassa/liberia2022/lassasurveillance_liberia_2022.csv") %>%
  dplyr::mutate(Pathogen = "Lassa arenavirus", DataType = "polygon",
                IfPolygon_AdminLevel=2, ADMSource="GADM") %>%
  dplyr::rename("Admin_name"=Admin2) %>%
  dplyr::mutate(LocationName = paste(Admin_name, Admin1, sep=", ")) %>%
  dplyr::select(-Admin1, -URL, -CasesMetric, -DiagnosticMethod, -Note)

# get centroids and adm codes
shp2 = sf::st_read("./data/spillovers/Lassa/liberia2022/gadm41_LBR_2.shp") %>%
  dplyr::filter(NAME_2 %in% las3$Admin_name) %>%
  dplyr::select(GID_2, NAME_2) %>%
  dplyr::rename("ADMcode"=1, "Admin_name"=2)
shp2 = shp2 %>%
  cbind(
    st_coordinates(st_centroid(shp2)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
  )
shp2$ADMarea = as.vector(sf::st_area(shp2)/10^6)

# add latlons and adm codes to disease data
las3 = left_join(las3, shp2 %>% dplyr::select(Admin_name, ADMcode, Longitude, Latitude) %>% st_drop_geometry())



# ----------- combine and save -------------

# combine spillovers
las = do.call(
  rbind.data.frame,
  list(las1, las2, las3)
)

# edit shapefile
shp = shp %>%
  left_join(dplyr::select(locs, LGA, Longitude, Latitude)) %>%
  dplyr::select(admin2Pcod, LGA, Longitude, Latitude) %>%
  dplyr::rename("ADMcode"=1, "Admin_name"=LGA) %>%
  dplyr::filter(ADMcode %in% las$ADMcode)
shp$ADMarea = as.vector(sf::st_area(shp)/10^6)
shp = rbind(shp, shp2)

# save data and shapefile
write.csv(las, "./output/spillovers_processed/spillovers_lassa.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_lassa_shp.R")


# viz
ws = rnaturalearth::ne_coastline(returnclass = "sf")
ws = st_crop(ws, extent(c(-18, 18, 0, 18)))
ggplot() + 
  geom_sf(data=ws, fill="grey80", color=NA) + 
  maptheme + 
  geom_point(data=las, aes(x=Longitude, y=Latitude), col="black", size=0.1)

