

# ================= Oropouche fever =================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")

# shapefile for red
shp_ref = sf::st_read("./data/shapefiles/africa_shp.shp")



# ------------ point data with uncertainty bounds -------------------

# data
oro = read.csv("./data/spillovers/Oropouche/v2_DRA/full_occs2.csv") %>%
  dplyr::filter(Reservoir == "humans") %>%
  dplyr::select(Longitude, Latitude, Country, Province.State, Locality, Years.of.detection, Uncertainty) %>%
  dplyr::rename("Admin_name" = 4, "LocationName" = 5, "Year" = 6, "buffer_radius" = 7) %>%
  dplyr::mutate(
    Year = replace(Year, Year == "1960-1981", 1981),
    Year = replace(Year, Year == "1960-1961_1968-69_1979-80_2008", 2008),
    Year = replace(Year, Year == "1961_1968-69_1979-80_2008", 2008),
    Year = replace(Year, Year == "1967_1979-80", 1980),
    Year = replace(Year, Year == "1974-75", 1975),
    Year = replace(Year, Year == "1980, 2006", 2006),
    Year = replace(Year, Year == "1980_2006", 2006),
    Year = replace(Year, Year == "1980_2009 ", 2009),
    Year = replace(Year, Year == "1981_2007-2008", 2008),
    Year = replace(Year, Year == "1992_1998_2000-2007", 2007),
    Year = replace(Year, Year == "1994_1998_2000-2007_2016", 2016),
    Year = replace(Year, Year == "2011-2012", 2012),
    Year = replace(Year, Year == "2019-2020", 2020)
  ) %>%
  distinct() %>%
  dplyr::mutate(Disease = "Oropouche fever",
                Pathogen = "Oropouche virus",
                Source = "RomeroAlvarez",
                DataType = "polygon",
                IfPolygon_AdminLevel=NA,
                ADMcode = NA, 
                ADMSource=NA,
                NumCases = NA) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Latitude, Longitude, Admin_name, LocationName, NumCases, IfPolygon_AdminLevel,
                ADMcode, ADMSource, Source, buffer_radius)

# IF SEPARATING BY YEAR - BEST NOT TO
# oro = read.csv("./data/spillovers/Oropouche/v2_DRA/full_occs2.csv") %>%
#   dplyr::filter(Reservoir == "humans") %>%
#   dplyr::select(Longitude, Latitude, Country, Province.State, Locality, Years.of.detection, Uncertainty) %>%
#   dplyr::rename("Admin_name" = 4, "LocationName" = 5, "Year" = 6, "buffer_radius" = 7) %>%
#   dplyr::mutate(
#     Year = replace(Year, Year == "1960-1981", 1981),
#     Year = replace(Year, Year == "1960-1961_1968-69_1979-80_2008", "1961_1969_1980_2008"),
#     Year = replace(Year, Year == "1961_1968-69_1979-80_2008", "1961_1969_1980_2008"),
#     Year = replace(Year, Year == "1967_1979-80", "1967_1980"),
#     Year = replace(Year, Year == "1974-75", 1975),
#     Year = replace(Year, Year == "1980, 2006", "1980_2006"),
#     Year = replace(Year, Year == "1981_2007-2008", "1981_2008"),
#     Year = replace(Year, Year == "1992_1998_2000-2007", "1992_1998_2007"),
#     Year = replace(Year, Year == "1994_1998_2000-2007_2016", "1994_1998_2007_2016"),
#     Year = replace(Year, Year == "2011-2012", 2012),
#     Year = replace(Year, Year == "2019-2020", 2020)
#   ) %>%
#   distinct() %>%
#   dplyr::mutate(Disease = "Oropouche fever",
#                 Pathogen = "Oropouche virus",
#                 Source = "RomeroAlvarez",
#                 DataType = "point",
#                 IfPolygon_AdminLevel=NA,
#                 ADMcode = NA, 
#                 ADMSource=NA,
#                 NumCases = NA) %>%
#   dplyr::select(Disease, Pathogen, DataType, Country, Year, Latitude, Longitude, Admin_name, LocationName, NumCases, IfPolygon_AdminLevel,
#                 ADMcode, ADMSource, Source, buffer_radius) %>%
#   tidyr::separate_rows(Year, sep="_") %>%
#   dplyr::mutate(Year = as.numeric(substr(Year, 1, 4)))


# ----------------- generate buffers for points -----------------

oro$ADMcode = paste("ORObuffer", 1:nrow(oro), sep="_")
for_buf = oro %>% dplyr::select(ADMcode, Longitude, Latitude, buffer_radius)
coordinates(for_buf) = ~Longitude+Latitude
for_buf = st_as_sf(for_buf) 
st_crs(for_buf) = crs(shp_ref)
for_buf$radius_m = for_buf$buffer_radius
createBuffer = function(x){
  return(sf::st_buffer( for_buf[x, ], dist=for_buf$radius_m[x] ))
}
buf = data.frame()
for(i in 1:nrow(for_buf)){
  bx = createBuffer(i)
  buf = rbind(buf, bx)
}
buf = buf %>% dplyr::select(-buffer_radius)
shp = buf

# -------------- save ---------------

write.csv(oro, "./output/spillovers_processed/spillovers_oropouche.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_oropouche_shp.R")

