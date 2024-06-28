

# ================== Lyme disease ====================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")



# -------------------- 20 years of surveillance data, USA -------------------------

# occurrences from USA
ld = read.csv("./data/spillovers/Lyme/us_cdc/LD-Case-Counts-by-County-00-19.csv") %>%
  reshape2::melt(id.vars = 1:4) %>%
  dplyr::rename("Year"=variable, "NumCases"=value) %>%
  dplyr::mutate(Year = as.numeric(substr(Year, 6, 10))) %>%
  dplyr::filter(CTYCODE != 999) %>%
  dplyr::filter(NumCases > 0) 
ld = ld[ -which(ld$Ctyname == "Bedford city" & ld$Stname == "Virginia"), ]

# US shapefile from census bureau with matching codes, only one issue
shp = sf::st_read("./data/spillovers/Lyme/us_cdc/cb_2018_us_county_500k.shp") %>%
  dplyr::mutate(CTYCODE = as.integer(COUNTYFP),
                STCODE = as.integer(STATEFP)) 
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
)

# check match on county codes - all correct
# matched = ld %>% dplyr::select(1:4) %>%
#   distinct() %>%
#   left_join(
#     shp %>%
#       dplyr::select(CTYCODE, STCODE, NAME) %>% 
#       st_drop_geometry()
#   )

# create standardised ADMcode
ld$ADMcode = paste(ld$STCODE, ld$CTYCODE, sep="_")
shp$ADMcode = paste(shp$STCODE, shp$CTYCODE, sep="_")

# subset shp to required fields
shp = shp %>% 
  dplyr::select(ADMcode, NAME, Longitude, Latitude) %>%
  dplyr::rename("Admin_name"=NAME)

# combine info with ld
ld = ld %>%
  dplyr::left_join(
    shp %>% st_drop_geometry()
  ) %>%
  dplyr::select(-Ctyname, -Stname, -STCODE, -CTYCODE) %>%
  dplyr::mutate(Disease = "Lyme disease",
                Pathogen = "Borrelia burgdorferi",
                ADMSource = "USCB_shapefile",
                DataType = "polygon",
                Country = "USA",
                IfPolygon_AdminLevel = 2,
                Source = "CDC") %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, 
              Longitude, Latitude, Admin_name, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)


# save
write.csv(ld, "./output/spillovers_processed/spillovers_lyme.csv", row.names = FALSE)
save(shp, file="./output/spillovers_processed/spillovers_lyme_shp.R")

