

# ==================== Dengue from global IHME compendium =======================

# IHME data up to 2015
# https://ghdx.healthdata.org/record/ihme-data/global-geo-referenced-dengue-occurrence-database-1960-2015

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")


# ======================= Data wrangling ====================

# 1. Points
den1 = read.csv("./data/spillovers/Dengue/ihme/IHME_DENGUE_1960_2015_POINTS_STANDARD_CHECKED_EC_2019M06D10.CSV") %>% 
  dplyr::select(Year, Longitude, Latitude) %>%
  dplyr::mutate(Disease = "Dengue", 
                Pathogen = "Dengue virus",
                DataType = "point", 
                IfPolygon_AdminLevel=NA, 
                Country = "",
                ADMcode = NA, 
                Source = "IHME",
                CasesMetric = "Outbreak_locations",
                DiagnosticMetric = "Not specified",
                SourceType = "Compiled from literature (IHME)") %>%
  dplyr::filter(Year >= 1985)



# 2. Admin1 polygons
den2 = read.csv("./data/spillovers/Dengue/ihme/IHME_DENGUE_1960_2015_POLY_STANDARD_CHECKED_EC_2019M06D10.CSV") %>% 
  dplyr::filter(Admin == 1) %>%
  dplyr::select(Year, Longitude, Latitude) %>%
  dplyr::mutate(Disease = "Dengue", 
                Pathogen = "Dengue virus",
                DataType = "polygon", 
                IfPolygon_AdminLevel=1, 
                Country = "",
                Source = "IHME",
                CasesMetric = "Outbreak_locations",
                DiagnosticMetric = "Not specified",
                SourceType = "Compiled from literature (IHME)") %>%
  dplyr::filter(Year >= 1985)

# extract adm1 polys from gadm
gadm = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level1")
to_ext = den2; coordinates(to_ext) = ~Longitude+Latitude
to_ext = st_as_sf(to_ext)
st_crs(to_ext) = st_crs(gadm)
sf::sf_use_s2(FALSE)
ii = st_intersects(to_ext, gadm)
gadm1 = gadm[ as.numeric(ii), ]
den2$ADMcode = gadm1$ID_1



# 3. Admin2 polygons
den3 = read.csv("./data/spillovers/Dengue/ihme/IHME_DENGUE_1960_2015_POLY_STANDARD_CHECKED_EC_2019M06D10.CSV") %>% 
  dplyr::filter(Admin == 2) %>%
  dplyr::select(Year, Longitude, Latitude) %>%
  dplyr::mutate(Disease = "Dengue", 
                Pathogen = "Dengue virus",
                DataType = "polygon", 
                IfPolygon_AdminLevel=2, 
                Country = "",
                Source = "IHME",
                CasesMetric = "Outbreak_locations",
                DiagnosticMetric = "Not specified",
                SourceType = "Compiled from literature (IHME)") %>%
  dplyr::filter(Year >= 1985)

# extract adm2 polys from gadm
gadm = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/gadm404-levels/gadm404-levels.gpkg", layer="level2")
to_ext = den3; coordinates(to_ext) = ~Longitude+Latitude
to_ext = st_as_sf(to_ext)
st_crs(to_ext) = st_crs(gadm)
sf::sf_use_s2(FALSE)
ii = st_intersects(to_ext, gadm)
gadm2 = gadm[ as.numeric(ii), ]
den3$ADMcode = gadm2$ID_2



# ------------ combine point and polygon data and save -------------

# combine and save
den = do.call(
  rbind.data.frame,
  list(den1, den2, den3)
)

# combine shapefiles
shp1 = gadm1 %>%
  dplyr::select(ID_1, NAME_1, COUNTRY, geom) %>%
  dplyr::rename("ADMcode"=1, "Admin_name"=2, "Country"=3, "geometry"=geom) %>%
  dplyr::mutate(Admin_level = 1)
shp2 = gadm2 %>%
  dplyr::select(ID_2, NAME_2, COUNTRY, geom) %>%
  dplyr::rename("ADMcode"=1, "Admin_name"=2, "Country"=3, "geometry"=geom) %>%
  dplyr::mutate(Admin_level = 2)
shp = rbind(shp1, shp2)
row.names(shp) = c()

# save data and shapefile
write.csv(den, "./output/spillovers_processed/spillovers_dengue.csv", row.names=FALSE)

# save shapefile as two halves (size issue)
shp1 = shp[ 1:2000, ]
shp2 = shp[ 2001:nrow(shp), ]
save(shp1, file="./output/spillovers_processed/spillovers_dengue_shp1.R")
save(shp2, file="./output/spillovers_processed/spillovers_dengue_shp2.R")



# # viz
# library(maptools)
# data("wrld_simpl")
# ws = st_as_sf(wrld_simpl)
# ggplot() + 
#   geom_sf(data=ws, fill="grey80", color=NA) + 
#   maptheme + 
#   geom_point(data=cchf, aes(x=Longitude, y=Latitude), col="black", size=0.25)

