
# ------------- Junin virus spillovers ----------

# from Argentinian Junin virus surveillance system, digitised from annual reports
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)
source("./scripts/00_plot_themes.R")



# -------------- pre-prep: cross-reference to locality lat-lons -------------------

# # create lookup table
# jun = read.csv("./data/spillovers/Junin/junin_extracted.csv") %>%
#   dplyr::select(Province, Localidad) %>%
#   distinct()
# 
# # geojson of localities
# gg = rgdal::readOGR("./data/spillovers/Junin/localidades.gjson.geojson")
# gg$province = substr(gg$provincia, 14, 50)
# gg = gg[ grep("Buenos|Córdoba|Santa Fe|La Pampa", gg$province), ]
# gg = as.data.frame(gg)
# gg$provincename = NA
# gg$provincename[ grep("Buenos", gg$province) ] = "Buenos Aires"
# gg$provincename[ grep("doba", gg$province) ] = "Cordoba"
# gg$provincename[ grep("Santa", gg$province) ] = "Santa Fe"
# gg$provincename[ grep("Pampa", gg$province) ] = "La Pampa"
# 
# gg = gg %>% dplyr::select(provincename, nombre, coords.x1, coords.x2, departamento) %>%
#   dplyr::arrange(provincename, nombre) %>%
#   dplyr::mutate(ID = paste("arg", 1:length(nombre), sep="_"))
# write.csv(gg, "./data/spillovers/Junin/latlon_lookup_codes.csv", row.names=FALSE)

# read back in junin data and latlons and combine
# jun = read.csv("./data/spillovers/Junin/data_extraction/junin_extracted.csv")
# ll = read.csv("./data/spillovers/Junin/data_extraction/junin_lookup_table.csv") %>% dplyr::select(-ID)
# jun = left_join(jun, ll) %>%
#   dplyr::filter(!is.na(Longitude)) # only 1 record can't be georef'd
# 
# # save
# write.csv(jun, "./data/spillovers/Junin/junin_georeferenced.csv", row.names=FALSE)



# ----------------------------------------------------------------------------------------

# read point data and add additional cols
jun = read.csv("./data/spillovers/Junin/junin_georeferenced.csv") %>% 
  dplyr::rename("LocationName"=Localidad,
                "Admin_name"=Province,
                "NumCases"=Cases_NotificationsMinusNegatives) %>%
  dplyr::mutate(DataType = "point",
                Country = "Argentina",
                IfPolygon_AdminLevel = NA,
                ADMcode = NA,
                ADMSource = NA,
                DiagnosticMetric = "Notifications minus negatives (considered the most accurate indicator of true cases)",
                CasesMetric = "Confirmed/probable case counts",
                SourceType = "National case surveillance system") %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Latitude, Longitude, Admin_name, LocationName, 
                NumCases, DiagnosticMetric, CasesMetric, IfPolygon_AdminLevel,
                ADMcode, ADMSource, Source, SourceType) %>%
  dplyr::arrange(Year, LocationName)

write.csv(jun, "./output/spillovers_processed/spillovers_junin.csv", row.names=FALSE)


# some quick viz

# ar = sf::st_read("C:/Users/roryj/Downloads/gadm41_ARG_shp/gadm41_ARG_1.shp") %>%
#   dplyr::filter(NAME_1 %in% c("Buenos Aires", "Córdoba", "La Pampa", "Santa Fe"))
# ggplot() + 
#   geom_sf(data=ar, col="grey50", fill="grey95") +
#   maptheme + 
#   geom_point(data=jun, aes(Longitude, Latitude, size=NumCases), alpha=0.4, col="coral2") +
#   facet_wrap(~Year, nrow=2)

# jun %>%
#   dplyr::group_by(Admin_name, Year) %>%
#   dplyr::summarise(NumCases = sum(NumCases)) %>%
#   ggplot() + 
#   geom_line(aes(Year, NumCases)) + 
#   facet_wrap(~Admin_name)


jun = read.csv("./output/spillovers_processed/spillovers_junin.csv")
ar = sf::st_read("C:/Users/roryj/Downloads/gadm41_ARG_shp/gadm41_ARG_1.shp") %>%
  dplyr::filter(NAME_1 %in% c("Buenos Aires", "Córdoba", "La Pampa", "Santa Fe"))
pp = ggplot() +
  geom_sf(data=ar, col="grey50", fill="grey95") +
  maptheme +
  geom_point(data=jun, aes(Longitude, Latitude, size=NumCases), alpha=0.4, col="coral2") +
  facet_wrap(~Year, nrow=3)
ggsave(pp, file="./junin_spatiotemp.png", device="png", units="in", width=10, height=8, dpi=600)
