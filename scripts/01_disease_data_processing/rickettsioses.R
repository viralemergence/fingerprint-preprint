

# ================ Spotted fever group rickettsioses (Brazil) ======================

# Brazilian surveillance data
# https://www.frontiersin.org/articles/10.3389/fcimb.2013.00027/full

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")



# ----------- 1. Surveillance reports, Brazil, 2006 to 2017 --------------------

# disease reports from municipality of residence (not notification)
rick = read.csv("./data/spillovers/Rickettsioses/brazil_moh/tickfever_rickettsia_febremaculosa_monthlyTS_2001_2020.csv") %>%
  dplyr::filter(Disease == "Tickfever_Rickettsia_FebreMaculosa" & Cases_metric == "Confirmed") %>%
  dplyr::filter(Location_type == "Residence_municipality") %>%
  dplyr::group_by(IBGE6, Year) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE),
                   MunipName = head(MunipName, 1)) %>%
  dplyr::mutate(Disease = "Brazilian spotted fever", 
                Pathogen = "Rickettsia rickettsii / Rickettsia parkerii",
                DataType = "polygon",
                Country = "Brazil",
                ADMSource = "Brazil_shapefile",
                IfPolygon_AdminLevel = 2,
                Source = "Brazil DATASUS system") %>%
  dplyr::rename("Admin_name" = MunipName,
                "ADMcode" = IBGE6) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Admin_name, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

shp = sf::st_read("./data/spillovers/Rickettsioses/brazil_moh/Brazil_shp_harm_2022.shp") %>%
  dplyr::select(IBGE6, name_mn) %>%
  dplyr::rename("ADMcode"=IBGE6, "Admin_name"=name_mn) %>%
  sf::st_transform(crs = 4326)
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) 

# combine with lat-lon info
rick = rick %>%
  dplyr::left_join(
    shp %>% st_drop_geometry() %>% dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% dplyr::mutate(ADMcode = as.numeric(ADMcode)) %>% dplyr::select(-Admin_name)
  ) %>%
  dplyr::arrange(Year, ADMcode)

# additional metadata
Encoding(rick$Admin_name) = "latin1"
rick$CasesMetric = "Confirmed case counts (municipality of residence)"
rick$DiagnosticMetric = "Not specified"
rick$SourceType = "National case surveillance system"



# ------------ combine datasets and save -----------------

# save
data.table::fwrite(rick, "./output/spillovers_processed/spillovers_rickettsia.csv", row.names=FALSE)

# combine shapefiles
shp = shp %>%
  dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% 
  dplyr::mutate(ADMcode = as.numeric(ADMcode), Admin_level=2) %>%
  dplyr::filter(ADMcode %in% rick$ADMcode)
save(shp, file="./output/spillovers_processed/spillovers_rickettsia_shp.R")



# viz
# library(maptools)
# data("wrld_simpl")
# ws = sf::st_as_sf(wrld_simpl)
# ggplot() +
#   geom_sf(data=ws, fill="grey80", color=NA) +
#   maptheme +
#   geom_point(data=rick, aes(x=Longitude, y=Latitude), col="black", size=0.25)
