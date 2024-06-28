


# ================ Plague ======================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")


# ----------------- 1. Western USA (1950-2005, Carlson 2021) ---------------

# counties shapefile with cases attached, get centroid
shp1 = sf::st_read("./data/spillovers/Plague/carlson2021_usa/CountiesNoHole.shp") %>%
  dplyr::mutate(ADMcode = paste(STATEFP_1, COUNTYFP_1, sep="")) %>%
  sf::st_transform(crs = 4326)
shp1 = cbind(
  shp1, 
  st_coordinates(st_centroid(shp1)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) 

# split out cases and process
pl1 = shp1 %>% st_drop_geometry() 
pl1 = cbind(pl1[ , c("ADMcode", "Longitude", "Latitude", "NAME_1") ], pl1[ , grep("YEAR", names(pl1)) ]) %>%
  reshape2::melt(id.vars = 1:4) %>%
  dplyr::rename("NumCases"=value,
                "Year"=variable,
                "Admin_name"=NAME_1) %>%
  dplyr::mutate(Year = substr(Year, 6, 11),
                DataType = "polygon", 
                Disease = "Plague", 
                Pathogen = "Yersinia pestis",
                Country = "USA", 
                ADMSource = "US counties shapefile",
                IfPolygon_AdminLevel = 2,
                Source = "US CDC",
                CasesMetric = "Confirmed cases (county of infection)",
                DiagnosticMetric = "Not specified",
                SourceType = "National case surveillance system") %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, Admin_name, NumCases, CasesMetric, DiagnosticMetric, IfPolygon_AdminLevel, ADMcode, ADMSource, Source, SourceType)

# shapefile
shp1 = shp1 %>%
  dplyr::select(ADMcode, NAME_1) %>% dplyr::rename("Admin_name"=NAME_1) %>% distinct()

# for now only save >0 cases (but in full fingerprint DB, keep zeroes)
pl1 = pl1 %>% dplyr::filter(NumCases > 0)





# -------------- 2. Plague global (Boris Schmid) ----------------

pl2 = read.csv("./data/spillovers/Plague/schmid_global/plague_verena.csv") %>%
  dplyr::mutate(Year = substr(sampling.date, 1, 4)) %>%
  dplyr::filter(!is.na(longitude)) %>%
  dplyr::mutate(ID = 1:length(Year))

# remove GBIF and Carlson (Western USA) records to avoid overlaps with above
pl2 = pl2[ -grep("Carlson", pl2$source), ]
pl2 = pl2[ -grep("GBIF", pl2$source), ]

# overlay with global countries shapefile to extract country info
cc = sf::st_read("./data/shapefiles/world-administrative-boundaries.shp")
pl_sf = sf::st_as_sf(pl2, coords = c("longitude", "latitude"))
st_crs(pl_sf) = crs(cc)
intersection = as.data.frame(sf::st_intersects(pl_sf, cc))
intersection$Country = cc$name[ intersection$col.id ]
intersection = dplyr::rename(intersection, "ID"=row.id)
intersection = intersection %>% dplyr::select(-col.id)

# combine and NAs are Madagascar coastal points
pl2 = left_join(pl2, intersection)
pl2$Country[ is.na(pl2$Country) ] = "Madagascar"

# pl2 %>%
#   dplyr::filter(Country %in% c("Madagascar")) %>%
#   ggplot() + 
#   geom_sf(data=cc[ cc$name == "Madagascar", ], fill="grey90", color="grey80") +
#   geom_point(aes(longitude, latitude, color=Country)) 

# standardise column names
pl2 = pl2 %>%
  dplyr::rename(
    "Longitude"=longitude,
    "Latitude"=latitude,
    "Source"=source
  ) %>%
  dplyr::mutate(
    NumCases = 1,
    Admin_name = NA,
    ADMcode=NA,
    DataType = "point",
    Disease = "Plague", 
    Pathogen = "Yersinia pestis",
    ADMSource = NA,
    IfPolygon_AdminLevel = 2,
    CasesMetric = "Confirmed cases",
    DiagnosticMetric = "Sequencing (NCBI)",
    SourceType = "Boris Schmid global data"
  ) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, Admin_name, NumCases, CasesMetric, DiagnosticMetric, IfPolygon_AdminLevel, ADMcode, ADMSource, Source, SourceType)



# combine and save full dataset
pl1 = rbind(pl1, pl2)
data.table::fwrite(pl1, "./output/spillovers_processed/spillovers_plague.csv", row.names=FALSE)
save(shp1, file="./output/spillovers_processed/spillovers_plague_shp.R")
