

# ================= Anthrax (Jason Blakburn data) =================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")


# ================ ProMed data ================

# human infections
thrax1 = sf::st_read("./data/spillovers/Anthrax/PromedPoints.shp") %>%
  dplyr::filter(No_Human_I > 0)
thrax2 = sf::st_read("./data/spillovers/Anthrax/PromedPointsOther.shp") %>%
  dplyr::filter(No_Human_I > 0)
thrax = rbind(thrax1, thrax2) %>%
  dplyr::mutate(Date = as.Date(Archive_Dt, "%Y-%m-%d"),
                Year = lubridate::year(Date))

# remove americas (too data sparse)
#table(thrax$Country)
thrax = thrax %>% dplyr::filter(!Country %in% c("USA", "Uruguay", "Canada", "Peru"))

# format in fingerprint format
thrax = thrax %>%
  st_drop_geometry() %>%
  dplyr::mutate(Disease = "Anthrax", 
                Pathogen = "Bacillus anthracis",
                DataType = "point",
                Admin_name = NA,
                LocationName = City_Town,
                Source = "BlackburnProMed",
                IfPolygon_AdminLevel = NA,
                ADMcode = NA,
                ADMSource = NA,
                NumCases = No_Human_I,
                DiagnosticMetric = "Not specified",
                CasesMetric = "Confirmed cases",
                SourceType = "ProMED outbreak locations") %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Latitude, Longitude, Admin_name, LocationName, NumCases, DiagnosticMetric, CasesMetric, IfPolygon_AdminLevel,
                ADMcode, ADMSource, Source, SourceType)

write.csv(thrax, "./output/spillovers_processed/spillovers_anthrax.csv", row.names=FALSE)

