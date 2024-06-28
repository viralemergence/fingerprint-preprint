

# ================= Melioidosis =================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")



# ------------ point data from Direk Limmathurotsakul mapping paper 2016 -------------------

# occurrences from 1910-2014 in original data
mel = read.csv("./data/spillovers/Melioidosis/Data_melioidosis_1910_2014.csv") %>%
  dplyr::filter(OCCURRENCE_TYPE == "Human cases") %>%
  dplyr::select(YEAR, COUNTRY, LOCATION, LONG, LAT) %>%
  distinct() %>%
  dplyr::rename("LocationName"=LOCATION,
                "Longitude" = LONG,
                "Latitude" = LAT,
                "Country" = COUNTRY, 
                "Year" = YEAR) %>%
  dplyr::mutate(Admin_name = NA,
                Disease = "Melioidosis",
                Pathogen = "Burkholderia pseudomallei",
                DataType = "point",
                Source = "Limmathurotsakul2016",
                IfPolygon_AdminLevel = NA,
                ADMcode = NA,
                ADMSource = NA,
                NumCases = NA) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Latitude, Longitude, Admin_name, LocationName, NumCases, IfPolygon_AdminLevel,
                ADMcode, ADMSource, Source)

write.csv(mel, "./output/spillovers_processed/spillovers_melioidosis.csv", row.names=FALSE)

