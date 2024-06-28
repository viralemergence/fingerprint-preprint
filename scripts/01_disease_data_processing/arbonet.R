

# ================== Arbonet arboviral infections (USA) ====================

# n.b. data stored outside GitHub repo due to sharing constraints

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(readxl)
source("./fingerprint/scripts/00_plot_themes.R")



# -------------------- Arbonet surveillance data -------------------------

# US shapefile from census bureau with matching codes, only one issue
shp = sf::st_read("./arbonet/arbonet_countyarboviral/cb_2018_us_county_500k.shp") %>%
  dplyr::mutate(CTYCODE = as.integer(COUNTYFP),
                STCODE = as.integer(STATEFP)) 
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) %>%
  dplyr::mutate(ADMcode = paste(STATEFP, COUNTYFP, sep=""))
sf::st_write(shp, "./fingerprint/data/shapefiles/usa_counties.shp")

# arbonet files
ff = list.files("./arbonet/arbonet_countyarboviral/", pattern=".xls", full.names = TRUE)

# dataframe
arbo = data.frame()

# read in
for(f in ff){
  xl = readxl::read_xlsx(f)
  if(grepl("EEE", f)){ xl$Disease = "Eastern equine encephalitis"; xl$Pathogen = "Eastern equine encephalitis virus"}
  if(grepl("LAC", f)){ xl$Disease = "LaCrosse encephalitis"; xl$Pathogen = "LaCrosse virus" }
  if(grepl("POW", f)){ xl$Disease = "Powassan virus disease"; xl$Pathogen = "Powassan virus"}
  if(grepl("SLE", f)){ xl$Disease = "St. Louis encephalitis"; xl$Pathogen = "St. Louis encephalitis virus"}
  if(grepl("WNV", f)){ xl$Disease = "West Nile virus disease"; xl$Pathogen = "West Nile virus"}
  if(grepl("JC", f)){ xl$Disease = "Jamestown Canyon virus disease"; xl$Pathogen = "Jamestown Canyon virus"}
  if(grepl("CHIK", f)){ xl$Disease = "Chikungunya"; xl$Pathogen = "Chikungunya virus"}
  if(grepl("ZIK", f)){ xl$Disease = "Zika"; xl$Pathogen = "Zika virus"}
  
  if(grepl("NonNeuro", f)){ 
    xl$Presentation = "Non-neuroinvasive"
  } else{ 
    xl$Presentation = "Neuroinvasive" 
    }
  arbo = rbind(arbo, xl)
}

# summarise overall as well as neuro/nonneuro
a1 = arbo %>%
  dplyr::filter(!Disease %in% c("Chikungunya", "Zika")) %>%
  dplyr::group_by(Disease, County, Year) %>%
  dplyr::summarise(State = head(State, 1), 
                   Pathogen = head(Pathogen, 1),
                   ADMcode = head(fipscode, 1),
                   COUNT = sum(COUNT),
                   Presentation = "All cases")
a1 = a1 %>%
  rbind(
    arbo %>%
      dplyr::filter(Disease %in% c("Chikungunya", "Zika")) %>%
      dplyr::mutate(Presentation = "All cases") %>%
      dplyr::rename("ADMcode"=fipscode)
  )

# combine
arbo = arbo %>%
  dplyr::rename("ADMcode"=fipscode) %>%
  rbind(a1) %>%
  dplyr::rename("NumCases"=COUNT) %>%
  dplyr::mutate(ADMSource = "USCB_shapefile",
                DataType = "polygon",
                Country = "USA",
                IfPolygon_AdminLevel = 2,
                Source = "ArboNet") %>%
  dplyr::filter(ADMcode %in% shp$ADMcode)
  
# shapefile
shp = shp %>% 
  dplyr::select(ADMcode, NAME, Longitude, Latitude) %>%
  dplyr::rename("Admin_name"=NAME)

# combine with centroids
arbo = arbo %>%
  dplyr::left_join(
    shp %>% st_drop_geometry()
  ) %>%
  dplyr::select(Disease, Pathogen, Presentation, DataType, Country, Year, 
                Longitude, Latitude, Admin_name, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

# save
write.csv(arbo, "./arbonet/fingerprint_formatted/spillovers_arbonet.csv", row.names = FALSE)
save(shp, file="./arbonet/fingerprint_formatted/spillovers_arbonet_shp.R")




# # =========== viz ==============
# 
# arbo_all = arbo[ arbo$Presentation == "All cases", ] %>%
#   dplyr::filter(!Disease %in% c("Chikungunya", "Zika"))
# ext = extent(c(range(arbo_all$Longitude), range(arbo_all$Latitude))) + 5
# 
# library(rnaturalearth)
# ws = rnaturalearth::ne_countries(continent="North America")
# ws = crop(ws, ext)
# ws = st_as_sf(ws)
# px = ggplot() + 
#   geom_sf(data=ws, fill="grey80", color=NA) + 
#   maptheme + 
#   geom_point(data=arbo_all, aes(x=Longitude, y=Latitude, col=Disease), size=0.4, alpha=0.8) + 
#   facet_wrap(~Disease, nrow=3) +
#   theme(legend.position = "none")
# ggsave(px, file="./arbonet/arbonet_USAarbos_map.png", device="png", units="in", width=9, height=10, dpi=600)
