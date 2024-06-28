
# ================ Hantavirus cardiopulmonary syndrome ======================

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")



# ----------- 1. Surveillance reports, Brazil, 2001 to 2020 (RG digitised) --------------------

# disease reports
ha1 = read.csv("./data/spillovers/Hantavirus_SouthAmerica/brazil_moh/hantavirus_monthlyTS_2001_2020.csv") %>%
  dplyr::filter(Disease == "Hantavirus" & Cases_metric == "Confirmed") %>%
  dplyr::filter(Location_type == "Infection_municipality") %>%
  dplyr::group_by(IBGE6, Year) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE),
                   MunipName = head(MunipName, 1)) %>%
  dplyr::mutate(Disease = "Hantavirus cardiopulmonary syndrome", 
                Pathogen = "Orthohantavirus",
                DataType = "polygon",
                Country = "Brazil",
                ADMSource = "Brazil_shapefile",
                IfPolygon_AdminLevel = 2,
                Source = "Brazil DATASUS system") %>%
  dplyr::rename("Admin_name" = MunipName,
                "ADMcode" = IBGE6) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Admin_name, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

shp = sf::st_read("./data/spillovers/Hantavirus_SouthAmerica/brazil_moh/Brazil_shp_harm_2022.shp") %>%
  dplyr::select(IBGE6, name_mn) %>%
  dplyr::rename("ADMcode"=IBGE6, "Admin_name"=name_mn) %>%
  sf::st_transform(crs = 4326)
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) 

# combine with lat-lon info
ha1 = ha1 %>%
  dplyr::left_join(
    shp %>% st_drop_geometry() %>% dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% dplyr::mutate(ADMcode = as.numeric(ADMcode)) %>% dplyr::select(-Admin_name)
  ) %>%
  dplyr::arrange(Year, ADMcode)

# additional data
Encoding(ha1$Admin_name) = "latin1"
ha1$CasesMetric = "Confirmed case counts (municipality of infection)"
ha1$DiagnosticMetric = "Not specified"
ha1$SourceType = "National case surveillance system"

shp = shp %>%
  dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% dplyr::mutate(ADMcode = as.numeric(ADMcode)) %>%
  dplyr::filter(ADMcode %in% ha1$ADMcode)




# ------------ 2. Argentinian surveillance system 2015-present ----------------

# anonymised individual level data; subset to confirmed cases only
ha2 = read.csv("./data/spillovers/rodent-borne ARG/hantavirus_argentina.csv") %>%
  dplyr::mutate(Date = as.Date(Fecha_Inicio_Sintomas, "%d/%m/%Y")) %>%
  dplyr::filter(lubridate::year(Date) >= 2015) %>%
  dplyr::filter(Nombre_Etiologia == "Hantavirus" & Prueba_Resultado == "Positivo") %>%
  dplyr::group_by(Id_sivila_Estudio_Individual) %>%
  dplyr::summarise(Date = head(Date, 1),
                   Provincio = tolower(head(Provincia_Residencia, 1)),
                   Departamento = tolower(head(Departamento_Domicilio_Txt, 1)),
                   Localidad = tolower(head(Localidad_Domicilio_Txt, 1)),
                   DiagnosticMetrics = paste(unique(Pruebas), collapse=", ")) %>%
  dplyr::mutate(Departamento = replace(Departamento, Departamento == "sin dato", NA),
                Localidad = replace(Localidad, Localidad == "sin dato", NA),
                Localidad = replace(Localidad, Localidad == Departamento, NA)) # where only the 2nd admin level is known
ha2$IfPolygon_AdminLevel = 1
ha2$IfPolygon_AdminLevel[ !is.na(ha2$Departamento) ] = 2
ha2$IfPolygon_AdminLevel[ !is.na(ha2$Localidad) ] = 3

# unique province/department/locality and admin level
locs = ha2 %>% dplyr::select(Provincio, Departamento, Localidad, IfPolygon_AdminLevel) %>% distinct()

# N.B. could geolocate to a more specific location
# identify polygons
locs_1 = locs %>% dplyr::filter(IfPolygon_AdminLevel == 1)
locs_2 = locs %>% dplyr::filter(IfPolygon_AdminLevel == 2)
locs_3 = locs %>% dplyr::filter(IfPolygon_AdminLevel == 3) 
locs_2 = rbind(locs_2, locs_3) # currently only geoloc 3 to 2 becuase no shaepfiles easily avialble; could be extended to l3 divisions or geolocated

# gaul shapefile l1
gaul1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")
g1 = gaul1 %>% 
  dplyr::filter(tolower(ADM1_NAME) %in% locs_1$Provincio)
locs_1 = locs_1 %>% 
  left_join(
  g1 %>% dplyr::select(ADM1_CODE, ADM1_NAME) %>% st_drop_geometry %>% dplyr::rename("ADMcode"=ADM1_CODE, "Provincio"=ADM1_NAME) %>% dplyr::mutate(Provincio = tolower(Provincio))
)
g1 = gaul1 %>% dplyr::filter(tolower(ADM1_CODE) %in% locs_1$ADMcode)
locs_1$ADMcode = paste("GAUL1", locs_1$ADMcode, sep="_")

# gaul shapefile l2
gaul2 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")
g2 = gaul2 %>% 
  dplyr::filter(ADM0_NAME == "Argentina" & tolower(ADM2_NAME) %in% locs_2$Departamento)
testx = g2 %>% 
  dplyr::select(ADM2_CODE, ADM2_NAME, ADM1_NAME) %>%
  dplyr::rename("ADMcode" = ADM2_CODE, "Provincio" = ADM1_NAME, "Departamento" = ADM2_NAME) %>%
  dplyr::mutate(Provincio = tolower(Provincio), Departamento = tolower(Departamento)) %>% st_drop_geometry()
locs_2 = locs_2 %>% 
  left_join(testx)

# manual fix missing 
locs_2 %>% dplyr::filter(is.na(ADMcode)) %>% distinct()
locs_2 = locs_2 %>% dplyr::mutate(ADMcode = replace(ADMcode, Provincio == "buenos aires" & Departamento == "cañuelas", 82725),
                                  ADMcode = replace(ADMcode, Provincio == "buenos aires" & Departamento == "general lamadrid", 4691),
                                  ADMcode = replace(ADMcode, Provincio == "jujuy" & Departamento == "gral manuel belgrano", 4645),
                                  ADMcode = replace(ADMcode, Provincio == "salta" & Departamento == "general san martin", 4773),
                                  ADMcode = replace(ADMcode, Provincio == "santa fe" & Departamento == "general lopez", 4831),
                                  ADMcode = replace(ADMcode, Provincio == "entre rios" & Departamento == "islas de ibicuy", 4628))
g2 = gaul2 %>% dplyr::filter(tolower(ADM2_CODE) %in% locs_2$ADMcode)
locs_2$ADMcode = paste("GAUL2", locs_2$ADMcode, sep="_")

# combine together case data
ha2 = left_join(ha2, rbind(locs_1, locs_2))

# create combined shapefile
p1 = g1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = g2 %>% 
  dplyr::mutate(ADMcode = paste("GAUL2", ADM2_CODE, sep="_"),
                Admin_name = ADM2_NAME,
                Admin_level = 2) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
shp_a = rbind(p1, p2)
all(shp_a$ADMcode %in% ha2$ADMcode)

# summarise cases by location and year
ha3 = ha2 %>%
  dplyr::mutate(Year = lubridate::year(Date)) %>%
  dplyr::group_by(ADMcode, Year) %>%
  dplyr::summarise(NumCases = length(Id_sivila_Estudio_Individual), 
                   DiagnosticMetric = paste(DiagnosticMetrics, collapse=", "), 
                   DiagnosticMetric = paste(unique(unlist(strsplit(DiagnosticMetric, ", "))), collapse=", "),
                   IfPolygon_AdminLevel = head(IfPolygon_AdminLevel, 1), 
                   LocationName = head(Localidad, 1)) %>%
  dplyr::mutate(Country = "Argentina", 
                Disease = "Hantavirus cardiopulmonary syndrome", 
                Pathogen = "Orthohantavirus",
                DataType = "polygon",
                ADMSource = "GAUL",
                Source = "Argentina national surveillance system",
                CasesMetric = "Confirmed case counts",
                SourceType = "National case surveillance system")

# add shapefile lat lons and admin name
shp_a = cbind(
  shp_a, 
  st_coordinates(st_centroid(shp_a)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) 
ha3 = left_join(ha3, shp_a %>% st_drop_geometry()) %>%
  dplyr::select(-Admin_level)




# ------------ combine datasets and save -----------------

hanta = rbind(
  ha1 %>% dplyr::mutate(LocationName = Admin_name, ADMcode = as.character(ADMcode)),
  ha3 %>% dplyr::filter(ADMcode != "GAUL2_NA")
)

shp = rbind(
  shp %>% dplyr::mutate(Admin_level = 2), 
  shp_a
)
data.table::fwrite(hanta, "./output/spillovers_processed/spillovers_hantavirus.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_hantavirus_shp.R")






# # viz
# library(maptools)
# data("wrld_simpl")
# ws = sf::st_as_sf(wrld_simpl)
# ggplot() +
#   geom_sf(data=ws, fill="grey80", color=NA) +
#   maptheme +
#   geom_point(data=hanta, aes(x=Longitude, y=Latitude), col="black", size=0.25)
