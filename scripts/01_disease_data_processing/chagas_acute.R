

# ================ Acute Chagas disease ======================

# use Google Earth Engine to access GAUL polygons to create associated polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf); library(reticulate); library(rgee)
source("./scripts/00_plot_themes.R")

# # install miniconda / setup python dependencies
# # numpy and ee (earth engine api for python)
#rgee::ee_install()
#rgee::ee_install_upgrade()
#ee_check()

# # 1. Initialize the Python Environment  
ee_Initialize()



# ------------ 1. Browne et al 2017 (IHME) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5387921/ ---------------

# chagas acute cases (Browne 2017)
cha1 = read.csv("./data/spillovers/Chagas/browne2017/2_acute_humans.csv") %>%
  dplyr::select(Start_year, Country, Geometry_type, Site_name_1, Longitude_1, Latitude_1, Coordinate_source_1, Admin2_1, Admin1_1, Diagnosis)

# chagas occurrence cases
# cha2 = read.csv("./data/spillovers/Chagas/browne2017/3_occurrence_humans.csv") %>%
#   dplyr::select(Start_year, Country, Geometry_type, Site_name_1, Longitude_1, Latitude_1, Coordinate_source_1, Admin2_1, Admin1_1, Diagnosis)

# admin1 and admin2 polygons
adm1_polys = unique(cha1$Admin1_1); adm1_polys = adm1_polys[ !is.na(adm1_polys) ]
adm2_polys = unique(cha1$Admin2_1); adm2_polys = adm2_polys[ !is.na(adm2_polys) ]

# adm1 polys from GEE
gee1 = ee$FeatureCollection("FAO/GAUL/2015/level1")$
  filter( 
    ee$Filter$Or(
      ee$Filter$inList('ADM1_CODE', adm1_polys)
    )
  )
gee1 = ee_as_sf(gee1)

# adm2 polys from GEE
gee2 = ee$FeatureCollection("FAO/GAUL/2015/level2")$
  filter( 
    ee$Filter$Or(
      ee$Filter$inList('ADM2_CODE', adm2_polys)
    )
  )
gee2 = ee_as_sf(gee2)

# combine into chagas polygons
p1 = gee1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = gee2 %>% 
  dplyr::mutate(ADMcode = paste("GAUL2", ADM2_CODE, sep="_"),
                Admin_name = ADM2_NAME,
                Admin_level = 2) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
shp = rbind(p1, p2)

# shp %>%
#   ggplot() +
#   geom_sf(aes(fill=Admin_level))

# standardise polygon reference in chagas data
cha1 = do.call(
  rbind.data.frame,
  list(
    cha1 %>% dplyr::filter(Geometry_type == "polygon") %>% dplyr::filter(!is.na(Admin2_1)) %>% dplyr::mutate(ADMcode = paste("GAUL2", Admin2_1, sep="_"), IfPolygon_AdminLevel=2),
    cha1 %>% dplyr::filter(Geometry_type == "polygon") %>% dplyr::filter(!is.na(Admin1_1)) %>% dplyr::mutate(ADMcode = paste("GAUL1", Admin1_1, sep="_"), IfPolygon_AdminLevel=1),
    cha1 %>% dplyr::filter(Geometry_type == "point") %>% dplyr::mutate(ADMcode = NA, IfPolygon_AdminLevel=NA)
  )
) %>%
  dplyr::select(-Admin2_1, -Admin1_1, -Coordinate_source_1) %>%
  dplyr::rename("Year"=Start_year,
                "DataType"=Geometry_type,
                "Admin_name"=Site_name_1,
                "Longitude" = Longitude_1,
                "Latitude" = Latitude_1,
                "DiagnosticMetric"=Diagnosis) %>%
  dplyr::mutate(Disease = "Chagas disease (acute)", 
                Pathogen = "Trypanosoma cruzi",
                Source = "Browne2017", 
                ADMSource = "GAUL", 
                NumCases = NA)

# access lat-lon centroids from polygons
shp = shp %>%
  cbind(
    st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
  )
cha1 = rbind(
  cha1 %>% 
    dplyr::filter(!is.na(ADMcode)) %>% 
    dplyr::select(-Longitude, -Latitude) %>% 
    dplyr::left_join(shp %>% st_drop_geometry() %>% dplyr::select(ADMcode, Longitude, Latitude)),
  cha1 %>% 
    dplyr::filter(is.na(ADMcode))
)

# arrange columns
cha1 = cha1 %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Longitude, Latitude, Admin_name, DiagnosticMetric, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)
shp1 = shp %>% sf::st_cast()

# add column on data format
cha1$CasesMetric = "Outbreak location"
cha1$SourceType = "Compiled from literature (IHME)"

# encoding
Encoding(cha1$Admin_name) = "latin1"



# ----------- 2. Acute Chagas surveillance reports, Brazil, 2006 to 2020 --------------------

# disease reports from municipality of infection (not notification)
ch2 = read.csv("./data/spillovers/Chagas/brazil_moh/chagas_acute_monthlyTS_2001_2020.csv") %>%
  dplyr::filter(Disease == "Chagas_acute" & Cases_metric == "Confirmed") %>%
  dplyr::filter(Location_type == "Infection_municipality") %>%
  dplyr::group_by(IBGE6, Year) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE),
                   MunipName = head(MunipName, 1)) %>%
  dplyr::mutate(Disease = "Chagas disease (acute)", 
                Pathogen = "Trypanosoma cruzi",
                DataType = "polygon",
                Country = "Brazil",
                ADMSource = "Brazil_shapefile",
                IfPolygon_AdminLevel = 2,
                Source = "Brazil DATASUS system") %>%
  dplyr::rename("Admin_name" = MunipName,
                "ADMcode" = IBGE6) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Admin_name, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)
Encoding(ch2$Admin_name) = "latin1"

shp = sf::st_read("./data/spillovers/Chagas/brazil_moh/Brazil_shp_harm_2022.shp") %>%
  dplyr::select(IBGE6, name_mn) %>%
  dplyr::rename("ADMcode"=IBGE6, "Admin_name"=name_mn)
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
)

# combine with lat-lon info
ch2 = ch2 %>%
  dplyr::left_join(
    shp %>% st_drop_geometry() %>% dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% dplyr::mutate(ADMcode = as.numeric(ADMcode)) %>% dplyr::select(ADMcode, Longitude, Latitude)
  ) %>%
  dplyr::arrange(Year, ADMcode)

# add data format
ch2$CasesMetric = "Confirmed case counts (municipality of infection)"
ch2$DiagnosticMetric = "Not specified"
ch2$SourceType = "National case surveillance system"



# ------------ combine datasets and save -----------------

# combine chagas data
cha = rbind(cha1, ch2)
write.csv(cha, "./output/spillovers_processed/spillovers_chagas.csv", row.names=FALSE)

# combine shapefiles
shp2 = shp %>%
  dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% 
  dplyr::mutate(ADMcode = as.numeric(ADMcode), Admin_level=2) %>%
  dplyr::filter(ADMcode %in% ch2$ADMcode) %>%
  st_transform(crs = crs(shp1))
  
shp = rbind(shp1, shp2)
#sf::st_write(shp, "./output/spillovers_processed/spillovers_chagas.shp", append=FALSE) # issue with st_write
save(shp, file="./output/spillovers_processed/spillovers_chagas_shp.R")



# viz
# library(maptools)
# data("wrld_simpl")
# ws = sf::st_as_sf(wrld_simpl)
# ggplot() +
#   geom_sf(data=ws, fill="grey80", color=NA) +
#   maptheme +
#   geom_point(data=ch2, aes(x=Longitude, y=Latitude), col="black", size=0.25)
