

# ==================== Yellow fever =======================

# matches records to GAUL polygons

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr); library(ggplot2); library(sf)
source("./scripts/00_plot_themes.R")


# ======================= Data wrangling ====================

# Global data from Shearer 2018

# Number 1 - "abraid" - mostly extracted from admin level reports
yf1 = read.csv("./data/spillovers/YellowFever/Shearer2018/yf_abraid_20160420_clean.csv") %>%
  dplyr::select(longitude, latitude, occurrence_date, location_name, precision, admin_unit_qc_gaul_code) %>%
  dplyr::rename("Longitude"=longitude, "Latitude"=latitude,
                "Year"=occurrence_date, "LocationName"=location_name) %>%
  dplyr::mutate(Disease = "Yellow fever",
                Pathogen = "Yellow fever virus",
                DataType = ifelse(precision=="PRECISE", "point", "polygon"),
                ADM1code = replace(admin_unit_qc_gaul_code, precision != "ADMIN1", NA),
                ADM2code = replace(admin_unit_qc_gaul_code, precision != "ADMIN2", NA),
                ADM3code = NA,
                ADMSource = "GAUL") %>%
  dplyr::filter(precision != "COUNTRY") %>%
  dplyr::select(-precision, -admin_unit_qc_gaul_code) %>%
  dplyr::mutate(IfPolygon_AdminLevel = NA, 
                IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, !is.na(ADM1code), 1),
                IfPolygon_AdminLevel = replace(IfPolygon_AdminLevel, !is.na(ADM2code), 2),
                Source = "Shearer2018",
                Admin_name = "",
                NumCases = NA,
                Country="", 
                DiagnosticMetric = "Not specified",
                CasesMetric = "Outbreak location",
                SourceType = "Compiled from administrative level surveillance reports (IHME)")

# Number 2 - "wdm" - extracted from literature mainly
yf2 = read.csv("./data/spillovers/YellowFever/Shearer2018/yf_wdm_JL_FS_additions.csv") %>%
  dplyr::filter(abundance_id != "absence" & host_type_id == "human") %>%
  dplyr::select(longitude, latitude, start_year, location_type_id, full_name, GAUL_CODE, admin_level, diagnostic_method) %>%
  dplyr::rename("Longitude"=longitude, "Latitude"=latitude,
                "Year"=start_year, "LocationName"=full_name, "IfPolygon_AdminLevel"=admin_level,
                "DataType"=location_type_id, "DiagnosticMetric"=diagnostic_method) %>%
  dplyr::mutate(Disease = "Yellow fever",
                Pathogen = "Yellow fever virus",
                ADM1code = replace(GAUL_CODE, IfPolygon_AdminLevel != 1, NA),
                ADM2code = replace(GAUL_CODE, IfPolygon_AdminLevel != 2, NA),
                ADM3code = replace(GAUL_CODE, IfPolygon_AdminLevel != 3, NA),
                ADMSource = "GAUL",
                Year = replace(Year, Year==2054, 2005),  # replace erroneous start year
                Source = "Shearer2018",
                Admin_name = "",
                NumCases = NA,
                CasesMetric = "Outbreak location",
                SourceType = "Compiled from literature (IHME)",
                Country="") %>%
  dplyr::select(-GAUL_CODE) 

# combine
yf1 = rbind(yf1, yf2)


# ------ polygons --------

# GAUL polys (ADM1)
adm1_polys = yf1 %>% dplyr::filter(IfPolygon_AdminLevel == 1) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_1/g2015_2014_1/g2015_2014_1.shp")

# missing with lat-lons - extract polygon directly
missing = adm1_polys %>% dplyr::filter(!ADM1code %in% g1$ADM1_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADM1code = g1_missing$ADM1_CODE

# missing without latlons
missing2 = adm1_polys %>% dplyr::filter(!ADM1code %in% g1$ADM1_CODE) %>% dplyr::filter(is.na(Longitude))

# combine
adm1_polys = rbind(adm1_polys %>% dplyr::filter(ADM1code %in% g1$ADM1_CODE),
                   missing, 
                   missing2)

# extract (121/125)
gps = unique(adm1_polys$ADM1code)
gaul1 = g1 %>% dplyr::filter(ADM1_CODE %in% gps)
print(paste("Found", nrow(gaul1), "of", length(gps), "GAUL1 polygons", sep=" "))


# GAUL polys (ADM2)
adm2_polys = yf1 %>% dplyr::filter(IfPolygon_AdminLevel == 2) 
g1 = sf::st_read("C:/Users/roryj/Dropbox/Research/fingerprint/data/shapefiles/g2015_2014/g2015_2014_2/g2015_2014_2/g2015_2014_2.shp")

# missing with lat-lons - extract polygon directly
missing = adm2_polys %>% dplyr::filter(!ADM2code %in% g1$ADM2_CODE) %>% dplyr::filter(!is.na(Longitude))
missingx = missing; coordinates(missingx) = ~Longitude+Latitude
missingx = st_as_sf(missingx)
st_crs(missingx) = st_crs(g1)
sf::sf_use_s2(FALSE)
ii = st_intersects(missingx, g1)
g1_missing = g1[ as.numeric(ii), ]
missing$ADM2code = g1_missing$ADM2_CODE

# missing without latlons
missing2 = adm2_polys %>% dplyr::filter(!ADM2code %in% g1$ADM2_CODE) %>% dplyr::filter(is.na(Longitude))

# combine
adm2_polys = rbind(adm2_polys %>% dplyr::filter(ADM2code %in% g1$ADM2_CODE),
                   missing, 
                   missing2)

# extract (44/46)
gps = unique(adm2_polys$ADM2code)
gaul2 = g1 %>% dplyr::filter(ADM2_CODE %in% gps)
print(paste("Found", nrow(gaul2), "of", length(gps), "GAUL2 polygons", sep=" "))

# combine all data
yf1 = adm1_polys %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1code, sep="_")) %>%
  rbind(adm2_polys %>% dplyr::mutate(ADMcode = paste("GAUL2", ADM2code, sep="_"))) %>%
  rbind(yf1 %>% dplyr::filter(DataType == "point") %>% dplyr::mutate(ADMcode = "")) %>%
  rbind(yf1 %>% dplyr::filter(IfPolygon_AdminLevel == 3) %>% dplyr::mutate(DataType = "point", ADMcode = paste("GAUL3_", ADM3code, sep="_"))) %>%
  dplyr::select(-ADM1code, -ADM2code, -ADM3code)

# create composite ADMcode

# combine and harmonise shapefiles 
p1 = gaul1 %>% 
  dplyr::mutate(ADMcode = paste("GAUL1", ADM1_CODE, sep="_"),
                Admin_name = ADM1_NAME,
                Admin_level = 1) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
p2 = gaul2 %>% 
  dplyr::mutate(ADMcode = paste("GAUL2", ADM2_CODE, sep="_"),
                Admin_name = ADM2_NAME,
                Admin_level = 2) %>%
  dplyr::select(ADMcode, Admin_name, Admin_level) 
shp1 = rbind(p1, p2)
all(shp1$ADMcode %in% yf1$ADMcode)



# ---------- 2. Brazil surveillance data -----------

yf3 = read.csv("./data/spillovers/YellowFever/brazil_moh/yellowfever_monthlyTS_2001_2016.csv") %>%
  dplyr::filter(Disease == "Yellowfever" & Cases_metric == "Confirmed" & Location_type == "Infection_municipality") %>%
  dplyr::group_by(IBGE6, Year) %>%
  dplyr::summarise(NumCases = sum(NumCases, na.rm=TRUE),
                   MunipName = head(MunipName, 1)) %>%
  dplyr::mutate(Disease = "Yellow fever", 
                Pathogen = "Yellow fever virus",
                DataType = "polygon",
                Country = "Brazil",
                ADMSource = "Brazil_shapefile",
                IfPolygon_AdminLevel = 2,
                Source = "Brazil DATASUS system") %>%
  dplyr::rename("Admin_name" = MunipName,
                "ADMcode" = IBGE6) %>%
  dplyr::mutate(LocationName = Admin_name) %>%
  dplyr::select(Disease, Pathogen, DataType, Country, Year, Admin_name, LocationName, NumCases, IfPolygon_AdminLevel, ADMcode, ADMSource, Source)

shp = sf::st_read("./data/spillovers/YellowFever/brazil_moh/Brazil_shp_harm_2022.shp") %>%
  dplyr::select(IBGE6, name_mn) %>%
  dplyr::rename("ADMcode"=IBGE6, "Admin_name"=name_mn) %>%
  sf::st_transform(crs = 4326)
shp = cbind(
  shp, 
  st_coordinates(st_centroid(shp)) %>% as.data.frame %>% dplyr::rename("Longitude"=1, "Latitude"=2)
) 

# combine with lat-lon info
yf3 = yf3 %>%
  dplyr::left_join(
    shp %>% st_drop_geometry() %>% dplyr::filter(!is.na(ADMcode) & ADMcode != 0) %>% dplyr::mutate(ADMcode = as.numeric(ADMcode)) %>% dplyr::select(-Admin_name)
  ) %>%
  dplyr::arrange(Year, ADMcode)

Encoding(yf3$Admin_name) = "latin1"
yf3$CasesMetric = "Confirmed case counts (municipality of infection)"
yf3$DiagnosticMetric = "Not specified"
yf3$SourceType = "National case surveillance system"


# combine
yfv = do.call(
  rbind.data.frame, 
  list(yf1, yf3)
)

# combine shapefiles
shp = rbind(
  shp1, 
  shp %>% dplyr::select(-Longitude, -Latitude) %>% dplyr::mutate(Admin_level = 2)
)

# save data and shapefile
data.table::fwrite(yfv, "./output/spillovers_processed/spillovers_yfv.csv", row.names=FALSE)
save(shp, file="./output/spillovers_processed/spillovers_yfv_shp.R")


# # viz
# library(maptools)
# data("wrld_simpl")
# ws = st_as_sf(wrld_simpl)
# ggplot() + 
#   geom_sf(data=ws, fill="grey80", color=NA) + 
#   maptheme + 
#   geom_point(data=yfv, aes(x=Longitude, y=Latitude), col="black", size=0.25)

