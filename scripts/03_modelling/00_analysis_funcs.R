

# ===================== FUNCTIONS FOR RUNNING ANALYSIS PIPELINE ===========================


# =================== functions for building dataset ==========================


# ------------ creates buffered convex polygon around occurrence points accounting for coastline -----------

#' @param lon vector of point longitudes
#' @param lat vector of point latitudes
#' @param custom_boundary sf polygon to use as custom boundary (for country-level surveillance data)
createStudyAreaPolygon = function(lat, lon, custom_boundary=NA, buffer=TRUE, 
                                  extent_border = 5, buf_size = 180000){

  if(class(custom_boundary)[1]=="sf"){
    
    study_area = custom_boundary
    
  } else{
    
    # study_area = sf::st_read("data/shapefiles/ne_10m_land.shp") %>%
    #   sf::as_Spatial() %>%
    #   crop(extent(c(xmin=min(lon), xmax=max(lon), ymin=min(lat), ymax=max(lat)))+5) %>%
    #   rgeos::gUnaryUnion() %>%
    #   sf::st_as_sf()
    
    sf::sf_use_s2(FALSE)
    study_area = sf::st_read("data/shapefiles/ne_10m_land.shp") %>%
      sf::st_crop(extent(c(xmin=min(lon), xmax=max(lon), ymin=min(lat), ymax=max(lat)))+extent_border) %>%
      sf::st_combine()
    sf::sf_use_s2(TRUE)
    
  }
  
  if(buffer == TRUE){
    
    # create buffered convex hull around points to clip study area
    buf_sa = st_as_sf(
      data.frame(Longitude=lon, Latitude=lat), 
      coords = c("Longitude", "Latitude"), 
      crs = 4326
    ) %>%
      st_union() %>%
      st_convex_hull() %>%
      st_buffer(buf_size) %>% # buffer around the outside
      rmapshaper::ms_simplify(keep = 0.03) %>%
      smoothr::smooth(method = "chaikin")
    
    # clip coastlines 
    sf::sf_use_s2(FALSE)
    study_area = sf::st_intersection(study_area, buf_sa) 
    if(!"data.frame" %in% class(study_area)){
      study_area = st_as_sf(data.frame(id = 1), geometry=study_area)
    }
    sf::sf_use_s2(TRUE)
    
    
  } else{
    
    if(!"data.frame" %in% class(study_area)){
      study_area = st_as_sf(data.frame(id = 1), geometry=study_area)
    }
    
  }

  # print and return    
  sf::sf_use_s2(FALSE)
  print(
    ggplot() + 
      geom_sf(data=study_area, fill="grey90", col="grey80") + 
      geom_point(data=data.frame(Longitude=lon, Latitude=lat), aes(Longitude, Latitude), col="red", size=1) +
      maptheme
  )
  sf::sf_use_s2(TRUE)
  
  return(study_area)
}


# ----------- creates buffers of radius buf_size_km around each lat-lon point (for points data) ------------

createPointBuffers = function(points_data, buf_size_km){
  
  # point data and give ADMcode
  pts = points_data
  pts$ADMcode = paste("Pts", 1:nrow(pts), sep="_")
  
  # use specified buffer size if provided, otherwise buf_size in metres
  if("BufferRadius" %in% names(pts)){
    pts$BufferRadius = replace(pts$BufferRadius, is.na(pts$BufferRadius), buf_size_km) # any not specified to standard size
  } else{
    pts$BufferRadius = buf_size_km
  }
  buf_radius = pts$BufferRadius * 1000 # convert to metres
  
  # create circular radius buffers around each point
  pts_sf = sf::st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326)
  buf = sf::st_buffer(pts_sf, buf_radius) %>%
    dplyr::left_join(pts %>% dplyr::select(record_id, Longitude, Latitude)) %>%
    dplyr::select(-BufferRadius)
  
  return(buf)
}



# ------- generates harmonised sf dataset of polygons (either original polygons, or buffers placed around points) ----------

#' @param dz disease dataframe containing point and polygon information
#' @param dz_shp shapefile associated with disease dataframe referenced via ADMcode column
#' @param point_buffer_km specify buffer radius to place around each point location
#' @param filter_area_thresh if want to remove datapoints associated with too large a geographical area, remove any with area > this km

harmoniseDataTypes = function(dz, dz_shp, point_buffer_km=5, filter_area_thresh=NA){
  
  # if dataset contains points data, convert to circular buffers of specified radius
  print("Harmonising datasets")
  if(any(dz$DataType == "point")){
    
    # create point buffers using function "createPointBuffers"
    buf_rad = point_buffer_km
    pts = createPointBuffers(dz %>% dplyr::filter(DataType == "point"), buf_rad)
    
    # combine with polygons if any contained within dataset
    if(any(dz$DataType == "polygon")){
      
      if(dz$Disease[1] %in% c("Rift Valley Fever")){ sf::sf_use_s2(FALSE) } # fix for spherical geom issues for some pathogens
      polys = dz_shp %>%
        dplyr::filter(ADMcode %in% dz$ADMcode) %>%
        dplyr::select(ADMcode) %>%
        distinct() %>%
        dplyr::full_join(dz[ dz$DataType == "polygon", ])
      if("BufferRadius" %in% names(polys)){
        polys = polys %>% dplyr::select(-BufferRadius)
      }
      dzx = rbind(pts, polys)
      sf::sf_use_s2(TRUE)
      
    } else{
      dzx = pts
    }
    
    # reporting
    print(
      paste(sum(dzx$DataType=="point"), "points,", sum(dzx$DataType == "polygon"), "polygons", sep=" ")
    )
    print(paste("Processed rows is equal to original rows:", nrow(dz) == nrow(dzx), sep=" "))
    dz = dzx
    
    # if no points, only process polygons
  } else{
    
    if(dz$Disease[1] %in% c("Rift Valley Fever")){ sf::sf_use_s2(FALSE) } # fix for spherical geom issues for some pathogens
    dzx = dz_shp %>%
      dplyr::filter(ADMcode %in% dz$ADMcode) %>%
      dplyr::select(ADMcode) %>%
      distinct() %>%
      dplyr::full_join(dz[ dz$DataType == "polygon", ])
    print(
      paste(sum(dzx$DataType=="point"), "points,", sum(dzx$DataType == "polygon"), "polygons", sep=" ")
    )
    print(paste("Processed rows is equal to original rows:", nrow(dz) == nrow(dzx), sep=" "))
    dz = dzx
    sf::sf_use_s2(TRUE)
    
  }
  
  # filter out polygons above a threshold geographical size (bc too imprecise)
  # currently 5000km2 (~70km square)
  if(!is.na(filter_area_thresh)){
    print("Excluding areas above threshold size")
    if(dz$Disease[1] %in% c("Rift Valley Fever")){ sf::sf_use_s2(FALSE) } # fix for spherical geom issues for some pathogens
    dz$poly_area = as.vector(sf::st_area(dz)/10^6)
    sf::sf_use_s2(TRUE)
    dz = dz %>% dplyr::filter(poly_area <= filter_area_thresh)
  }
  
  return(dz)
}



# --------- generate background points across study area either randomly or population weighted -----------

#' @param num_points number of points to generate
#' @param study_area polygon covering extent of model area
#' @param buffer_radius radius of polygon buffers to generate
#' @param method method for spatially generating background points, either popweight (weighted by log population) or random (random across study area)
generateBackgroundPolygons = function(num_points, study_area, buffer_radius, method="popweight", seed=NA, pop_res="1km"){
  
  print(paste("Generating", num_points, "background points using method:", method, sep=" "))
  
  # set random seed for reproducibility
  if(!is.na(seed)){ set.seed(seed) }
  if(!method %in% c("popweight", "random")) return("Error: method must be 'popweight' or 'random")
  if(!pop_res %in% c("1km", "10km")) return("Error: pop_res must be 1km or 10km")
  
  if(method == "popweight"){
    
    # create log population raster
    if(pop_res == "1km"){
      pd = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_1km_Aggregated.tif") %>% crop(study_area)
    }
    if(pop_res == "10km"){
      pd = raster("C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/population_wp/ppp_2010_10km_agg.tif") %>% crop(study_area)
    }
    pd = raster::mask(pd, fasterize::fasterize(study_area, raster(pd), field = "id"))
    pd = log(pd + 1)
    vv = values(pd)
    vv[ vv == 0 ] = median(vv[ vv != 0 & !is.na(vv)])
    values(pd) = vv
    
    # sample random points (num_points) across log pop raster (could have enmSdm dependency removed)
    pas = enmSdm::sampleRast(pd, n=num_points, replace=TRUE, prob=TRUE)
    if(is.null(pas)){ pas = enmSdm::sampleRast(pd, n=npts, replace=TRUE, prob=TRUE) } # replicate to deal with buggy
    pas = as.data.frame(pas) %>%
      dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
      dplyr::mutate(record_id = paste("bg", 1:num_points, sep="_"))
  }
  
  if(method == "random"){
    
    pas = sp::spsample(as_Spatial(study_area), n=num_points, type="random") %>%
      as.data.frame() %>%
      dplyr::rename("Longitude"=1, "Latitude"=2)  %>%
      dplyr::mutate(record_id = paste("bg", 1: num_points, sep="_"))
  }
  
  # generate buffers
  print("Generating background point buffers")
  bg = st_as_sf(pas, coords = c("Longitude", "Latitude"), crs = 4326)
  bg = st_buffer(bg, buffer_radius) %>%
    dplyr::left_join(pas %>% dplyr::select(record_id, Longitude, Latitude)) %>%
    dplyr::mutate(ADMcode = paste("bg", 1:num_points, sep="_"))
  return(bg)
  
  # remove seed 
  rm(.Random.seed, envir=.GlobalEnv)
}



# ========================= functions to extract covariates ============================

# ------------ extract Hansen data from Google Earth Engine ------------

#' @param shp sf shapefile object
extractHansenGEE = function(shp, type='loss'){
  
  # initialise ee python env
  library(reticulate)
  library(rgee)
  ee_Initialize()
  
  #https://developers.google.com/earth-engine/tutorials/tutorial_forest_03
  # hansen loss: use "mean" function to extract proportion of forest loss pixel within polygon/buffer region around point
  #hansen <- ee$Image('UMD/hansen/global_forest_change_2021_v1_9')$select(type)
  
  # forest base year = % of cell covered by forest in 2000; # loss year = did a loss event occur? (complete stand removal)
  # multiply together to get amount of cell that was forest that was lost (e.g. if 25% was lost)
  # then extract with mean == gets % of total area that is lost forest area
  baseyear <- ee$Image('UMD/hansen/global_forest_change_2021_v1_9')$select('treecover2000')
  loss <- ee$Image('UMD/hansen/global_forest_change_2021_v1_9')$select('loss')
  loss_pc <- baseyear$multiply(loss)
  
  # iterate and extract in batches of ~100 (aovids over-requesting for such a dense raster)
  shp_unique = shp %>% dplyr::select(ADMcode) %>% distinct()
  n_grps = round(nrow(shp_unique)/100)
  shp_unique$batch = sample(1:n_grps, size=nrow(shp_unique), replace=TRUE)
  
  # run iterative extraction
  shp_hans = data.frame()
  for(i in 1:n_grps){
    
    cat(paste(i, "...", sep=""))
    shp_i = shp_unique %>% dplyr::filter(batch == i)
    extr_i <- ee_extract(x = loss_pc, y = shp_i["ADMcode"], sf = FALSE, fun=ee$Reducer$mean())
    extr_i <- extr_i %>%
      dplyr::rename("forest_loss"=2) %>%
      dplyr::mutate(forest_loss = forest_loss/100)
    shp_hans <- rbind(shp_hans, extr_i)
  }
  
  shp = left_join(shp, shp_hans)
  ext = shp %>% dplyr::select(forest_loss) %>% st_drop_geometry()
    
  return(ext)
}




# ------------ extract specified covariates ------------

# march 2023 - RG to double check sf projection correctly harmonising to raster projecton

extractCovariateVals = function(sf_data, covariate_name){
  
  # location where all covariates are stored
  predictors_location = "C:/Users/roryj/Dropbox/Research/fingerprint/data/predictors/"
  
  # era5 land temperature anomaly 1991-2020 (cf baseline 1950-1990)
  if(covariate_name == "tmean_anomaly"){
    ras = raster::raster(paste(predictors_location, "era5_land/temperature_metrics/temperature_meananomaly_19912020.tif", sep=""))
    names(ras) = "tmean_anomaly"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # era5 land precip anomaly 1991-2020 (cf baseline 1950-1990)
  if(covariate_name == "precip_anomaly"){
    ras = raster::raster(paste(predictors_location, "era5_land/precip_metrics/precip_meananomaly_19912020.tif", sep=""))
    names(ras) = "precip_anomaly"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # mean annual temp
  if(covariate_name == "tmean"){
    ras = raster::raster(paste(predictors_location, "era5_land/temperature_metrics/temperature_annualmean_20002020.tif", sep=""))
    names(ras) = "tmean"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # mean annual precip
  if(covariate_name == "precip"){
    ras = raster::raster(paste(predictors_location, "era5_land/precip_metrics/precip_meanannual_20002020.tif", sep=""))
    names(ras) = "precip"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  
  # era5 land mean annual temperature difference 2001-2020 compared to 1950-1970 baseline
  if(covariate_name == "tmean_change"){
    ras = raster::raster(paste(predictors_location, "era5_land/temperature_metrics/temperature_meanchange_baselinetopresent.tif", sep=""))
    names(ras) = "tmean_change"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # era5 land mean annual precip difference 2001-2020 compared to 1950-1970 baseline
  if(covariate_name == "precip_change"){
    ras = raster::raster(paste(predictors_location, "era5_land/precip_metrics/precip_meanchange_baselinetopresent.tif", sep=""))
    names(ras) = "precip_change"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # copernicus 100m tree cover
  if(covariate_name == "forest_cover"){
    ras = raster::raster(paste(predictors_location, "landcover_copernicus/PROBAV_LC100_global_v3.0.1_2015-base_Tree-CoverFraction-layer_EPSG-4326.tif", sep=""))
    names(ras) = "forest_cover"
    ext = exactextractr::exact_extract(ras, sf_data)
    cfunc = function(x, ...){
      x = x %>% dplyr::filter(value != 255)
      return( sum((x$value/100) * x$coverage_fraction, na.rm=TRUE) / sum(x$coverage_fraction, na.rm=TRUE))
    }
    ext = sapply(ext, cfunc)
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # hansen "global forest change" dataset 2000-2019
  if(covariate_name == "forest_loss"){
    if(disease_name %in% c("rvf")){ sf::sf_use_s2(FALSE) } # fix for spherical geom issues for some pathogens
    ext = extractHansenGEE(sf_data)
    sf::sf_use_s2(TRUE)
  }
  
  # copernicus 100m crops cover
  if(covariate_name == "crop_cover"){
    ras = raster::raster(paste(predictors_location, "landcover_copernicus/PROBAV_LC100_global_v3.0.1_2015-base_Crops-CoverFraction-layer_EPSG-4326.tif", sep=""))
    names(ras) = "crop_cover"
    ext = exactextractr::exact_extract(ras, sf_data)
    cfunc = function(x, ...){
      x = x %>% dplyr::filter(value != 255)
      return( sum((x$value/100) * x$coverage_fraction, na.rm=TRUE) / sum(x$coverage_fraction, na.rm=TRUE))
    }
    ext = sapply(ext, cfunc)
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # crops expansion (% area with net cropland gain 2000-2019)
  if(covariate_name == "crop_expansion"){
    ras = raster::raster(paste(predictors_location, "cropland_expansion/Global_cropland_3km_netgain.tif", sep=""))
    names(ras) = "crop_expansion"
    ext = exactextractr::exact_extract(ras, sf_data)
    cfunc = function(x, ...){
      x = x %>% dplyr::filter(value != 255)
      return( sum((x$value/100) * x$coverage_fraction, na.rm=TRUE) / sum(x$coverage_fraction, na.rm=TRUE))
    }
    ext = sapply(ext, cfunc)
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # copernicus 100m urban cover
  if(covariate_name == "urban_cover"){
    ras = raster::raster(paste(predictors_location, "landcover_copernicus/PROBAV_LC100_global_v3.0.1_2015-base_BuiltUp-CoverFraction-layer_EPSG-4326.tif", sep=""))
    names(ras) = "urban_cover"
    ext = exactextractr::exact_extract(ras, sf_data)
    cfunc = function(x, ...){
      x = x %>% dplyr::filter(value != 255)
      return( sum((x$value/100) * x$coverage_fraction, na.rm=TRUE) / sum(x$coverage_fraction, na.rm=TRUE))
    }
    ext = sapply(ext, cfunc)
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # urban expansion from ESA-CCIs
  if(covariate_name == "urban_expansion"){
    ras = raster::raster(paste(predictors_location, "esa_cci/esacci_urbanexpansion_masked_20002018.tif", sep=""))
    names(ras) = "urban_expansion"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # mining (% cover per grid cell) https://www.nature.com/articles/s41597-022-01547-4
  if(covariate_name == "mining"){
    ras = raster::raster(paste(predictors_location, "mining/global_miningarea_v2_30arcsecond.tif", sep=""))
    names(ras) = "mining"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # local Biodiversity Intactness Index
  if(covariate_name == "biodiv_intact"){
    ras = raster::raster(paste(predictors_location, "biodiversity_intactness/local_bii_index.tif", sep=""))
    names(ras) = "biodiv_intact"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # protected areas
  if(covariate_name == "protected_areas"){
    ras = raster::raster(paste(predictors_location, "protected_areas/WDPA_Sep2022_1km_ras_final.tif", sep=""))
    names(ras) = "protected_areas"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # hunting defaunation index
  if(covariate_name == "hunting"){
    ras = raster::raster(paste(predictors_location, "defaunation (benitez lopez)/DefInd.tif", sep=""))
    #ref = raster::raster(paste(predictors_location, "protected_areas/WDPA_Sep2022_1km_ras_final.tif", sep=""))
    proj4string(ras) = CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
    #ras = raster::projectRaster(from=ras, to=CRS(ref))
    names(ras) = "hunting"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # EVI first-order coefficient of variation (variance) at 1km2 res (derived from 250m base layers) http://www.earthenv.org/texture
  if(covariate_name == "evi_coefvar"){
    ras = raster::raster(paste(predictors_location, "tuanmu_heterogeneity/cv_01_05_1km_uint16.tif", sep=""))
    names(ras) = "evi_coefvar"
    ras = crop(ras, extent(sf_data)+0.25) * 0.0001 # scale to real values
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # EVI second-order dissimilarity http://www.earthenv.org/texture; sensitive to heterogeneity in human-dominated land classes
  if(covariate_name == "evi_dissimilarity"){
    ras = raster::raster(paste(predictors_location, "tuanmu_heterogeneity/Dissimilarity_01_05_1km_uint32.tif", sep=""))
    names(ras) = "evi_dissimilarity"
    ras = crop(ras, extent(sf_data)+0.25) * 0.0001 # scale to real values
    values(ras)[ values(ras) > 200000 ] = NA 
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # if(covariate_name == "livestock_keepers"){
  #   ras = raster::raster(paste(predictors_location, "fao_livestockkeepers/rplk.grd", sep=""))
  #   names(ras) = "livestock_keepers"
  #   ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
  # }
  
  # combined density of livestock (cattle, sheep, goats, buffalo, pigs)
  if(covariate_name == "livestock_all"){
    ras = raster::stack(list.files(paste(predictors_location, "livestock (glw3)/unmodelled_aw/", sep=""), pattern=".tif", full.names=TRUE))
    ras = ras[[ grep("Ct|Gt|Sh|Bf|Pg", names(ras)) ]]
    ras = raster::crop(sum(ras), extent(sf_data)+0.25)
    ras = ras / raster::area(ras)
    names(ras) = "livestock_all"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # livestock ruminants (cattle, goats, sheep, buffalo)
  if(covariate_name == "livestock_ruminants"){
    ras = raster::stack(list.files(paste(predictors_location, "livestock (glw3)/unmodelled_aw/", sep=""), pattern=".tif", full.names=TRUE))
    ras = ras[[ grep("Ct|Gt|Sh|Bf", names(ras)) ]]
    ras = crop(sum(ras), extent(sf_data)+0.25)
    ras = ras / raster::area(ras)
    names(ras) = "livestock_ruminants"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # livestock poultry (chicken/duck)
  if(covariate_name == "livestock_poultry"){
    ras = raster::stack(list.files(paste(predictors_location, "livestock (glw3)/unmodelled_aw/", sep=""), pattern=".tif", full.names=TRUE))
    ras = ras[[ grep("Ch|Dk", names(ras)) ]]
    ras = crop(sum(ras), extent(sf_data)+0.25)
    ras = ras / raster::area(ras)
    names(ras) = "livestock_poultry"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  if(covariate_name == "livestock_cattle"){
    ras = raster::raster(paste(predictors_location, "livestock (glw3)/unmodelled_aw/6_Ct_2010_Aw.tif", sep=""))
    names(ras) = "livestock_cattle"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  if(covariate_name == "livestock_horses"){
    ras = raster::raster(paste(predictors_location, "livestock (glw3)/unmodelled_aw/6_Ho_2010_Aw.tif", sep=""))
    names(ras) = "livestock_horses"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  if(covariate_name == "livestock_pigs"){
    ras = raster::raster(paste(predictors_location, "livestock (glw3)/unmodelled_aw/6_Pg_2010_Aw.tif", sep=""))
    names(ras) = "livestock_pigs"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # if(covariate_name == "livestock_camels"){
  #   # ras = raster::raster(paste(predictors_location, "", sep=""))
  #   # names(ras) = "livestock_camels"
  #   # ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
  # }
  
  if(covariate_name == "health_travel"){
    ras = raster::raster(paste(predictors_location, "healthcare_travel/2020_motorized_travel_time_to_healthcare/2020_motorized_travel_time_to_healthcare.geotiff", sep=""))
    names(ras) = "health_travel"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  if(covariate_name == "social_vulnerability"){
    ras = raster::raster(paste(predictors_location, "socialvulnerability_global/socialvul_mean_fact20.tif", sep=""))
    names(ras) = "social_vulnerability"
    ext = exactextractr::exact_extract(ras, sf_data, fun='mean')
    ext = data.frame(x = ext); names(ext)[1] = names(ras)
  }
  
  # return
  return(ext)
}



