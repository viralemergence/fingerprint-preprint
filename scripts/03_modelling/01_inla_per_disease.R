

# ================ Runs individual disease models using INLA geospatial regression ===================

# Wrapper script to run separate models for all individual diseases (geospatial logistic regression using R-INLA)
# Individual disease scripts stored in "./scripts/03_modelling/inla_per_disease/"

# Per-disease pipeline:
# 1. Study area polygon generated (convex hull buffer around all presence points except for certain pan global diseases)
# 2. Point and polygon data standardised into shared sf format (points 5km radius buffer unless specified)
# 3. Background (pseudoabsence) points generated, spatially weighted by population
# 4. Socio-environmental and climatic covariates extracted and extraction checked
# 5. "Univariate" models fitted (each driver individually + geospatial random effect)
# 6. Hypothesis-driven models fitted for two hypothesis "types" generated using coauthor hypothesis exercise
# 6a = "Majority rule" (broader; more authors stated an effect of driver (either +/-) than no effect)
# 6b = "Top ranked" (stricter; all drivers ranked among top 3 drivers by at least 1 author)
# 7. Model output objects saved and outputted to "./output/model_results/"

# Individual disease scripts follow a standardised pipeline but a few specific changes per disease:
# 1. The number of background points was changed to anywhere betwen 2 and 8 times the number of presence points
# (varied depending on presence points number and geographical coverage), with the largest multiple for 
# diseases with a low number of presence points over a very wide geographic area  (e.g. Ebola, Marburg, Nipah),
# to ensure adequate sampling of background socio-environmental conditions
# 2. For any zoonoses with a known major livestock type affiliation, the "livestock" covariate was changed to represent
# density of that specific type (e.g. H5N1 poultry; Rift Valley fever ruminants), otherwise it represents density
# of all livestock within grid cell (from Gridded Livestock of the World)
# 3. Geospatial random effect mesh and hyperpriors tuned to ensure a good visual SPDE fit to data and
# relatively low computational cost

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr) 
library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)

# generate covariate dataset param
# TRUE = re-extract covariates for each disease
# FALSE = just re-run models (if covars already extracted?)
generate_dataset = FALSE

# scripts per disease
scr = list.files("./scripts/03_modelling/inla_per_disease/", full.names=TRUE)

# run each script in turn, print details, clear workspace after each
for(scr_i in scr){
  
  print(scr_i)
  source(scr_i)
  
  # clear workspace (keeping generate_dataset and scr)
  wspace = ls()
  wspace = wspace[ -which(wspace %in% c("scr", "generate_dataset"))]
  rm(list=wspace)
  
}

# ends