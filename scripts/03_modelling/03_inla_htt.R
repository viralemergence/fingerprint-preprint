

# ================ Runs individual disease models using INLA geospatial regression ===================

# Wrapper script to run separate models for all individual diseases (geospatial logistic regression using R-INLA)
# Fitting with scaled linear effect of HTT; uses datasets built in 02_inla_per_disease (for Figure 5)
# Individual disease scripts stored in "./scripts/03_modelling/inla_htt/"

# setup objects
setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/fingerprint/")
library(raster); library(rgdal); library(dplyr); library(magrittr) 
library(ggplot2); library(sf); library(reticulate); library(rgee)
library(INLA); library(enmSdm)

# scripts per disease
scr = list.files("./scripts/03_modelling/inla_htt/", full.names=TRUE)

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