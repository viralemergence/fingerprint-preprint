# Code pipeline for Gibb et al, "The anthropogenic fingerprint on emerging infectious diseases".

This repository contains all the publicly-shareable data, analysis scripts, and model outputs for the manuscript Gibb et al. "The anthropogenic fingerprint on emerging infectious diseases". This version of the code base was updated to the repository in October 2025 and covers the revised version of the manuscript (original version published on medRxiv in June 2024). The code base also contains a folder called "example_pipeline" which contains a subset of the code that can be quickly run to demonstrate functionality, for 3 diseases whose data are contained within this repo (acute Chagas, dengue, H5N1 influenza). Much of the disease data was not obtained from public-domain sources, so these datasets are not included within the repository; we have detailed a full list of these data and their sources in the manuscript.

## data

 ↳ shapefiles → Shapefiles for data processing and modelling
 
 ↳ spillovers → Disease data in their original format (13 diseases for which the data are already in the public domain)

## example_pipeline

'run_example_models.R' → Runs individual disease models for three diseases (acute Chagas, dengue, H5N1 influenza) to demonstrate the modelling pipeline functionality. Plots of fixed effects and fitted spatial fields are outputted into the "output" subfolder. 

## output 

 ↳ output → model_outputs → figures → Figures as included in the manuscript
 
 ↳ output → model_outputs → hypotheses → Outputs from the participatory hypothesis exercise

 ↳ output → model_df → Dataframes of outbreaks, background points and covariates used for model fitting

 ↳ output → model_outputs → Contains fitted fixed effect parameters ("disease_params" and "disease_params_htt"), posterior samples from individual disease models ("disease_pooled" and "disease_pooled_arboregions"), and results of the framework sensitivity check for US arboviruses ("wnv_testmodels")

 ↳ output → spillovers_processed → Processed disease datasets: each disease has a CSV of outbreaks plus an accompanying shapefile of geographic uncertainty in outbreak location

## scripts
 
 ↳ 00_hypothesis_exercise → Scripts to process coauthor-completed hypothesis exercise
 
 ↳ 01_disease_data_processing → Scripts to pre-process and harmonize all disease outbreak event data sources
 
 ↳ 02_covariate_data_processing → Scripts to pre-process certain covariate data when required
 
 ↳ 03_modelling → Scripts to run geospatial modelling pipelines (global all diseases; individual pipelines for each disease)
 
 ↳ 04_results → Scripts to summarise across models and output figures and key results


## Contact information

For enquiries contact Rory Gibb, University College London (rory.gibb.14@ucl.ac.uk).
