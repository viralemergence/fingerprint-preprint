
# --------------- Estimate pooled mean effect size across individual disease models (Figure 4) -------------

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(INLA)
source("./scripts/00_plot_themes.R")
source("./scripts/03_modelling/00_inla_funcs.R")
source("./scripts/03_modelling/00_analysis_funcs.R")

# disease types
dz = read.csv("scripts/04_results/dz_abbrevs.csv")

# recover parameter samples across all disease models 
params = list.files("./output/model_outputs/disease_pooled/", pattern=".csv", full.names=TRUE)
params = params[ grep("samp", params) ]

# function to get pooled estimates for subset of diseases
getPooledEstimates = function(hyp_type, disease_subset="All"){
  
  # subset to disease group
  if(disease_subset == "ZVB"){
    params = params[ grep("cchf|chagas|eee|jac|jev|lacrosse|lyme|mayaro|oropouche|pknowlesi|plague|powassan|rickettsia|rvf|sle|wnv|yfv", params) ]
  } else if(disease_subset == "ZD"){
    params = params[ grep("anthrax|ebola|h5n1|hantavirus|junin|lassa|marburg|mers|mpx|nipah", params) ]
  } 
  
  # subset to hypothesis type
  params = params[ grep(hyp_type, params)]
  
  # read in posterior samples from all diseases
  pp = lapply(params, read.csv)
  
  # get parameter names
  pnames = unique(unlist(lapply(pp, names)))
  pnames = pnames[ -grep("Disease|type|id", pnames)]
  
  # for storing results
  result = data.frame()
  
  # for each parameter, recover overall estimate
  # take N sets of 1 sample per disease, calculate mean each time
  for(pname in pnames){
    
    # reporting
    print(pname)
    
    # get just that parameter from each
    getParam = function(x){
      x = x[ , which(names(x) == pname), drop=FALSE ]
      if(nrow(x)>0){ return(x) }
    }
    samples = lapply(pp, getParam)
    samples = samples[ which( sapply(samples, ncol) > 0 ) ]
    
    # cbind into a data frame
    samples_c = do.call(cbind.data.frame, samples)
    
    # calculate the mean across N unique combinations of 1 sample per disease
    mean_ests = apply(
      #apply(samples_c, 2, function(x) sample(x, 5000, replace=TRUE)),
      samples_c,
      1,
      mean
    )
    
    # plot histogram
    hist(mean_ests, 100)
    
    if(hyp_type=="majorityrule"){ hn = "Majority rule" }
    if(hyp_type=="topranked"){ hn = "Top-ranked" }
    if(hyp_type=="_any_"){ hn = "Any author" }

    # summarise
    res = data.frame(
      hyp = hn,
      disease_subset = disease_subset,
      param = substr(pname, 1, nchar(pname)-2),
      num_diseases = length(samples), 
      median = median(mean_ests),
      lower_95 = quantile(mean_ests, 0.025),
      lower_67 = quantile(mean_ests, 0.1666),
      upper_67 = quantile(mean_ests, 0.8333),
      upper_95 = quantile(mean_ests, 0.975)
    )
    result = rbind(result, res)
    
  }
  
  return(result)
}

# get pooled estimates for majority rule
res1 = do.call(
  rbind.data.frame,
  list(getPooledEstimates(hyp_type = "majorityrule", disease_subset = "All"),
       getPooledEstimates(hyp_type = "majorityrule", disease_subset = "ZVB"),
       getPooledEstimates(hyp_type = "majorityrule", disease_subset = "ZD"))
)

# get pooled estimates for top-ranked
res2 = do.call(
  rbind.data.frame,
  list(getPooledEstimates(hyp_type = "topranked", disease_subset = "All"),
       getPooledEstimates(hyp_type = "topranked", disease_subset = "ZVB"),
       getPooledEstimates(hyp_type = "topranked", disease_subset = "ZD"))
)

# get pooled estimates for any author
res3 = do.call(
  rbind.data.frame,
  list(getPooledEstimates(hyp_type = "_any_", disease_subset = "All"),
       getPooledEstimates(hyp_type = "_any_", disease_subset = "ZVB"),
       getPooledEstimates(hyp_type = "_any_", disease_subset = "ZD"))
)

# plot figures for each hypothesis type
for(i in 1:3){
  
  if(i == 1){
    result = res1
  } 
  if(i == 2){
    result = res2
  }
  if(i == 3){
    result = res3
  }
  
  result = result %>%
    dplyr::mutate(param = replace(param, param == "health_travel_log", "Healthcare travel time (log)"),
                  param = replace(param, param == "livestock_log", "Livestock density (log)"),
                  param = replace(param, param == "livestock_ruminants_log", "Livestock density (log)"),
                  param = replace(param, param == "hunting", "Hunting Pressure Index"),
                  param = replace(param, param == "forest_cover", "Forest cover"),
                  param = replace(param, param == "forest_loss", "Forest loss"),
                  param = replace(param, param == "crop_expansion", "Cropland expansion"),
                  param = replace(param, param == "evi_dissimilarity", "Landscape heterogeneity"),
                  param = replace(param, param == "urban_expansion", "Urban expansion"),
                  param = replace(param, param == "urban_cover", "Urban cover"),
                  param = replace(param, param == "crop_cover", "Cropland cover"),
                  param = replace(param, param == "mining", "Mining cover"),
                  param = replace(param, param == "tmean_change", "Temperature change"),
                  param = replace(param, param == "precip_change", "Precipitation change"),
                  param = replace(param, param == "social_vulnerability", "Social vulnerability"),
                  param = replace(param, param == "biodiv_intact", "Biodiversity Intactness Index"),
                  param = replace(param, param == "protected_areas", "Protected area cover"),
                  param = replace(param, param == "popdens_log", "Population density (log)"))
  
  result$vartype = "Land use intensity"
  result$vartype[ result$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
  result$vartype[ result$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
  result$vartype[ result$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
  result$vartype[ result$param %in% c("Forest cover", "Cropland cover", "Landscape heterogeneity", "Biodiversity Intactness Index") ] = "Ecosystem structure"
  result$vartype = factor(result$vartype, levels=c("Detection", "Ecosystem structure",  "Socioeconomic",  "Climate change", "Land use intensity"))
  
  # exclude any with <5 estimates to synthesis
  x_est = 5
  result_underx = result %>% dplyr::filter(num_diseases < x_est)
  result_underx$median = NA
  result_underx$lower_95 = NA
  result_underx$lower_67 = NA
  result_underx$upper_67 = NA
  result_underx$upper_95 = NA
  result = rbind(
    result %>% dplyr::filter(num_diseases >= x_est),
    result_underx
  )
  
  # missing variables
  msg = as.data.frame(table(result$disease_subset, result$param)) %>%
    dplyr::filter(Freq == 0)
  for(i in 1:nrow(msg)){
    ms_i = result %>% dplyr::filter(disease_subset == msg$Var1[i]) %>% dplyr::slice_head(n=1)
    ms_i$param = msg$Var2[i]
    ms_i$num_diseases = 0
    ms_i$median = NA
    ms_i$lower_95 = NA
    ms_i$lower_67 = NA
    ms_i$upper_67 = NA
    ms_i$upper_95 = NA
    ms_i$vartype = as.character(result$vartype[ result$param == ms_i$param][1])
    result = rbind(result, ms_i)
  }
  
  # remove variables that are NA for all subsets for viz purposes
  all_nas = result %>% 
    dplyr::group_by(param) %>%
    dplyr::summarise(all_na = all(is.na(median))) %>%
    dplyr::filter(all_na == TRUE)
  result = result %>%
    dplyr::filter(!param %in% all_nas$param)
  
  fac_order = result %>% 
    dplyr::filter(disease_subset == "All") %>%
    dplyr::arrange(vartype, desc(median))
  
  result$param = factor(result$param, levels=rev(fac_order$param), ordered=TRUE)
  result$param_cd = as.numeric(result$param)
  
  result = result %>%
    dplyr::mutate(
      disease_subset2 = replace(disease_subset, disease_subset=="All", "All diseases [n=31]"),
      disease_subset2 = replace(disease_subset2, disease_subset2=="ZD", "Zoonotic (direct) [n=10]"),
      disease_subset2 = replace(disease_subset2, disease_subset2=="ZVB", "Zoonotic (vector-borne) [n=17]")
    )
  result = result %>%
    dplyr::mutate(
      disease_subset2 = factor(disease_subset2, 
                               levels=c("All diseases [n=31]",
                                        "Zoonotic (vector-borne) [n=17]",
                                        "Zoonotic (direct) [n=10]"),
                               ordered=TRUE)
    )
  
  # facetted
  
  # hide 95% CI on high uncertainty variables
  result = result %>%
    dplyr::mutate(
      lower_95 = replace(lower_95, lower_95 < -2 | upper_95 > 2, NA),
      upper_95 = replace(upper_95, lower_95 < -2 | upper_95 > 2, NA)
    )
  

  px = result %>%
    ggplot() + 
    geom_rect(aes(ymin=param_cd - 0.5, ymax=param_cd + 0.5, xmin=-Inf, xmax=Inf, fill=vartype), alpha=0.4) +
    geom_hline(yintercept=seq(0.5, nrow(fac_order)+0.5), linewidth=0.2, color="white") +
    geom_vline(xintercept = 0, lty=2, color="grey30") +
    geom_linerange(aes(y=param_cd, xmin=lower_95, xmax=upper_95), color="grey25", linewidth=0.5) + 
    geom_linerange(aes(y=param_cd, xmin=lower_67, xmax=upper_67), color="grey25", linewidth=1.2) +
    geom_point(aes(x=median, y=param_cd), color="grey20", size=2.2) + 
    geom_text(aes(x=max(upper_95, na.rm=TRUE)+0.07, y=param_cd, label=paste("(", num_diseases, ")", sep="")), size=2, color="grey30", fontface="italic") + 
    theme_classic() +
    scale_y_continuous(breaks=1:length(fac_order$param), labels=rev(fac_order$param), limits = c(0.5, length(fac_order$param)+0.5)) +
    scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
    xlab("Mean marginal effect size across all tested diseases (log odds)") +
    facet_wrap(~disease_subset2, nrow=1) +
    theme(strip.background = element_blank(), strip.text = element_text(size=12)) +
    scale_color_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
    theme(legend.position="bottom") +
    theme(axis.title.y = element_blank(),
          axis.line.y = element_blank()) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.6) ) ) 
  
  if(result$hyp[1]=="Majority rule"){ hn = "majorityrule" }
  if(result$hyp[1]=="Top-ranked"){ hn = "topranked" }
  if(result$hyp[1]=="Any author"){ hn = "anyauthor" }

  ggsave(px, file=paste("./output/plots/Figure4_PooledIndivEffects_", hn, "_CS.jpg", sep=""), device="jpg", units="in", width=10, height=4.5, dpi=600)

}





