

# ============ Summary figure of effect sizes and directions across all individual disease models ===============

# Generates Figure 3 from manuscript

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

library(raster); library(dplyr); library(magrittr); library(ggplot2); library(sf)
library(INLA)
source("./scripts/00_plot_themes.R")

params = list.files("./output/model_outputs/disease_pooled/", pattern=".csv", full.names=TRUE)
params = params[ -grep("infocriteria", params)]
params = params[ grep("summary", params)]
pp = do.call(rbind.data.frame, lapply(params, read.csv)) 
pp$mean = round(pp$mean, 3)
pp$lower = round(pp$lower, 3)
pp$upper = round(pp$upper, 3)
pp$sig = pp$lower < 0 & pp$upper < 0 | pp$lower > 0 & pp$upper > 0

# add abbrevs
pp = pp %>% left_join( read.csv("scripts/04_results/dz_abbrevs.csv") ) %>%
  dplyr::mutate(Disease = abbrev2)

pp = pp %>%
  dplyr::filter(param != "Intercept") %>%
  dplyr::filter(!param %in% c("evi_coefvar", "precip_anomaly", "tmean_anomaly")) %>%
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

# # add NAs for non-modelled vars
# pp = pp %>% dplyr::select(param, Disease, type, mean, sig, Num_observations)
# add_nas = expand.grid(unique(pp$param), unique(pp$Disease), unique(pp$type)) %>%
#   as.data.frame() %>%
#   dplyr::rename("param" = 1, "Disease"=2, "type"=3) %>%
#   dplyr::mutate(sig = NA, mean = NA) %>% 
#   dplyr::left_join(pp %>% dplyr::select(Disease, Num_observations) %>% distinct())
# add_nas = add_nas[ ! paste(add_nas$param, add_nas$Disease, add_nas$type, sep="_") %in% paste(pp$param, pp$Disease, pp$type, sep="_"), ]
# pp = rbind(pp, add_nas)

# order by n points
np = pp %>% dplyr::select(Disease, Num_observations) %>% dplyr::arrange(Num_observations) %>% distinct()
pp$Disease = factor(pp$Disease, levels = np$Disease, ordered=TRUE)

# comparison of # of disease-driver combinations actually tested between diseases
# after removing collinear variables, including a priori confounders etc
# 78% with any; 70% with majority rule; 50% with top-ranked
combs = pp %>% 
  dplyr::filter(param != "Intercept") %>%
  dplyr::select(type) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::arrange(desc(Freq))
potential_n = n_distinct(pp$param[ pp$param != "Intercept" ]) * n_distinct(pp$Disease)
combs$prop_total = combs$Freq / potential_n

# variable type
pp$vartype = "Land use intensity"
pp$vartype[ pp$param %in% c("Urban cover", "Healthcare travel time (log)") ] = "Detection"
pp$vartype[ pp$param %in% c("Temperature change", "Precipitation change") ] = "Climate change"
pp$vartype[ pp$param %in% c("Social vulnerability", "Livestock density (log)") ] = "Socioeconomic"
pp$vartype[ pp$param %in% c("Forest cover", "Cropland cover", "Landscape heterogeneity", "Biodiversity Intactness Index") ] = "Ecosystem structure"



# ============== specify criterion =============

library(patchwork)



# ========== create combined figure across all diseases =============

# run for each hypothesis set
hyps = unique(pp$type)

for(i in 1:length(hyps)){
  
  hyp = hyps[i]
  
  rr = pp %>%
    dplyr::filter(type == hyp) %>%
    dplyr::group_by(param) %>%
    dplyr::rename("signif"=sig) %>%
    dplyr::summarise(
      sig = n_distinct(Disease[ signif == TRUE]),
      pos = n_distinct(Disease[ signif == TRUE & mean > 0 ]),
      neg = n_distinct(Disease[ signif == TRUE & mean < 0 ]),
      tested = n_distinct(Disease),
      prop_tested = tested / n_distinct(pp$Disease),
      null = tested - sig,
      vartype = head(vartype, 1)
    ) %>%
    dplyr::arrange(desc(sig)) %>%
    dplyr::mutate(
      rank = rank(-sig, ties.method = "min")
    )
  
  rank_drivers =  rr$param
  
  # plot 0 = rank
  p0 = rr %>%
    dplyr::mutate(grp = "a",
                  param = factor(param, levels=rev(rank_drivers), ordered=TRUE),
                  param_cd = as.numeric(param)) %>%
    ggplot() + 
    geom_hline(yintercept=seq(0.5, length(rank_drivers)+0.5), linewidth=0.3, color="grey85") +
    geom_text(aes(grp, param_cd, label=rank), fontface="bold") + 
    theme_classic() + 
    xlab("Rank") + 
    ylab("Driver") +
    theme(axis.title.y = element_text(size=14),
          panel.border=element_blank(),
          axis.text.y = element_text(size=12),
          plot.title = element_text(size=12.25, hjust=0.5, face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(color="grey40"),
          axis.title.x = element_text(size=12, color="white"),
          axis.text.x = element_text(color="white", size=12),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color="grey40")) + 
    theme(panel.background = element_rect(fill="white")) + 
    scale_y_continuous(breaks=1:length(rank_drivers), labels=rev(rank_drivers), limits = c(0.5, length(rank_drivers)+0.5)) +
    ggtitle("Rank")
    
  # plot 1 = directionality - including "not tested"
  p1 = rr %>%
    dplyr::mutate(not_tested = n_distinct(pp$Disease) - tested) %>%
    tidyr::pivot_longer(cols=c("pos", "neg", "null", "not_tested")) %>%
    dplyr::mutate(param = factor(param, levels=rev(rank_drivers), ordered=TRUE),
                  param_cd = as.numeric(param),
                  name = replace(name, name == "pos", "Pos"),
                  name = replace(name, name == "neg", "Neg"),
                  name = replace(name, name == "null", "None"),
                  name = replace(name, name == "not_tested", "NT"),
                  value = replace(value, value==0, NA),
                  name = factor(name, levels=c("NT", "Neg", "None", "Pos"), ordered=TRUE)) %>%
    ggplot() + 
    geom_hline(yintercept=seq(0.5, length(rank_drivers)+0.5), linewidth=0.3, color="grey85") +
    geom_point(aes(name, param_cd, fill=name, size=value, group=name), pch=21, alpha=1) + 
    theme_classic() + 
    xlab("Directionality") + 
    theme(axis.title.y = element_blank(),
          plot.title = element_text(size=12.25, hjust=0.5, face = "bold"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(color="grey40"),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11),
          #panel.border = element_rect(color="grey75", fill=NA),
          panel.border=element_blank(),
          axis.ticks.x = element_line(color="grey40"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12, color="white"),
          axis.text.x = element_text(color="black", size=11)) +
    scale_y_continuous(breaks=1:length(rank_drivers), labels=rank_drivers, limits = c(0.5, length(rank_drivers)+0.5)) +
    scale_fill_manual(values=c("Pos"="violetred3", "Neg"="steelblue2", "None"="white", "NT"="grey85"), name="Direction", guide="none") + 
    scale_size(name = "Num. diseases", breaks=c(1, 5, 10, 15, 20, 25, 30), range=c(1, 6)) +
    theme(panel.background = element_rect(fill="white")) +
    ggtitle("Effect direction")
  
  
  # plot 2 = effect sizes
  cc = pp %>%
    dplyr::filter(type == hyp)
  cc$param = factor(cc$param, levels=rev(rank_drivers), ordered=TRUE)
  cc$vartype = factor(cc$vartype, levels=c("Detection", "Ecosystem structure", "Socioeconomic", "Climate change", "Land use intensity"))
  
  # param numeric
  cc$param_cd = as.numeric(cc$param)
  
  # set seed for jitter
  set.seed(1200)
  p2 = ggplot() + 
    geom_rect(data=cc[ !duplicated(cc$param), ], aes(ymin=param_cd - 0.5, ymax=param_cd + 0.5, xmin=-Inf, xmax=Inf, fill=vartype), alpha=0.4) +
    geom_hline(yintercept=seq(0.5, length(rank_drivers)+0.5), linewidth=0.3, color="white") +
    geom_vline(xintercept=0, lty=2, color="grey50", linewidth=0.5) + 
    #geom_linerange(data=cc, aes(y=param, xmin=lower, xmax=upper, group=Disease, color=vartype), alpha=0.25, size=0.8) +
    #geom_point(data=cc[ cc$sig == FALSE, ], aes(y=param_cd, x=mean, group=Disease), fill="grey20", pch=21, size=1.5, position=position_dodge(width=0.2), alpha=0.7) +
    #geom_point(data=cc[ cc$sig == TRUE, ], aes(y=param_cd, x=mean, group=Disease), color="grey20", pch=21, size=1.5, position=position_dodge(width=0.2), alpha=0.7) + 
    #ggforce::geom_sina(data=cc, aes(y=param_cd, x=mean, group=paste(Disease, param_cd), pch=sig), fill="grey97", color="grey10", size=3, alpha=0.75) + 
    geom_point(data=cc[ cc$sig == FALSE, ], aes(y=param_cd, x=mean, group=Disease), pch=21, fill="grey97", color="grey10", position=position_jitter(height=0.16), size=1.7, alpha=0.8) + 
    geom_point(data=cc[ cc$sig == TRUE, ], aes(y=param_cd, x=mean, group=Disease), pch=21, fill="grey20", color="grey10", position=position_jitter(height=0.16), size=2.2, alpha=0.7) + 
    theme_classic() + 
    scale_x_continuous(breaks=seq(-3, 3, by=1), labels=seq(-3, 3, by=1), limits=c(-max(abs(range(cc$mean)))-0.15, max(abs(range(cc$mean)))+0.15)) +
    scale_y_continuous(breaks=1:length(rank_drivers), labels=rank_drivers, limits = c(0.5, length(rank_drivers)+0.5)) +
    #xlim(-max(abs(range(c(cc$lower, cc$upper)))), max(abs(range(c(cc$lower, cc$upper))))) +
    #xlim(-max(abs(range(cc$mean)))-0.15, max(abs(range(cc$mean)))+0.15) +
    xlab("Z-score scaled effect on outbreak event risk (log odds)") + 
    scale_fill_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") +
    scale_color_manual(values = colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(5)[c(1, 2, 5, 4, 3)], name="Driver type") + 
    theme(axis.title.y = element_blank(),
          plot.title = element_text(size=12.25, hjust=0.5, face = "bold"),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(color="grey40"),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11),
          #panel.border = element_rect(color="grey75", fill=NA),
          panel.border = element_blank(),
          axis.ticks.x = element_line(color="grey40"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12),
          axis.text.x = element_text(color="black", size=11)) +
    scale_shape_manual(values=c(21, 19)) +
    theme(panel.background = element_rect(fill="white")) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.6) ) ) +
    ggtitle("Posterior mean effect size per disease")
  
  # create combined plot with shared legend
  combined = p0 + p2 + p1 & theme(legend.position = "right")
  
  # including/excluding "not tested"
  combined = combined + plot_layout(guides = "collect", widths=c(0.2, 1.8, 0.65)) # incl
  #combined = combined + plot_layout(guides = "collect", widths=c(0.2, 0.55, 2.5)) # excl
  combined
  
  # save
  if(hyp=="Majority rule"){ hn = "majorityrule" }
  if(hyp=="Top-ranked"){ hn = "topranked" }
  if(hyp=="Any"){ hn = "anyauthor" }

  filename = paste("./output/plots/Figure3_", hn, "_CS.jpg", sep="")
  ggsave(combined, file=filename, device="jpg", units="in", width=12, height=5.5, dpi=600)
  
}




