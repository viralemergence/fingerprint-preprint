
# =================== visualise results of hypothesis exercise ==================

library(tidyverse)
library(magrittr)
library(patchwork)

PATH = dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd(PATH)

# user summary and consensus
users = read.csv("./output/hypotheses/exercise_responsesummary.csv")
hyp = read.csv("./output/hypotheses/consensus_hypotheses_matched.csv")
hyp$dz_name[ hyp$Disease == "Monkeypox (mpox)"] = "Monkeypox (mpox)"

# summary of responses (including disease names as stated in the exercise)
p1 = users %>%
  dplyr::mutate(Disease = replace(Disease, Disease == "Monkeypox", "Monkeypox (mpox)")) %>%
  dplyr::mutate(Disease = factor(Disease, levels=rev(Disease), ordered=TRUE)) %>%
  ggplot() + 
  #geom_bar(aes(Num_Respondents, Disease), stat="identity", fill=colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(9)[ 5 ], color="grey30") + 
  geom_segment(aes(x = 0, xend=Num_Respondents, y=Disease, yend=Disease), color=colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(9)[ 5 ], show.legend = FALSE, size=0.7) + 
  geom_point(aes(Num_Respondents, Disease), pch=21, size=3, fill=colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(9)[ 5 ], color="grey30") +
  theme_classic() + 
  xlab("Number of respondents (out of 25)") + 
  ylab("Disease") +
  ggtitle("Response rate per disease") +
  theme(axis.text = element_text(size=9),
        plot.title = element_text(size=14, hjust=0.5),
        axis.title = element_text(size=13))

# drivers that were top ranked across diseases
p2 = hyp %>%
  dplyr::group_by(Driver) %>%
  dplyr::summarise(numdiseases1 = n_distinct(Disease[ TestLib == "Test" ]),
                   numdiseases2 = n_distinct(Disease[ TestCon == "Test" ]),
                   propdis2 = numdiseases2 / n_distinct(hyp$Disease)) %>%
  dplyr::arrange(desc(numdiseases2)) %>%
  dplyr::mutate(Driver = factor(Driver, levels = rev(Driver), ordered=TRUE)) %>%
  ggplot() + 
  #geom_bar(aes(propdis2, Driver), stat="identity", fill=colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(9)[ 8 ], color="grey30") +
  geom_segment(aes(x = 0, xend=propdis2, y=Driver, yend=Driver), color=colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(9)[ 8 ], show.legend = FALSE, size=0.7) + 
  geom_point(aes(propdis2, Driver), pch=21, size=4.5, fill=colorRampPalette(MetBrewer::met.brewer("Archambault", 9))(9)[ 8 ], color="grey30") +
  theme_classic() + 
  xlab("Proportion of diseases") + 
  ylab("Driver") +
  ggtitle("Inclusion in top 3 ranked drivers") +
  theme(axis.text = element_text(size=10),
        plot.title = element_text(size=14, hjust=0.5),
        axis.title = element_text(size=13))


pc1 = gridExtra::grid.arrange(p1, p2, ncol=2)  


# matrix summarising disease by driver (including abbreviations)

abb = read.csv("./scripts/04_results/dz_abbrevs.csv")
abb$Disease[12] = "Influenza A H5N1"
abb$Disease[23] = "Monkeypox (mpox)"

# p3 = hyp %>%
#   dplyr::ungroup() %>%
#   dplyr::select(dz_name, Driver, TestLib, TestCon) %>%
#   dplyr::rename("Disease"=dz_name) %>%
#   reshape2::melt(id.vars = 1:2) %>%
#   dplyr::mutate(value = ifelse(value == "Test", 1, NA),
#                 variable = ifelse(variable == "TestLib", "1: Majority rule", "2: Top-ranked"),
#                 Driver = factor(Driver, levels = rev(unique(hyp$Driver)), ordered=TRUE)) %>%
#   dplyr::left_join(abb[ , c("Disease", "abbrev2")]) %>%
#   dplyr::filter(!is.na(value)) %>%
#   ggplot() + 
#   geom_tile(aes(abbrev2, Driver, fill=variable), alpha=0.8) + 
#   #geom_point(aes(abbrev2, Driver, fill=variable), alpha=0.7, size=4, pch=21) +
#   theme_classic() + 
#   facet_wrap(~variable, ncol=2) + 
#   scale_fill_viridis_d(na.value=NA, begin=0.35, end=0.6) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle=90, vjust=0.4), 
#         strip.text = element_text(size=15), 
#         strip.background = element_blank(),
#         axis.title = element_text(size=14),
#         panel.grid.major.x = element_line(color="grey90")) + 
#   xlab("Disease") +   
#   theme(axis.text = element_text(size=10))

p3 = hyp %>%
  dplyr::ungroup() %>%
  dplyr::select(dz_name, Driver, TestLib, TestCon, TestAny) %>%
  dplyr::rename("Disease"=dz_name) %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::mutate(value = ifelse(value == "Test", 1, NA),
                variable = as.character(variable),
                variable = replace(variable, variable == "TestAny", "1: Any author"),
                variable = replace(variable, variable == "TestLib", "2: Majority rule"),
                variable = replace(variable, variable == "TestCon", "3: Top-ranked"),
                Driver = factor(Driver, levels = rev(unique(hyp$Driver)), ordered=TRUE)) %>%
  dplyr::left_join(abb[ , c("Disease", "abbrev2")]) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot() + 
  geom_tile(aes(abbrev2, Driver, fill=variable), alpha=0.8) + 
  #geom_point(aes(abbrev2, Driver, fill=variable), alpha=0.7, size=4, pch=21) +
  theme_classic() + 
  facet_wrap(~variable, ncol=3) + 
  scale_fill_viridis_d(na.value=NA, begin=0.35, end=0.6) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, vjust=0.4, size=9), 
        strip.text = element_text(size=15), 
        strip.background = element_blank(),
        axis.title = element_text(size=14),
        panel.grid.major.x = element_line(color="grey90")) + 
  xlab("Disease") +   
  theme(axis.text = element_text(size=10))

# save combined plot
pc = gridExtra::grid.arrange(pc1, p3, nrow=2, heights=c(0.8, 1))
pc = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"), 
                           fontface = "plain", size = 28, 
                           x = c(0.001, 0.5, 0.001), y = c(0.99, 0.99, 0.55))

ggsave(pc, file="./output/plots/SuppFigure_SpillDriver.jpg", device="jpg", units="in", dpi=300, width=15, height=10, scale=0.9)

