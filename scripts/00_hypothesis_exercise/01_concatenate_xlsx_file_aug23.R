library(readxl)
library(tidyverse)
library(magrittr)
library(patchwork)

setwd("C:/Users/roryj/Documents/PhD/202011_fingerprint/")

# file has user names so stored externally
file <- 'SpillDriver_20230829.xlsx'

# read one by one due to an active memory issue
sheets <- lapply(c(1:25), function(j) {read_xlsx(file,
                                                 sheet = j,
                                                 range = "A1:AI19")})

lapply(c(1:length(sheets)), function(j) {
  sheets[[j]] %>%
    group_by(Driver) %>%
    pivot_longer(!Driver, names_to = "Disease", values_to = "Raw") %>%
    mutate(Author = j)
}) %>% bind_rows() -> df

df %<>% dplyr::select(Author, Disease, Driver, Raw) %>% arrange(Author, Disease)

# some name fixes
df$Disease[ grep("Monkeypox", df$Disease) ] = "Monkeypox"
df = df %>% dplyr::filter(Disease != "...35")

# blanks for any disease with at least 1 non-blank response are set to "don't know"
ddc = expand.grid(unique(df$Disease), unique(df$Author))
df2 = data.frame()
for(i in 1:nrow(ddc)){
  
  df_i = df[ df$Disease == ddc$Var1[ i ] & df$Author == ddc$Var2[i], ]
  df_i$Raw[ is.na(df_i$Raw) ] = "don't know"
  df2 = rbind(df2, df_i)
  
}
df = df2

df %<>%
  mutate(Sign = str_extract(Raw, "[A-Za-z+-]"),
         Rank = str_extract(Raw, "[0-9]")) %>%
  select(-Raw)

# summary by disease and driver
# two methods of generating hypotheses:
# 1. Majority rule (TestLib): all drivers for which more respondents stated any effect (positive or negative) than no effect (none) excluding 'don't know's
# 2. Top ranked (TestCon): all drivers that were ranked in the top 3 by at least 1 respondent
df %>%
  group_by(Disease, Driver) %>%
  summarize(Pos = sum(Sign == '+', na.rm = TRUE),
            Unk = sum(Sign == 'd', na.rm = TRUE),
            Neg = sum(Sign == '-', na.rm = TRUE),
            None = sum(Sign == 'n', na.rm = TRUE),
            Any = sum(Pos, Neg),
            TestLib = ifelse(Any >= None, "Test", "Omit"),
            Top = sum(str_detect(Rank, "[1-3]"), na.rm = TRUE),
            TestCon = ifelse(Top > 0, "Test", "Omit")) -> consensus

table(consensus$TestLib) %>% prop.table()
table(consensus$TestCon) %>% prop.table()

# combine with covariate and disease names as referenced throughout project
hyp = consensus
hyp$Driver[ grep("Precipitation", hyp$Driver)] = "Precipitation change (increase)"
hyp$Driver[ grep("Temperature", hyp$Driver)] = "Temperature change (increase)"
hyp$Driver[ grep("hunting", hyp$Driver)] = "Wildlife hunting"

# match up with covariate and disease names as referenced throughout project
hyp = hyp %>%
  dplyr::left_join(
    read.csv("./fingerprint/output/hypotheses/disease_mapping_table.csv")
  ) %>%
  dplyr::left_join(
    read.csv("./fingerprint/output/hypotheses/driver_mapping_table.csv")
  )

# save for use in modelling dataframes
write.csv(hyp, "./fingerprint/output/hypotheses/consensus_hypotheses_matched.csv", row.names=FALSE)

# summary statistics about number of users per disease
df %>%
  dplyr::filter(Sign != "d") %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(Num_Respondents = n_distinct(Author)) %>%
  dplyr::arrange(desc(Num_Respondents)) %>% 
  write.csv("./fingerprint/output/hypotheses/exercise_responsesummary.csv", row.names=FALSE)

# 
# 
# 
# 
# # ======================================================
# 
# # visualise for MS
# 
# 
# 
# # summary of responses (including disease names as stated in the exercise)
# p1 = resp_summary %>%
#   dplyr::mutate(Disease = factor(Disease, levels=rev(Disease), ordered=TRUE)) %>%
#   ggplot() + 
#   geom_bar(aes(Num_Respondents, Disease), stat="identity", fill="grey80", color="grey30") + 
#   theme_classic() + 
#   xlab("Number of respondents (out of 24)") + 
#   ylab("Disease") +
#   theme(axis.text = element_text(size=10),
#         axis.title = element_text(size=14))
# 
# # summarising disease by driver including abbreviations
# 
# abb = read.csv("./fingerprint/scripts/04_results/dz_abbrevs.csv")
# abb$Disease[12] = "Influenza A H5N1"
# 
# p2 = hyp %>%
#   dplyr::ungroup() %>%
#   dplyr::select(dz_name, Driver, TestLib, TestCon) %>%
#   dplyr::rename("Disease"=dz_name) %>%
#   reshape2::melt(id.vars = 1:2) %>%
#   dplyr::mutate(value = ifelse(value == "Test", 1, NA),
#                 variable = ifelse(variable == "TestLib", "1: Consensus (majority of respondents)", "2: Top-ranked by at least 1 respondent"),
#                 Driver = factor(Driver, levels = rev(unique(hyp$Driver)), ordered=TRUE)) %>%
#   dplyr::left_join(abb[ , c("Disease", "abbrev2")]) %>%
#   dplyr::filter(!is.na(value)) %>%
#   ggplot() + 
#   #geom_tile(aes(abbrev2, Driver, fill=value), alpha=0.8) + 
#   geom_point(aes(abbrev2, Driver, fill=variable), alpha=0.7, size=4, pch=21) +
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
# 
# 
# pc = gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(0.8, 1))
# ggsave(pc, file="./fingerprint/output/plots/SuppFigure_SpillDriver.jpg", device="jpg", units="in", dpi=300, width=14, height=10)
# 
