library(readxl)
library(tidyverse)
library(magrittr)
library(patchwork)

setwd("C:/Users/roryj/Documents/Research/202011_fingerprint/")

# file has user names so stored externally
file <- 'SpillDriver.xlsx'

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
df$Disease[ grep("Monkeypox", df$Disease) ] = "Monkeypox (mpox)"
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
# two methods of generating hypotheses (initial MS)
# 1. Majority rule (TestLib): all drivers for which more authors stated any effect (positive or negative) than no effect (none) excluding 'don't know's
# 2. Top ranked (TestCon): all drivers that were ranked in the top 3 by at least 1 author
# 3. Any hypothesised (TestAny): any drivers that were hypothesised to have a non-zero effect by at least 1 author
df %>%
  group_by(Disease, Driver) %>%
  summarize(Pos = sum(Sign == '+', na.rm = TRUE),
            Unk = sum(Sign == 'd', na.rm = TRUE),
            Neg = sum(Sign == '-', na.rm = TRUE),
            None = sum(Sign == 'n', na.rm = TRUE),
            Any = sum(Pos, Neg),
            TestLib = ifelse(Any >= None, "Test", "Omit"),
            Top = sum(str_detect(Rank, "[1-3]"), na.rm = TRUE),
            TestCon = ifelse(Top > 0, "Test", "Omit"),
            TestAny = ifelse(Any > 0, "Test", "Omit"),
            TestAll = ifelse(Any > 0 & None == 0, "Test", "Omit")) -> consensus

table(consensus$TestAny) %>% prop.table()
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

