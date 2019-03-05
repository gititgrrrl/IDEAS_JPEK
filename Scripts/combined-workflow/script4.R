# fourth script in workflow
# get data ready for analysis - tally response variable (par species rich per host, type, etc. )
# also tally number of citations per host

library(tidyverse)
library(magrittr)

# --- load combined data from all sources (cleaned and joined in scripts 0-3) ---

rm(list=ls())
dat <- read_csv("~/Desktop/IDEAS_JPEK/Data/JPEK/script3.csv")

# --- adding variables for citation/sampling bias ---

dat %>%
  select(hostName, citation) %>%
  distinct() %>%                       # distinct host citation pairs
  group_by(hostName) %>%
  tally() %>% rename(numHostCitations = n)  -> hCitDf

dat <- left_join(dat, hCitDf, by="hostName")

# --- final cleaning to gather response variables and log transform variables for model --- 

tmp <- dat %>%
  select(hostName, parasiteName, parTransClose, parType) %>%
  distinct() # get distinct host - parasite pairs with transmission and type columns

rich_trans <- tmp %>%  # data of parasite richness in total and # parasites that have close transmission
  select(hostName, parasiteName, parTransClose) %>%
  group_by(hostName) %>%
  summarize(parRich = sum(!is.na(parasiteName)),
            parRich_close = sum(parTransClose, na.rm = TRUE)) %>%
  mutate(propClose=parRich_close/parRich)

rich_type <- tmp %>%  # data of parasite richness by type
  select(hostName, parType) %>%
  group_by(hostName, parType) %>%
  summarize(parRich = n()) %>%
  distinct() %>%
  group_by(hostName) %>%
  spread(parType, parRich) 

colnames(rich_type)[-1] <- paste0(tolower(colnames(rich_type)[-1]), "Rich" )
rich_type[is.na(rich_type)] <- 0

datFinal <- dat %>%
  left_join(rich_trans, by = "hostName") %>%
  left_join(rich_type, by = "hostName") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(logParRich = log(parRich),
         logHostSpeciesRange = log(hostSpeciesRange),
         logHostMass = log(massG),
         logNumHostCitations = log(numHostCitations))

datFinal$groupSizePriUng[is.nan(datFinal$groupSizePriUng)] <- NA

# ---

write_csv(datFinal, "./Data/JPEK/script4.csv")
