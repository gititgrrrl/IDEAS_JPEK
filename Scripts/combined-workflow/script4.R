# fourth script in workflow
# get data ready for analysis - tally response variable (par species rich per host, type, etc. )
# also tally number of citations per host

library(tidyverse)
library(magrittr)

# --- load combined data from all sources (cleaned and joined in scripts 0-3) ---

rm(list=ls())
dat <- read_csv("Data/JPEK/script3.csv")

# --- adding variables for citation/sampling bias ---

dat %>%
  select(hostName, citation) %>%
  distinct() %>%                       # distinct host citation pairs
  group_by(hostName) %>%
  tally() %>% rename(numHostCitations = n)  -> hCitDf

dat <- left_join(dat, hCitDf, by="hostName")

# --- final cleaning to gather response variables and log transform variables for model --- 
# (Ellen note) Once we add in the parasites coded only to genus, then 31% of the records don't have information on parasite transmission (whether it's direct transmission only or not). 
tmp <- dat %>%
  select(hostName, parasiteName, parTransCloseOnly, parType) %>%
  distinct() # get distinct host - parasite pairs with transmission and type columns

rich_trans <- tmp %>%  # data of parasite richness in total and # parasites that have ONLY close transmission
  select(hostName, parasiteName, parTransCloseOnly) %>%
  group_by(hostName) %>%
  summarize(parRich = sum(!is.na(parasiteName)),
            parRichCloseOnly = sum(parTransCloseOnly==1, na.rm = TRUE),
            parRichNonClose = sum(parTransCloseOnly==0, na.rm = TRUE),
            parRichTransKnown = sum(!is.na(parTransCloseOnly))) # <<<<<<<<<< (Ellen) 'CloseOnly' means the only way it's transmitted is by close contact; 'NonClose' means it can be transmitted by means other than close contact

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
         absHostMeanLat = abs(hostMeanLat),
         logHostSpeciesRange = log(hostSpeciesRange),
         logHostMass = log(massKG),
         logNumHostCitations = log(numHostCitations))
# ---

write_csv(datFinal, "./Data/JPEK/script4.csv")
