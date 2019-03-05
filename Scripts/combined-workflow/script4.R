# fourth script in workflow
# get data ready for analysis
# add bias column

library(tidyverse)
library(magrittr)

# --- load combined data from all sources (cleaned and joined in scripts 0-3)

rm(list=ls())
dat <- read_csv("~/Desktop/IDEAS_JPEK/Data/JPEK/script3.csv")

# adding variables for citation/sampling bias
dat %>%
  group_by(hostName) %>%
  distinct(Citation) %>%  
  tally() %>% rename(numHostCitations=n) -> hCitDf

dat %>%
  group_by(parasiteName) %>%
  distinct(Citation) %>%  
  tally() %>% rename(numParCitations=n) -> pCitDf  

dat <- left_join(dat, hCitDf, by="hostName")
dat <- left_join(dat, pCitDf, by="parasiteName")

# final cleaning to unify variable names etc
LH <- read_delim("~/Google Drive/GMPD/GlobalParasites/Data/LH/LH.txt", "\t", escape_double = FALSE, na = "-999", trim_ws = TRUE) # Set NA's to -999

LH_sub <- LH[, c(5, 6, 18, 23, 29, 32, 41, 44, 47)] # create a new data frame with a subset of the columns
names(LH_sub) <- c("HostName", "HostActivCycle", "HostHomeRange", "HostMaxLifespan", "HostGroupSize", "HostTrophic", "HostSpeciesRange", "HostMeanLat", "HostMeanLong")
LH_sub$HostName <- gsub(" ", "_", LH_sub$HostName) # replace space with underscore
LH_sub %<>% 
  mutate(HostActivCycle = factor(HostActivCycle, labels = c("nocturnal", "crepuscular", "diurnal")),
         HostTrophic = factor(HostTrophic, labels=c("herbivore", "omnivore", "carnivore")), 
         HostSpeciesRange = HostSpeciesRange/1000,
         HostMaxLifespan = HostMaxLifespan/12)

Dat <- read_csv("Data/JPEK/dat-residual-corrected.csv") # errors on import, but unrelated to the columns we are interested in
TempDat <- Dat %>%
  rename(IUCN = IUCN.Status.1.2,
         HostName = hostName,
         HostGroup = Group,
         ParName = parasiteName,
         ParTrans_Close = close,
         NumHostCitations = numHostCitations,
         UngPrimMeanGroupSize = comGrpSize, 
         CarnGroupSize = carnGrp) %>%
  mutate(HostMass = Mass.g/1000,
         HostIUCNcomb = ifelse(IUCN %in% c("CR", "EN", "VU", "NT"), "threatened", ifelse(IUCN=="LC", "not_threatened", NA))) %>% # had 1 data deficient & 1 EP--convert those to NA
  select(HostName, HostMass, HostIUCNcomb, HostGroup, ParName, ParType, ParTrans_Close, NumHostCitations, UngPrimMeanGroupSize, CarnGroupSize) %>%
  distinct()
TempDat$CarnGroupSize <- tolower(TempDat$CarnGroupSize)

# NOTE: ~10% of host-parasite combinations don't have parasite transmission information
RichnessByTrans <- TempDat %>%
  select(HostName, ParName, ParTrans_Close) %>%
  group_by(HostName) %>%
  summarize(ParRich = sum(!is.na(ParName)),
            ParRich_Close = sum(ParTrans_Close, na.rm = TRUE))
RichnessByParType <- TempDat %>%
  select(HostName, ParType) %>%
  group_by(HostName, ParType) %>%
  summarize(ParRich = n()) %>%
  distinct() %>%
  group_by(HostName) %>%
  spread(ParType, ParRich) %>%
  select(-Prion, -Fungus) # these categories don't have enough data
colnames(RichnessByParType)[-1] <- paste0("ParRich_", colnames(RichnessByParType)[-1])
RichnessByParType[is.na(RichnessByParType)] <- 0

FinalDat <- TempDat %>%
  select(-ParName, -ParType, -ParTrans_Close) %>%
  distinct() %>%
  left_join(RichnessByTrans, by = "HostName") %>%
  left_join(RichnessByParType, by = "HostName") %>%
  left_join(LH_sub, by = "HostName") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(logParRich = log(ParRich),
         logHostSpeciesRange = log(HostSpeciesRange),
         logHostMass = log(HostMass),
         logNumHostCitations = log(NumHostCitations),
         propParCloseTrans = ParRich_Close/ParRich)
FinalDat$UngPrimMeanGroupSize[is.nan(FinalDat$UngPrimMeanGroupSize)] <- NA
saveRDS(FinalDat, "~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen/FinalDat.RDS")


write_csv(dat, "./Data/JPEK/script4.csv")
