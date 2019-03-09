# third script
# add data for social group from sonia's collaborators for carnivores

library(tidyverse)
library(magrittr)

# --- load gmpd/pan/phyla data from script2.R and sonia collaborator data

rm(list=ls())
dat <- read_csv("Data/JPEK/script2.csv")
carnDat <- read_csv("Data/JPEK/carnivoretraitdata.csv") 

# --- select 2 columns from carnivore data --- 

carnDat %<>%
  rename(hostName=`MSW Bionomial`) %>%
  select(hostName, Groupsize) %>%
  mutate(hostName=gsub(pattern=" ", x=hostName, replacement = "_"))

# --- make new column to manipulate ---

carnDat$carnGrp <- carnDat$Groupsize
# carnDat$carnGrp[carnDat$carnGrp == "Groups"|carnDat$carnGrp == "Pairs to family clans"|carnDat$carnGrp == "Solitary or groups"] <- "group" >>>> (Ellen) see replacement line below, including change of capitalization
carnDat$carnGrp[carnDat$carnGrp %in% c("Groups", "Solitary or Groups", "Pairs to family clans", "Pairs to Groups")] <- "group" 
carnDat$carnGrp[carnDat$carnGrp %in% c("Solitary", "Pairs", "Solitary or Pairs", "Solitary?", "Solitary to Pairs", "Pairs?")] <- "non_group"
carnDat %<>% filter(!(is.na(carnGrp))) %>% select(-Groupsize)

# --- join with rest of data ---

datNew <- right_join(carnDat,dat, by="hostName") %>%
  rename(groupSizePriUng = comGrpSize,  # primates and ungulates
         groupSizePriUng_max = maxGrpSize, # <<<<<< (Ellen added)
         groupSizeCar = carnGrp)

# ---

write_csv(datNew, "./Data/JPEK/script3.csv")

