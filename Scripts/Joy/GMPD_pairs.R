#############################
## Author:  Joy Vaz       ###
## Date:    12 Feb 2019   ###
## Project: IDEAS_JPEK    ###
#############################

# set up library
library(tidyverse)
library(magrittr)
library(dplyr)
library(stringr)

# load in GMPD
gmpd <- read.csv("./Data/GMPD/GMPD_main.csv") # Rows are unique observations of parasite(ParasiteCorrectedName) occurance in a host(HostCorrectedName) for wild primates, carnivores and ungulates. Data from: Stephens et al. 2017, downloaded from https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.1799 on 9.11.2018

# extract species-level traits of 
gmpd_spp <- gmpd %>% filter(HasBinomialName == "yes") %>% 
  select(HostGroup=Group, 
         HostName=HostCorrectedName, 
         HostOrder, HostFamily, 
         HostEnvironment, 
         ParName=ParasiteCorrectedName, 
         ParType, ParPhylum, ParClass) %>% 
  mutate(HostParPair = paste(HostName, ParName, sep = ", ")) 

# create dataset of unique host-parasite pairs
gmpd_pairs <- gmpd_spp %>% 
  select(HostName, ParName, HostParPair) %>% 
  distinct()
#save as csv file
write.csv(gmpd_pairs, "./Data/JPEK/GMPD_pairs.csv")

#  calculate host range for each parasite spp
ParsHostRange <- as.data.frame(table(gmpd_pairs$ParName)) %>% 
  rename(ParName=Var1, NumHostSpp=Freq)
ParsHostRange_dist <- as.data.frame(table(ParsHostRange$NumHostSpp))
## plot
hist(ParsHostRange$NumHostSpp, breaks = "Scott")

# calculate parasite richness for each host spp
HostsParRange <- as.data.frame(table(gmpd_pairs$HostName)) %>% 
  rename(HostName=Var1, NumParSpp=Freq)
HostsParRange_dist <- as.data.frame(table(HostsParRange$NumParSpp))
## plot
hist(HostsParRange$NumParSpp, breaks = "FD")
#__________________________________________________________________________________
