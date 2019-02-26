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

# load in GMPD_main.csv
gmpd <- read.csv("./Data/GMPD/GMPD_main.csv") # Rows are unique observations of parasite(ParasiteCorrectedName) occurance in a host(HostCorrectedName) for wild primates, carnivores and ungulates. Data from: Stephens et al. 2017, downloaded from https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.1799 on 9.11.2018
eid2 <- read.csv("./Data/EID2/SpeciesInteractions_EID2.csv")

# restrict to species-level traits of hosts and parasites
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
## save as csv
#write.csv(gmpd_pairs, "./Data/JPEK/GMPD_pairs.csv")

# calculate host range for each parasite spp
ParsHostRange <- as.data.frame(table(gmpd_pairs$ParName)) %>% 
  rename(ParName=Var1, NumHostSpp=Freq)
ParsHostRange_dist <- as.data.frame(table(ParsHostRange$NumHostSpp))
## plot
hist(ParsHostRange$NumHostSpp, 
     breaks = "Scott", 
     main = "Histogram: Host Species Range of GMPD Parasites", 
     xlab = "Number of Host Species")

# calculate parasite richness for each host spp
HostsParRange <- as.data.frame(table(gmpd_pairs$HostName)) %>% 
  rename(HostName=Var1, NumParSpp=Freq)
HostsParRange_dist <- as.data.frame(table(HostsParRange$NumParSpp))
## plot
hist(HostsParRange$NumParSpp, 
     breaks = "FD", 
     main = "Histogram: Parasite Species Richness of GMPD Mammals",
     xlab = "Number of Parasite Species")
#__________________________________________________________________________________

# capitalise first letter of binomial names in eid2prot to match zooscore and gmpdprot
eid2 <- read.csv("./Data/EID2/SpeciesInteractions_EID2.csv")
eid2 %<>%  filter(Publications.count >= 1) %<>% select(ParName = Cargo, HostName = Carrier, Carrier.classification, Cargo.classification, Publications.count) 
eid2$ParName <- gsub("(^[a-z])", "\\U\\1", tolower(eid2$ParName), perl = T)
eid2$HostName <- gsub("(^[a-z])", "\\U\\1", tolower(eid2$HostName), perl = T)

#__________________________________________________________________________________

# list of all gmpd hosts
gmpd_hosts <- as.tbl(as.data.frame(table(gmpd_pairs$HostName)))
length(gmpd_hosts$Var1) # 462 hosts

# list of all eid2 hosts
eid2_hosts <- as.tbl(as.data.frame(table(eid2$HostName)))
length(eid2_hosts$Var1) # 286 hosts

# compare and contrast unique host spp of the two main datasets
eid2xtrahosts <- setdiff(unique(eid2_hosts$Var1), unique(gmpd_hosts$Var1))
gmpdxtrahosts <- setdiff(unique(gmpd_hosts$Var1), unique(eid2_hosts$Var1))
sharedhosts <- intersect(unique(gmpd_hosts$Var1), unique(eid2_hosts$Var1)) 
allhosts <- unique(union(gmpd_hosts, eid2_hosts))

# check if zero
sum(length(eid2xtra), length(gmpdxtra), length(sameprot)) - length(parnames$parname)
#__________________________________________________________________________________#__________________________________________________________________________________

# compare and contrast unique parasite spp of the two main datasets
eid2xtra <- setdiff(unique(eid2prot$parname), unique(gmpd$ParName))
gmpdxtra <- setdiff(unique(gmpdprot$parname), unique(eid2prot$parname))
sameprot <- intersect(unique(gmpdprot$parname), unique(eid2prot$parname)) # 30 prots shared by both databases, 3 of them in noscore
# check if zero
sum(length(eid2xtra), length(gmpdxtra), length(sameprot)) - length(parnames$parname)
#__________________________________________________________________________________