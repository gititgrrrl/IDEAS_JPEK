# initial GMPD clean script in workflow 
# but multiple inclusion exclusion criteria to be performed still...

# ---- 
# still to do in this script: 
# - exclude hosts/parasite not enough data (fungi prion)
# - addresses exclusion criteria from sonia
# ---- 

library(tidyverse)
library(magrittr)

# --- load super raw GMPD  ---

rm(list=ls())
GMPD_raw <- read_csv("Data/GMPD/GMPD_main.csv") # main file
GMPD_par_traits_raw <- read_csv("Data/GMPD/GMPD_parasite_traits.csv")


# --- initial GMPD Cleaning (NEEDS WORK SEE ABOVE) -----

GMPD <- GMPD_raw %>% 
  mutate(HostCorrectedName=gsub(pattern=" ", x=HostCorrectedName, replacement = "_")) %>%
  rename(hostName=HostCorrectedName, parasiteName=ParasiteCorrectedName,
         hostGroup=Group,
         parType=ParType, parPhylum=ParPhylum, prevalence=Prevalence, 
         hostsSampled=HostsSampled, samplingType=SamplingType, 
         hostEnvironment=HostEnvironment, citation=Citation) 

### --- filter data by exclusion criteria --- ###

isMissingDf <- GMPD %>% 
  filter(HasBinomialName=="no") %>%
  filter(parasiteName != "ABOLISHED") %>%
  filter(parasiteName != "not identified to genus") %>%
  filter(parasiteName != "SPLITTED in ICTV") %>%
  select(hostName, parasiteName) %>%
  distinct()
  

GMPD %<>% 
  filter(prevalence>0) %>% # parasite prevalence greater than 0
  filter(parType!="Fungus") 

### --- ADD IN INCLUSION OF PARASIETS WITH UNIQUE GENERA PER HOST --- ###

# --- select currently required columns gmpd data

GMPD %<>% select(hostGroup, hostName, parasiteName, 
                 parType, parPhylum, prevalence, 
                 hostsSampled, samplingType, 
                 hostEnvironment, citation)

parTraits <- GMPD_par_traits_raw %>% rename(parasiteName=ParasiteCorrectedName, 
                                            parTransClose = close)

GMPD <- left_join(GMPD, parTraits, by="parasiteName")

# ---- 

write_csv(GMPD, "./Data/JPEK/script0.csv")


