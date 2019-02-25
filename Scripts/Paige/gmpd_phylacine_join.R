# script to join GMPD with phylacine by host name to get host threat status of hosts in GMPD
library(tidyverse)
library(magrittr)
library(ggplot2) 

GMPD_raw <- read_csv("Data/GMPD/GMPD_main.csv")
phyla_raw <- read_csv("Data/phylacine/Trait_data.csv")

# The host matching variable in GMPD is HostCorrected (format = Genus species)
# in phyla it is Binomial.1.2 (formate = Genus_species)
# so we need to change those formats to match

GMPD <- GMPD_raw %>% 
  filter(HasBinomialName == "yes") %>% 
  mutate(HostCorrectedName=gsub(pattern=" ", x=HostCorrectedName, replacement = "_")) %>%
  rename(hostName=HostCorrectedName, parasiteName=ParasiteCorrectedName) 

# select currently required columns gmpd data
GMPD %<>% select(Group, hostName, parasiteName, 
                 ParType, ParPhylum, Prevalence, 
                 HostsSampled, SamplingType, 
                 HostEnvironment, Citation)

# select curretnly required columns phylacine data
phyla <- select(phyla_raw, "Binomial.1.2", "Order.1.2", "Family.1.2", "Genus.1.2", "IUCN.Status.1.2", "Mass.g") %>%
  rename(hostName=Binomial.1.2)

# Then we *should* be able to right join matching rows from phyla to GMPD
GMPD_threat <- right_join(phyla, GMPD, by="hostName")

write_csv(GMPD_threat, "./Data/JPEK/GMPD_threat.csv")
