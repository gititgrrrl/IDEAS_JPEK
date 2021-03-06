# first script 
# join GMPD with phylacine by host name to get host threat status of hosts in GMPD

library(tidyverse)
library(magrittr)

# --- load GMPD (cleaned by script 0) ---

rm(list=ls())
GMPD_raw <- read_csv("Data/JPEK/script0.csv") 
phyla_raw <- read_csv("Data/phylacine/Trait_data.csv")

# --- select required columns phylacine data ---

phyla <- select(phyla_raw, "Binomial.1.2", "Order.1.2", "Family.1.2", 
                "Genus.1.2", "IUCN.Status.1.2", "Mass.g") %>%
  rename(hostName=Binomial.1.2, 
         hostOrder=Order.1.2,
         hostFamily=Family.1.2,
         hostGenus=Genus.1.2,
         IUCN=IUCN.Status.1.2,
         massG=Mass.g) %>%
  mutate(massKG=massG/1000) %>% # (Ellen) changed colname to KG b/c dividing by 1000
  mutate(combIUCN = ifelse(IUCN %in% c("CR", "EN", "VU"), "threatened", ifelse(IUCN %in% c("NT", "LC"), "not_threatened", NA))) # 03/28/19, EC changed to conform with more typical division of threat vs. non-threat


# --- Then we right join matching rows from phyla to GMPD ---

GMPD_threat <- right_join(phyla, GMPD_raw, by="hostName")

# --- 
write_csv(GMPD_threat, "./Data/JPEK/script1.csv")
