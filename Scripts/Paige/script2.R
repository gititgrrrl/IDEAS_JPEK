# second script
# join threat, calculate naive spp richness, & GMPD data with pantheria

library(tidyverse)
library(magrittr)
library(ggplot2) 

# browseURL("http://esapubs.org/archive/ecol/e090/184/metadata.htm") # link to metadata for Pantheria
# download.file(url = "http://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR93_Aug2008.txt", destfile = "Data/pantheria/pan.txt")  

pan_raw <- read_delim("Data/pantheria/pan.txt", 
                  "\t", escape_double = FALSE, na = "-999", 
                  trim_ws = TRUE)
gmpd_iucn <- read_csv("Data/JPEK/script1.csv")

par_rich <- gmpd_iucn %>%
  group_by(hostName, parasiteName) %>%
  distinct() %>% # to get distinct host parasite pairs
  group_by(hostName) %>%
  tally() %>% # to get parasite richness for each host
  rename(parRich = n)

gmpd_iucn %<>% full_join(par_rich, by = "hostName")

pan <- pan_raw %>% 
  filter(!is.na(MSW93_Binomial)) %>% 
  mutate(MSW93_Binomial = gsub(pattern = " ", x = MSW93_Binomial, replacement = "_")) %>%
  rename(hostName = MSW93_Binomial) %>%
  rename(dietBreadth = `6-1_DietBreadth`, habitatBreadth = `12-1_HabitatBreadth`,
         maxLongevity = `17-1_MaxLongevity_m`, popGrpSize = `10-1_PopulationGrpSize`,
         socGrpSize = `10-2_SocialGrpSize`, popDenChange = `27-4_HuPopDen_Change`) %>%
  select(hostName, dietBreadth, habitatBreadth,
         maxLongevity, popGrpSize,
         socGrpSize, popDenChange) 

foo <- pan %>%
  select(hostName, popGrpSize, socGrpSize) %>%
  rowwise() %>%
  mutate(comGrpSize=mean(c(socGrpSize, popGrpSize), na.rm = TRUE))

pan <- right_join(pan, foo, by="hostName")

dat_combine <- right_join(pan, gmpd_iucn, by = "hostName")

dat_combine %<>% select(-popGrpSize.x, -popGrpSize.y, -socGrpSize.x, -socGrpSize.y)

write_csv(dat_combine, "./Data/JPEK/script2.csv")

