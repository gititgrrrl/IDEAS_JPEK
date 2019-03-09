# second script
# join threat, calculate naive spp richness, & GMPD data with pantheria for HOST TRAITS

# ---- 
# still to do in this script: 
# - consistency check with Ellen code
# - determine which host traits actually needed at this point for project
# - get rid of not required variables
# ---- 

library(tidyverse)
library(magrittr)

# --- Browse metadata about pantheria ---
# browseURL("http://esapubs.org/archive/ecol/e090/184/metadata.htm") 
# download.file(url = "http://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR93_Aug2008.txt", destfile = "Data/pantheria/pan.txt")  

# --- load pantheria and script 1 data 

rm(list=ls())
pan_raw <- read_delim("Data/pantheria/pan.txt", 
                  "\t", escape_double = FALSE, na = "-999", 
                  trim_ws = TRUE)
gmpd_iucn <- read_csv("Data/JPEK/script1.csv")

# --- Clean host trait data --- (NEEDS CHECKS SEE ABOVE)

pan <- pan_raw %>% 
  rename(hostName = MSW93_Binomial) %>%
  mutate(hostName = gsub(pattern = " ", x = hostName, replacement = "_")) %>%
  rename(hostMaxLifespan = `17-1_MaxLongevity_m`, 
         hostHomeRange=`22-1_HomeRange_km2`, 
         popGrpSize = `10-1_PopulationGrpSize`,
         socGrpSize = `10-2_SocialGrpSize`, 
         hostSpeciesRange = `26-1_GR_Area_km2`,
         popDenChange = `27-4_HuPopDen_Change`,
         hostMeanLat=`26-4_GR_MRLat_dd`,
         hostMeanLong=`26-7_GR_MRLong_dd`) %>%
  select(hostName, hostHomeRange, hostMaxLifespan, 
         socGrpSize, popGrpSize, hostSpeciesRange, 
         hostHomeRange, hostMeanLat, hostMeanLong) %>%
  mutate(hostSpeciesRange = hostSpeciesRange/1000,
         hostMaxLifespan = hostMaxLifespan/12) 

# --- get combined population/social group size data as recommended by Sonia --- 

tmp <- pan %>%
  select(hostName, popGrpSize, socGrpSize) %>%
  rowwise() %>%
  mutate(comGrpSize=mean(c(socGrpSize, popGrpSize), na.rm = TRUE), # leaves some NaN
         maxGrpSize=max(socGrpSize, popGrpSize, na.rm = TRUE)) # (Ellen) added this one as another meaningful way to consider group size--for double-check. Leaves some -Inf

# --- join and clean host trait data with GMPD data --- 

pan <- right_join(pan, tmp, by="hostName")

dat_combine <- right_join(pan, gmpd_iucn, by = "hostName")

dat_combine %<>% select(-popGrpSize.x, -popGrpSize.y, -socGrpSize.x, -socGrpSize.y)
dat_combine$maxGrpSize[is.infinite(dat_combine$maxGrpSize)] <- NA # (Ellen) just for consistent characterization of non-numbers
dat_combine$comGrpSize[is.nan(dat_combine$comGrpSize)] <- NA
# --- 

write_csv(dat_combine, "./Data/JPEK/script2.csv")

