# second script
# add data for social group from sonia's collaborators for carnivores

library(tidyverse)
library(magrittr)
library(ggplot2) 

dat <- read_csv("Data/JPEK/script2.csv")
carnDat <- read_csv("Data/JPEK/carnivoretraitdata.csv") 

carnDat %<>%
  rename(hostName=`MSW Bionomial`) %>%
  select(hostName, Groupsize) %>%
  mutate(hostName=gsub(pattern=" ", x=hostName, replacement = "_"))

# new column to manipulate 
carnDat$carnGrp <- carnDat$Groupsize
carnDat$carnGrp[carnDat$carnGrp == "Groups"|carnDat$carnGrp == "Pairs to family clans"|carnDat$carnGrp == "Solitary or groups"] <- "Group"
carnDat$carnGrp[carnDat$carnGrp!="Group"] <- "NonGroup"

carnDat %<>% filter(!(is.na(carnGrp))) %>% select(-Groupsize)

# join with rest of data
foo <- right_join(carnDat,dat, by="hostName")

write_csv(foo, "./Data/JPEK/script3.csv")

View(foo %>% filter(Group=="carnivores", is.na(carnGrp)) %>%
       distinct(hostName, .keep_all = TRUE)) 
