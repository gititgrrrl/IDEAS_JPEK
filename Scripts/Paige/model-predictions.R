# Model observations and predictions

library(magrittr)
library(cowplot)
library(tidyverse)
library(rstan)
library(brms)
library(broom)
library(tidybayes)
library(purrr)

### --- DATA --- ###
allDat <- read_csv("./Data/JPEK/script4.csv") # columns 28 & 29 = # host and parasite citations 

# SIMPLE RICHNESS
simpleDat <- allDat %>%
  select(hostName, parRich, logNumHostCitations, combIUCN, hostGroup) %>% # <<<<< (Ellen) grab the logNumHostCitations (not numHostCitations)
  distinct()
simpleDat <- simpleDat[complete.cases(simpleDat),]              # 362 hosts


# FULL RICHNESS
fullDat <- allDat %>%
  select(hostName, parRich, logNumHostCitations, combIUCN, 
         hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() # 392 records
fullDat <- fullDat[complete.cases(fullDat),] # only 232 records (before, it was 255 records--not sure why it's fewer now <<<<< ?????)

# FULL PARAS TRANS
fullDat_parastrans <- allDat %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, logNumHostCitations, combIUCN, hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 362

fullDat_parastrans <- fullDat_parastrans[complete.cases(fullDat_parastrans),]   # 220

# FULL PARA TYPE
fullDat_parastype <- allDat %>%
  select(hostName, logNumHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 393

fullDat_parastype <- fullDat_parastype[complete.cases(fullDat_parastype),] # 232

### --- MODEL PARAMETERS --- ###

simpleMu <- readRDS("./Data/JPEK/simple/simple_brm_all_mu.RDS")
simpleMu_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_mu_fulldat.RDS")
fullMu <- readRDS("./Data/JPEK/full/full_brm_all_mu.RDS")

### -- MODEL PREDICTIONS --- ###

simplePred <- readRDS("./Data/JPEK/simple/simple_brm_all_predict.RDS")
simplePred_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_predict_fulldat.RDS")
fullPred <- readRDS("./Data/JPEK/full/full_brm_all_predict.RDS")

### --- MODEL OBSERVATIONS PREDICTIONS --- ###

