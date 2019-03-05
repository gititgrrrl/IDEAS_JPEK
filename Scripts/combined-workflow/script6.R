# sixth script in workflow 
# "FULL MODEL": includes many variables of interest  

# if running code first time: 
# rm(list=ls())
# 
# ### LOAD PACKAGES ----
# packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "rstan", "brms", "broom", "tidybayes", "purrr", "glmmTMB")
# package.check <- lapply(packages, FUN = function(x) {
#   if (!require(x, character.only = TRUE)) {
#     install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
#     library(x, character.only = TRUE)
#   }
# })
# rstan_options (auto_write=TRUE)
# options (mc.cores=parallel::detectCores ())

# if running code on familiar computer: 
library(magrittr)
library(cowplot)
library(GGally)
library(scales)
library(tidyverse)
library(rstan)
library(brms)
library(broom)
library(tidybayes)
library(purrr)
library(glmmTMB)

# --- load cleaned data ---

rm(list=ls())
allDat <- read_csv("Data/JPEK/script4.csv") 

# --- subset data for full model ALL PARASITE TYPES --- 

fullDat <- allDat %>%
  select(hostName, parRich, propClose, numHostCitations, combIUCN, 
         hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, hostMeanLat)  
fullDat <- fullDat[complete.cases(fullDat),] # only 255 records b/c 120 missing species range ans 149 missing max lifespan. With host group size included, only 156 records left.

# --- run full model ALL PARASITE TYPES---

# **** BEST MODEL: MOD4 - ADD INTERACTION BTWN CITATIONS AND IUCN STATUS ----
# PP-check shows still slight problems with ungulates and a bit of problem with carnivores & non-threatened, but not bad
# Loo-pit shows S-shape

fullBrm <- brm(
  data = fullDat, 
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 log(numHostCitations) * combIUCN + # <<<<<<
                 logHostSpeciesRange + 
                 logHostMass * hostGroup +
                 hostMaxLifespan +
                 propClose*hostGroup), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

fullBrm <- add_ic(fullBrm, ic = "loo", reloo = TRUE)  
fullBrm <- add_ic(fullBrm, ic = "kfold")

# --- subset data for full model CLOSE PARASITES --- 

# --- run full model CLOSE PARASITES ---

# --- subset data for full model NON CLOSE PARASITES --- 

# --- run full model NON CLOSE PARASITES---

# --- subset data for full model MICROSPARASITES --- 

# --- run full model MICROSPARASITES---

# --- subset data for full model MACROSPARASITES --- 

# --- run full model MACROSPARASITES---

# --- 

saveRDS(fullBrm, "./Data/JPEK/full_brm_all.RDS") # model of parasite spp richness total
saveRDS(fullBrm_close, "./Data/JPEK/full_brm_close.RDS") #   model of closely transmitted parasite spp richness 
saveRDS(fullBrm_nonclose, "./Data/JPEK/full_brm_nonclose.RDS") #   model of NON closely transmitted parasite spp richness 
saveRDS(fullBrm_micro, "./Data/JPEK/full_brm_micro.RDS") #   model of micro spp richness 
saveRDS(fullBrm_macro, "./Data/JPEK/full_brm_macro.RDS") #   model of macro spp richness
