# Additional code for Paige to run on fast computers

# Get residuals for full and simple_full data models
# RICHNESS
simple_brm_all_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_fulldat.RDS") 
full_brm_all <- readRDS("./Data/JPEK/full/full_brm_all.RDS") 

simple_brm_all_resid_fulldat <- residuals(simple_brm_all_fulldat, type = "pearson")
full_brm_all_resid <- residuals(full_brm_all, type = "pearson")

saveRDS(simple_brm_all_resid_fulldat, "./Data/JPEK/simple/simple_brm_all_resid_fulldat.RDS")
saveRDS(full_brm_all_resid, "./Data/JPEK/full/full_brm_all_resid.RDS")

# PARAS TRANS
simple_brm_parastrans_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastrans_fulldat.RDS") 
full_brm_parastrans <- readRDS("./Data/JPEK/full/full_brm_parastrans.RDS") 

simple_brm_parastrans_resid_fulldat <- residuals(simple_brm_parastrans_fulldat, type = "pearson")
full_brm_parastrans_resid <- residuals(full_brm_parastrans, type = "pearson")

saveRDS(simple_brm_parastrans_resid_fulldat, "./Data/JPEK/simple/simple_brm_parastrans_resid_fulldat.RDS")
saveRDS(full_brm_parastrans_resid, "./Data/JPEK/full/full_brm_parastrans_resid.RDS")

# PARAS TYPE
simple_brm_parastype_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastype_fulldat.RDS") 
full_brm_parastype <- readRDS("./Data/JPEK/full/full_brm_parastype.RDS") 

simple_brm_parastype_resid_fulldat <- residuals(simple_brm_parastype_fulldat, type = "pearson")
full_brm_parastype_resid <- residuals(full_brm_parastype, type = "pearson")

saveRDS(simple_brm_parastype_resid_fulldat, "./Data/JPEK/simple/simple_brm_parastype_resid_fulldat.RDS")
saveRDS(full_brm_parastype_resid, "./Data/JPEK/full/full_brm_parastype_resid.RDS")

### -- GROUP SIZE SUB-ANALYSES ----
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

# allDat <- read_csv("Data/JPEK/script4.csv")
# # --- subset data for full model ALL PARASITES for CARNIVORES ---- 
# # (ELLEN ALREADY RAN THIS MODEL AND SAVED THE MODEL TO FILES)
# # It seems solitary animals have higher parasite richness
# 
# fullDat_carngroup <- allDat %>%
#   filter(hostGroup == "carnivores") %>%
#   select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() # 148 records
# fullDat_carngroup <- fullDat_carngroup[complete.cases(fullDat_carngroup),] # only 74 records
# 
# 
# fullBrm_carngroup <- brm(
#   data = fullDat_carngroup, 
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations +
#                  logHostSpeciesRange +
#                  combIUCN:logHostSpeciesRange +
#                  groupSizeCar +
#                  combIUCN:groupSizeCar +
#                  logHostMass +
#                  hostMaxLifespan +
#                  absHostMeanLat), 
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(fullBrm_carngroup)
# plot(fullBrm_carngroup)
# pp_check(fullBrm_carngroup, nsamples = 500)
# saveRDS(fullBrm_carngroup, "./Data/JPEK/full/full_brm_carngroup_all.RDS") # model

# FOR PAIGE TO RUN ON SUPER COMPUTER
fullBrm_carngroup <- readRDS("./Data/JPEK/full/full_brm_carngroup_all.RDS") # model
# add information criteria
fullBrm_carngroup <- add_ic(fullBrm_carngroup, ic = "loo", reloo = TRUE)
fullBrm_carngroup <- add_ic(fullBrm_carngroup, ic = "kfold")

# model fits, predictions and residuals
fullMu_carngroup <- fitted(fullBrm_carngroup)
fullPredict_carngroup <- predict(fullBrm_carngroup)
fullResid_carngroup <- residuals(fullBrm_carngroup, type = "pearson")

# save outputs
saveRDS(fullMu_carngroup, "./Data/JPEK/full/full_brm_carngroup_all_mu.RDS")
saveRDS(fullPredict_carngroup, "./Data/JPEK/full/full_brm_carngroup_all_predict.RDS")
saveRDS(fullResid_carngroup, "./Data/JPEK/full/full_brm_carngroup_all_resid.RDS")

# marginal effects
full_carngroup_me <- plot(marginal_effects(fullBrm_carngroup), method = "fitted", plot = FALSE)
saveRDS(full_carngroup_me, "./Data/JPEK/full/full_brm_carngroup_me.RDS")

# # --- subset data for full model ALL PARASITES for UNGULATES ----
# # (ELLEN ALREADY RAN THIS MODEL AND SAVED THE MODEL TO FILES)
# # It seems parasite richness declines as group size increases
# fullDat_unggroup <- allDat %>%
#   filter(hostGroup == "ungulates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 102 records
# fullDat_unggroup <- fullDat_unggroup[complete.cases(fullDat_unggroup),] # only 60 records
# 
# fullBrm_unggroup <- brm(
#   data = fullDat_unggroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations +
#                  logHostSpeciesRange +
#                  combIUCN:logHostSpeciesRange +
#                  logGroupSizePriUng +
#                  combIUCN:logGroupSizePriUng +
#                  logHostMass +
#                  hostMaxLifespan +
#                  absHostMeanLat),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(fullBrm_unggroup)
# plot(fullBrm_unggroup)
# pp_check(fullBrm_unggroup, nsamples = 500)
# saveRDS(fullBrm_unggroup, "./Data/JPEK/full/full_brm_unggroup_all.RDS") # model

# FOR PAIGE TO RUN ON SUPER COMPUTER
fullBrm_unggroup <- readRDS("./Data/JPEK/full/full_brm_unggroup_all.RDS") # model
# add information criteria
fullBrm_unggroup <- add_ic(fullBrm_unggroup, ic = "loo", reloo = TRUE)
fullBrm_unggroup <- add_ic(fullBrm_unggroup, ic = "kfold")

# model fits, predictions and residuals
fullMu_unggroup <- fitted(fullBrm_unggroup)
fullPredict_unggroup <- predict(fullBrm_unggroup)
fullResid_unggroup <- residuals(fullBrm_unggroup, type = "pearson")

# save outputs
saveRDS(fullMu_unggroup, "./Data/JPEK/full/full_brm_unggroup_all_mu.RDS")
saveRDS(fullPredict_unggroup, "./Data/JPEK/full/full_brm_unggroup_all_predict.RDS")
saveRDS(fullResid_unggroup, "./Data/JPEK/full/full_brm_unggroup_all_resid.RDS")

# marginal effects
full_unggroup_me <- plot(marginal_effects(fullBrm_unggroup), method = "fitted", plot = FALSE)
saveRDS(full_unggroup_me, "./Data/JPEK/full/full_brm_unggroup_me.RDS")

# # --- subset data for full model ALL PARASITES for PRIMATES ----
# # (ELLEN ALREADY RAN THIS MODEL AND SAVED THE MODEL TO FILES)
# # It seems group size is not important for primates
# fullDat_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 142 records
# fullDat_primgroup <- fullDat_primgroup[complete.cases(fullDat_primgroup),] # only 73 records
# 
# fullBrm_primgroup <- brm(
#   data = fullDat_primgroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations +
#                  logHostSpeciesRange +
#                  combIUCN:logHostSpeciesRange +
#                  logGroupSizePriUng +
#                  combIUCN:logGroupSizePriUng +
#                  logHostMass +
#                  hostMaxLifespan +
#                  absHostMeanLat),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(fullBrm_primgroup)
# plot(fullBrm_primgroup)
# pp_check(fullBrm_primgroup, nsamples = 500)
# saveRDS(fullBrm_primgroup, "./Data/JPEK/full/full_brm_primgroup_all.RDS") # model

# FOR PAIGE TO RUN ON SUPER COMPUTER
fullBrm_primgroup <- readRDS("./Data/JPEK/full/full_brm_primgroup_all.RDS") # model
# add information criteria
fullBrm_primgroup <- add_ic(fullBrm_primgroup, ic = "loo", reloo = TRUE)
fullBrm_primgroup <- add_ic(fullBrm_primgroup, ic = "kfold")

# model fits, predictions and residuals
fullMu_primgroup <- fitted(fullBrm_primgroup)
fullPredict_primgroup <- predict(fullBrm_primgroup)
fullResid_primgroup <- residuals(fullBrm_primgroup, type = "pearson")

# save outputs
saveRDS(fullMu_primgroup, "./Data/JPEK/full/full_brm_primgroup_all_mu.RDS")
saveRDS(fullPredict_primgroup, "./Data/JPEK/full/full_brm_primgroup_all_predict.RDS")
saveRDS(fullResid_primgroup, "./Data/JPEK/full/full_brm_primgroup_all_resid.RDS")

# marginal effects
full_primgroup_me <- plot(marginal_effects(fullBrm_primgroup), method = "fitted", plot = FALSE)
saveRDS(full_primgroup_me, "./Data/JPEK/full/full_brm_primgroup_me.RDS")
