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

# --- subset data for full model ALL PARASITES --- 

fullDat <- allDat %>%
  select(hostName, parRich, logNumHostCitations, combIUCN, 
         hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() # 392 records
fullDat <- fullDat[complete.cases(fullDat),] # only 232 records (before, it was 255 records--not sure why it's fewer now <<<<< ?????)

# quick check of how correlated the predictors are--no very strong correlations, so should be fine to use all
ggpairs(fullDat[, 3:9])

### --- run full model ALL PARASITES---
# Question is, which interactions should be included in full model?
fullBrm <- brm(
  data = fullDat, 
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN * hostGroup +
                 logNumHostCitations +
                 combIUCN:logNumHostCitations +
                 logHostSpeciesRange +
                 combIUCN:logHostSpeciesRange +
                 hostGroup:logHostSpeciesRange + 
                 combIUCN:hostGroup:logHostSpeciesRange +
                 logHostMass +
                 hostGroup:logHostMass +
                 hostMaxLifespan +
                 hostGroup:hostMaxLifespan +
                 absHostMeanLat +
                 hostGroup:absHostMeanLat), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(fullBrm)
plot(fullBrm)
pp_check(fullBrm, nsamples = 500)

# add information criteria
# fullBrm <- add_ic(fullBrm, ic = "loo", reloo = TRUE)  
# fullBrm <- add_ic(fullBrm, ic = "kfold")

# model fits and predictions
fullMu <- fitted(fullBrm)
fullPredict <- predict(fullBrm)

# save outputs
saveRDS(fullBrm, "./Data/JPEK/full/full_brm_all.RDS") # model of parasite spp richness total
saveRDS(fullMu, "./Data/JPEK/full/full_brm_all_mu.RDS")
saveRDS(fullPredict, "./Data/JPEK/full/full_brm_all_predict.RDS")

# marginal effects
full_me <- plot(marginal_effects(fullBrm), method = "fitted", plot = FALSE)
saveRDS(full_me, "./Data/JPEK/full/full_brm_me.RDS")

# --- run SIMPLE model ALL PARASITES on the same subset of data---
# doesn't make sense to call it fulldat b/c it's actually a smaller subset of the data, but the 'full' refers to the data used in the full model
simpleBrm_fulldat <- brm(
  data = fullDat, 
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 combIUCN * logNumHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# quick checks
summary(simpleBrm_fulldat)
plot(simpleBrm_fulldat)
pp_check(simpleBrm_fulldat, nsamples = 500)

# add information criteria
# simpleBrm_fulldat <- add_ic(simpleBrm_fulldat, ic = "loo", reloo = TRUE)
# simpleBrm_fulldat <- add_ic(simpleBrm_fulldat, ic = "kfold")

# model fits and predictions
simpleMu_fulldat <- fitted(simpleBrm_fulldat)
simplePredict_fulldat <- predict(simpleBrm_fulldat)

# save outputs
saveRDS(simpleBrm_fulldat, "./Data/JPEK/simple/simple_brm_all_fulldat.RDS") # model of parasite spp richness total
saveRDS(simpleMu_fulldat, "./Data/JPEK/simple/simple_brm_all_mu_fulldat.RDS")
saveRDS(simplePredict_fulldat, "./Data/JPEK/simple/simple_brm_all_predict_fulldat.RDS")

# marginal effects
simple_me_fulldat <- plot(marginal_effects(simpleBrm_fulldat), method = "fitted", plot = FALSE)
saveRDS(simple_me_fulldat, "./Data/JPEK/simple/simple_brm_me_fulldat.RDS")



### --- subset data for full model PARASITE TRANSMISSION --- 
fullDat_parastrans <- allDat %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, logNumHostCitations, combIUCN, hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 362

fullDat_parastrans <- fullDat_parastrans[complete.cases(fullDat_parastrans),]   # 220

# --- run full model PARASITE TRANSMISSION  ---

# BINOMIAL regression
fullBrm_parastrans <- brm(
  data = fullDat_parastrans,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN * hostGroup +
                 logNumHostCitations +
                 combIUCN:logNumHostCitations +
                 logHostSpeciesRange +
                 combIUCN:logHostSpeciesRange +
                 hostGroup:logHostSpeciesRange + 
                 combIUCN:hostGroup:logHostSpeciesRange +
                 logHostMass +
                 hostGroup:logHostMass +
                 hostMaxLifespan +
                 hostGroup:hostMaxLifespan +
                 absHostMeanLat +
                 hostGroup:absHostMeanLat), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# quick checks
summary(fullBrm_parastrans)
plot(fullBrm_parastrans)
pp_check(fullBrm_parastrans, nsamples = 500)

# add information criteria
# fullBrm_parastrans <- add_ic(fullBrm_parastrans, ic = "loo", reloo = TRUE)  
# fullBrm_parastrans <- add_ic(fullBrm_parastrans, ic = "kfold")

# model fits and predictions
fullMu_parastrans <- fitted(fullBrm_parastrans)
fullPredict_parastrans <- predict(fullBrm_parastrans)

# save outputs
saveRDS(fullBrm_parastrans, "./Data/JPEK/full/full_brm_parastrans.RDS") # model of parasite spp richness total
saveRDS(fullMu_parastrans, "./Data/JPEK/full/full_brm_parastrans_mu.RDS")
saveRDS(fullPredict_parastrans, "./Data/JPEK/full/full_brm_parastrans_predict.RDS")

# marginal effects
full_me_parastrans <- plot(marginal_effects(fullBrm_parastrans), method = "fitted", plot = FALSE)
saveRDS(full_me_parastrans, "./Data/JPEK/full/full_brm_parastrans_me.RDS")


# --- run SIMPLE model PARASITE TRANSMISSION on the same subset of data---
simpleBrm_parastrans_fulldat <- brm(
  data = fullDat_parastrans, 
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN * hostGroup + 
                 combIUCN * logNumHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  

# quick checks
summary(simpleBrm_parastrans_fulldat)
plot(simpleBrm_parastrans_fulldat)
pp_check(simpleBrm_parastrans_fulldat, nsamples = 500)

# add information criteria
simpleBrm_parastrans_fulldat <- add_ic(simpleBrm_parastrans_fulldat, ic = "loo", reloo = TRUE)
simpleBrm_parastrans_fulldat <- add_ic(simpleBrm_parastrans_fulldat, ic = "kfold")

# model fits and predictions
simpleMu_parastrans_fulldat <- fitted(simpleBrm_parastrans_fulldat)
simplePredict_parastrans_fulldat <- predict(simpleBrm_parastrans_fulldat)

# save outputs
saveRDS(simpleBrm_parastrans_fulldat, "./Data/JPEK/simple/simple_brm_parastrans_fulldat.RDS") # model of parasite spp by transmission
saveRDS(simpleMu_parastrans_fulldat, "./Data/JPEK/simple/simple_brm_parastrans_mu_fulldat.RDS")
saveRDS(simplePredict_parastrans_fulldat, "./Data/JPEK/simple/simple_brm_parastrans_predict_fulldat.RDS")

# marginal effects
simple_me_parastrans_fulldat <- plot(marginal_effects(simpleBrm_parastrans_fulldat), method = "fitted", plot = FALSE)
saveRDS(simple_me_parastrans_fulldat, "./Data/JPEK/simple/simple_brm_parastrans_me_fulldat.RDS")



### --- subset data for full model PARASITE TYPE --- 
fullDat_parastype <- allDat %>%
  select(hostName, logNumHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 393

fullDat_parastype <- fullDat_parastype[complete.cases(fullDat_parastype),] # 232

# --- run full model PARASITE TYPE  ---

# BINOMIAL regression
fullBrm_parastype <- brm(
  data = fullDat_parastype,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN * hostGroup +
                 logNumHostCitations +
                 combIUCN:logNumHostCitations +
                 logHostSpeciesRange +
                 combIUCN:logHostSpeciesRange +
                 hostGroup:logHostSpeciesRange + 
                 combIUCN:hostGroup:logHostSpeciesRange +
                 logHostMass +
                 hostGroup:logHostMass +
                 hostMaxLifespan +
                 hostGroup:hostMaxLifespan +
                 absHostMeanLat +
                 hostGroup:absHostMeanLat), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# quick checks
summary(fullBrm_parastype)
plot(fullBrm_parastype)
pp_check(fullBrm_parastype, nsamples = 500)

# add information criteria
fullBrm_parastype <- add_ic(fullBrm_parastype, ic = "loo", reloo = TRUE)  
fullBrm_parastype <- add_ic(fullBrm_parastype, ic = "kfold")

# model fits and predictions
fullMu_parastype <- fitted(fullBrm_parastype)
fullPredict_parastype <- predict(fullBrm_parastype)

# save outputs
saveRDS(fullBrm_parastype, "./Data/JPEK/full/full_brm_parastype.RDS") # model of parasite by type
saveRDS(fullMu_parastype, "./Data/JPEK/full/full_brm_parastype_mu.RDS")
saveRDS(fullPredict_parastype, "./Data/JPEK/full/full_brm_parastype_predict.RDS")

# marginal effects
full_me_parastype <- plot(marginal_effects(fullBrm_parastype), method = "fitted", plot = FALSE)
saveRDS(full_me_parastype, "./Data/JPEK/full/full_brm_parastype_me.RDS")


# --- run SIMPLE model PARASITE TYPE on the same subset of data---
simpleBrm_parastype_fulldat <- brm(
  data = fullDat_parastype,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN * hostGroup + 
                 combIUCN * logNumHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  

# quick checks
summary(simpleBrm_parastype_fulldat)
plot(simpleBrm_parastype_fulldat)
pp_check(simpleBrm_parastype_fulldat, nsamples = 500)

# add information criteria
simpleBrm_parastype_fulldat <- add_ic(simpleBrm_parastype_fulldat, ic = "loo", reloo = TRUE)
simpleBrm_parastype_fulldat <- add_ic(simpleBrm_parastype_fulldat, ic = "kfold")

# model fits and predictions
simpleMu_parastype_fulldat <- fitted(simpleBrm_parastype_fulldat)
simplePredict_parastype_fulldat <- predict(simpleBrm_parastype_fulldat)

# save outputs
saveRDS(simpleBrm_parastype_fulldat, "./Data/JPEK/simple/simple_brm_parastype_fulldat.RDS") # model of parasite spp by type
saveRDS(simpleMu_parastype_fulldat, "./Data/JPEK/simple/simple_brm_parastype_mu_fulldat.RDS")
saveRDS(simplePredict_parastype_fulldat, "./Data/JPEK/simple/simple_brm_parastype_predict_fulldat.RDS")

# marginal effects
simple_me_parastype_fulldat <- plot(marginal_effects(simpleBrm_parastype_fulldat), method = "fitted", plot = FALSE)
saveRDS(simple_me_parastype_fulldat, "./Data/JPEK/simple/simple_brm_parastype_me_fulldat.RDS")