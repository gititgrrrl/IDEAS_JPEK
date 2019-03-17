### PARASITE RICHNESS SCRIPT FOR MODELS INCLUDING GROUP SIZE

rm(list=ls())

### LOAD PACKAGES ----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "rstan", 
              "brms", "broom", "tidybayes", "purrr", "glmmTMB")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run chains on multiple cores

# # if running code on familiar computer: 
# library(magrittr)
# library(cowplot)
# library(GGally)
# library(scales)
# library(tidyverse)
# library(rstan)
# library(brms)
# library(broom)
# library(tidybayes)
# library(purrr)
# library(glmmTMB)
#
### CARNIVORES WITH PARASITE RICHNESS AS RESPONSE ----
# # Create data
# allDat <- read_csv("Data/JPEK/script4.csv")
# 
# fullDat_carngroup <- allDat %>%
#   filter(hostGroup == "carnivores") %>%
#   select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() # 148 records
# fullDat_carngroup <- fullDat_carngroup[complete.cases(fullDat_carngroup),] # only 74 records
# 
# # ...model with ALL threat interactions
# # > ONLY GOING AS FAR AS IC CALCS ON THESE--IF THEY ARE BETTER MODELS, THEN FOR FINAL PAPER NEED TO GET MARGINAL EFFECTS, ETC.
# allIntBrm_carngroup <- brm(
#   data = fullDat_carngroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations +
#                  combIUCN*logHostSpeciesRange +
#                  combIUCN*groupSizeCar +
#                  combIUCN*logHostMass +
#                  combIUCN*hostMaxLifespan +
#                  combIUCN*absHostMeanLat),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(allIntBrm_carngroup)
# plot(allIntBrm_carngroup)
# pp_check(allIntBrm_carngroup, nsamples = 500)
# 
# # add information criteria
# allIntBrm_carngroup <- add_ic(allIntBrm_carngroup, ic = "loo", reloo = TRUE)
# allIntBrm_carngroup <- add_ic(allIntBrm_carngroup, ic = "kfold")
# saveRDS(allIntBrm_carngroup, "./Data/JPEK/allInt/allInt_brm_carngroup_all.RDS")
#  
# # ...full model (some threat interactions omitted)
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
# 
# fullBrm_carngroup <- readRDS("./Data/JPEK/full/full_brm_carngroup_all.RDS") # model
# # add information criteria
# fullBrm_carngroup <- add_ic(fullBrm_carngroup, ic = "loo", reloo = TRUE)
# fullBrm_carngroup <- add_ic(fullBrm_carngroup, ic = "kfold")
# saveRDS(fullBrm_carngroup, "./Data/JPEK/full/full_brm_carngroup_all.RDS")
# 
# # model fits, predictions and residuals
# fullMu_carngroup <- fitted(fullBrm_carngroup)
# fullPredict_carngroup <- predict(fullBrm_carngroup)
# fullResid_carngroup <- residuals(fullBrm_carngroup, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_carngroup, "./Data/JPEK/full/full_brm_carngroup_all_mu.RDS")
# saveRDS(fullPredict_carngroup, "./Data/JPEK/full/full_brm_carngroup_all_predict.RDS")
# saveRDS(fullResid_carngroup, "./Data/JPEK/full/full_brm_carngroup_all_resid.RDS")
# 
# # marginal effects
# full_carngroup_me <- plot(marginal_effects(fullBrm_carngroup), method = "fitted", plot = FALSE)
# saveRDS(full_carngroup_me, "./Data/JPEK/full/full_brm_carngroup_me.RDS")

### UNGULATES WITH PARASITE RICHNESS AS RESPONSE ----
# # Create data
# fullDat_unggroup <- allDat %>%
#   filter(hostGroup == "ungulates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 102 records
# fullDat_unggroup <- fullDat_unggroup[complete.cases(fullDat_unggroup),] # only 60 records
# 
# # # ...model with ALL threat interactions
# # > ONLY GOING AS FAR AS IC CALCS ON THESE--IF THEY ARE BETTER MODELS, THEN FOR FINAL PAPER NEED TO GET MARGINAL EFFECTS, ETC.
# allIntBrm_unggroup <- brm(
#   data = fullDat_unggroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations +
#                  combIUCN*logHostSpeciesRange +
#                  combIUCN*logGroupSizePriUng +
#                  combIUCN*logHostMass +
#                  combIUCN*hostMaxLifespan +
#                  combIUCN*absHostMeanLat),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(allIntBrm_unggroup)
# plot(allIntBrm_unggroup)
# pp_check(allIntBrm_unggroup, nsamples = 500)
# 
# # add information criteria
# allIntBrm_unggroup <- add_ic(allIntBrm_unggroup, ic = "loo", reloo = TRUE)
# allIntBrm_unggroup <- add_ic(allIntBrm_unggroup, ic = "kfold") # didn't work with default K = 10 or with K = 5
# saveRDS(allIntBrm_unggroup, "./Data/JPEK/allInt/allInt_brm_unggroup_all.RDS")
# 
# # ...full model (some threat interactions omitted)
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
# fullBrm_unggroup <- readRDS("./Data/JPEK/full/full_brm_unggroup_all.RDS") # model
# # add information criteria
# fullBrm_unggroup <- add_ic(fullBrm_unggroup, ic = "loo", reloo = TRUE)
# fullBrm_unggroup <- add_ic(fullBrm_unggroup, ic = "kfold")
# saveRDS(fullBrm_unggroup, "./Data/JPEK/full/full_brm_unggroup_all.RDS")
# 
# # model fits, predictions and residuals
# fullMu_unggroup <- fitted(fullBrm_unggroup)
# fullPredict_unggroup <- predict(fullBrm_unggroup)
# fullResid_unggroup <- residuals(fullBrm_unggroup, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_unggroup, "./Data/JPEK/full/full_brm_unggroup_all_mu.RDS")
# saveRDS(fullPredict_unggroup, "./Data/JPEK/full/full_brm_unggroup_all_predict.RDS")
# saveRDS(fullResid_unggroup, "./Data/JPEK/full/full_brm_unggroup_all_resid.RDS")
# 
# # marginal effects
# full_unggroup_me <- plot(marginal_effects(fullBrm_unggroup), method = "fitted", plot = FALSE)
# saveRDS(full_unggroup_me, "./Data/JPEK/full/full_brm_unggroup_me.RDS")
# 
### PRIMATES WITH PARASITE RICHNESS AS RESPONSE ----
# # Create data
# fullDat_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 142 records
# fullDat_primgroup <- fullDat_primgroup[complete.cases(fullDat_primgroup),] # only 73 records
# 
# # ...model with ALL threat interactions
# > ONLY GOING AS FAR AS IC CALCS ON THESE--IF THEY ARE BETTER MODELS, THEN FOR FINAL PAPER NEED TO GET MARGINAL EFFECTS, ETC.
# allIntBrm_primgroup <- brm(
#   data = fullDat_primgroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations +
#                  combIUCN*logHostSpeciesRange +
#                  combIUCN*logGroupSizePriUng +
#                  combIUCN*logHostMass +
#                  combIUCN*hostMaxLifespan +
#                  combIUCN*absHostMeanLat),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(allIntBrm_primgroup)
# plot(allIntBrm_primgroup)
# pp_check(allIntBrm_primgroup, nsamples = 500)
# 
# # add information criteria
# allIntBrm_primgroup <- add_ic(allIntBrm_primgroup, ic = "loo", reloo = TRUE)
# allIntBrm_primgroup <- add_ic(allIntBrm_primgroup, ic = "kfold")
# saveRDS(allIntBrm_primgroup, "./Data/JPEK/allInt/allInt_brm_primgroup_all.RDS")
# 
# # ...full model (some threat interactions omitted)
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
# 
# fullBrm_primgroup <- readRDS("./Data/JPEK/full/full_brm_primgroup_all.RDS") # model
# # add information criteria
# fullBrm_primgroup <- add_ic(fullBrm_primgroup, ic = "loo", reloo = TRUE)
# fullBrm_primgroup <- add_ic(fullBrm_primgroup, ic = "kfold")
# saveRDS(fullBrm_primgroup, "./Data/JPEK/full/full_brm_primgroup_all.RDS")
# 
# # model fits, predictions and residuals
# fullMu_primgroup <- fitted(fullBrm_primgroup)
# fullPredict_primgroup <- predict(fullBrm_primgroup)
# fullResid_primgroup <- residuals(fullBrm_primgroup, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_primgroup, "./Data/JPEK/full/full_brm_primgroup_all_mu.RDS")
# saveRDS(fullPredict_primgroup, "./Data/JPEK/full/full_brm_primgroup_all_predict.RDS")
# saveRDS(fullResid_primgroup, "./Data/JPEK/full/full_brm_primgroup_all_resid.RDS")
# 
# # marginal effects
# full_primgroup_me <- plot(marginal_effects(fullBrm_primgroup), method = "fitted", plot = FALSE)
# saveRDS(full_primgroup_me, "./Data/JPEK/full/full_brm_primgroup_me.RDS")
 
 
 
### SIMPLE MODELS ARE DONE #### ----
### SIMPLE MODEL, CARNIVORES WITH PARASITE RICHNESS AS RESPONSE ----
# simpBrm_carngroup <- brm(
#   data = fullDat_carngroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_carngroup)
# plot(simpBrm_carngroup)
# pp_check(simpBrm_carngroup, nsamples = 500)
# saveRDS(simpBrm_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all.RDS") # model
# 
# # simpBrm_carngroup <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all.RDS") # model
# # add information criteria
# simpBrm_carngroup <- add_ic(simpBrm_carngroup, ic = "loo", reloo = TRUE)
# simpBrm_carngroup <- add_ic(simpBrm_carngroup, ic = "kfold")
# saveRDS(simpBrm_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all.RDS")
#
# # model fits, predictions and residuals
# simpMu_carngroup <- fitted(simpBrm_carngroup)
# simpPredict_carngroup <- predict(simpBrm_carngroup)
# simpResid_carngroup <- residuals(simpBrm_carngroup, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all_mu.RDS")
# saveRDS(simpPredict_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all_predict.RDS")
# saveRDS(simpResid_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all_resid.RDS")
# 
# # marginal effects
# simp_carngroup_me <- plot(marginal_effects(simpBrm_carngroup), method = "fitted", plot = FALSE)
# saveRDS(simp_carngroup_me, "./Data/JPEK/simple/simp_brm_carngroup_me.RDS")
#
### SIMPLE MODEL, UNGULATES WITH PARASITE RICHNESS AS RESPONSE ----
# simpBrm_unggroup <- brm(
#   data = fullDat_unggroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_unggroup)
# plot(simpBrm_unggroup)
# pp_check(simpBrm_unggroup, nsamples = 500)
# saveRDS(simpBrm_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all.RDS") # model
# 
# simpBrm_unggroup <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all.RDS") # model
# # add information criteria
# simpBrm_unggroup <- add_ic(simpBrm_unggroup, ic = "loo", reloo = TRUE)
# simpBrm_unggroup <- add_ic(simpBrm_unggroup, ic = "kfold")
# saveRDS(simpBrm_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all.RDS")
# 
# # model fits, predictions and residuals
# simpMu_unggroup <- fitted(simpBrm_unggroup)
# simpPredict_unggroup <- predict(simpBrm_unggroup)
# simpResid_unggroup <- residuals(simpBrm_unggroup, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all_mu.RDS")
# saveRDS(simpPredict_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all_predict.RDS")
# saveRDS(simpResid_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all_resid.RDS")
# 
# # marginal effects
# simp_unggroup_me <- plot(marginal_effects(simpBrm_unggroup), method = "fitted", plot = FALSE)
# saveRDS(simp_unggroup_me, "./Data/JPEK/simple/simp_brm_unggroup_me.RDS")
#
### SIMPLE MODEL, PRIMATES WITH PARASITE RICHNESS AS RESPONSE ----
# # fullDat_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 142 records
# fullDat_primgroup <- fullDat_primgroup[complete.cases(fullDat_primgroup),] # only 73 records
# 
# simpBrm_primgroup <- brm(
#   data = fullDat_primgroup,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_primgroup)
# plot(simpBrm_primgroup)
# pp_check(simpBrm_primgroup, nsamples = 500)
# saveRDS(simpBrm_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all.RDS") # model
# 
# simpBrm_primgroup <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all.RDS") # model
# # add information criteria
# simpBrm_primgroup <- add_ic(simpBrm_primgroup, ic = "loo", reloo = TRUE)
# simpBrm_primgroup <- add_ic(simpBrm_primgroup, ic = "kfold")
# saveRDS(simpBrm_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all.RDS") # model
# 
# # model fits, predictions and residuals
# simpMu_primgroup <- fitted(simpBrm_primgroup)
# simpPredict_primgroup <- predict(simpBrm_primgroup)
# simpResid_primgroup <- residuals(simpBrm_primgroup, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all_mu.RDS")
# saveRDS(simpPredict_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all_predict.RDS")
# saveRDS(simpResid_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all_resid.RDS")
# 
# # marginal effects
# simp_primgroup_me <- plot(marginal_effects(simpBrm_primgroup), method = "fitted", plot = FALSE)
# saveRDS(simp_primgroup_me, "./Data/JPEK/simple/simp_brm_primgroup_me.RDS")