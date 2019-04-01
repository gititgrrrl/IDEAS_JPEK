# rm(list=ls())
# 
### LOAD PACKAGES ----
# packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "rstan",
#               "brms", "broom", "tidybayes", "purrr", "glmmTMB")
# package.check <- lapply(packages, FUN = function(x) {
#   if (!require(x, character.only = TRUE)) {
#     install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
#     library(x, character.only = TRUE)
#   }
# })
# rstan_options (auto_write=TRUE)
# options (mc.cores=parallel::detectCores ()) # Run chains on multiple cores

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

### FORMAT THE DATA----
# # ...carnivores----
# fullDat_carngroup <- allDat %>%
#   filter(hostGroup == "carnivores") %>%
#   select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass) %>%
#   distinct() # 141 records
# fullDat_carngroup <- fullDat_carngroup[complete.cases(fullDat_carngroup),] 
# dim(fullDat_carngroup) # 77 records, mostly b/c missing group size data
# 
# fullDat_carngroup$groupSizeCar <- factor(fullDat_carngroup$groupSizeCar, levels = c("non_group", "group"))
# fullDat_carngroup$combIUCN <- factor(fullDat_carngroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_carngroup_c <- fullDat_carngroup %>%
#   mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE))
# write_csv(fullDat_carngroup_c, "./Data/JPEK/allDat_totrich_carngroup.csv")
        
# # # Quick check correlation:
# # ggpairs(data = subset(fullDat_carngroup_c, select = -c(hostName, parRich)))
# 
# #...primates----
# fullDat_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 130 records
# fullDat_primgroup <- fullDat_primgroup[complete.cases(fullDat_primgroup),] # only 87 records
# dim(fullDat_primgroup)
# 
# fullDat_primgroup$combIUCN <- factor(fullDat_primgroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_primgroup_c <- fullDat_primgroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE))
# write_csv(fullDat_primgroup_c, "./Data/JPEK/allDat_totrich_primgroup.csv")
# 
# # # Quick check correlation:
# # ggpairs(data = subset(fullDat_primgroup_c, select = -c(hostName, parRich)))
# 
# #...ungulates----
# fullDat_unggroup <- allDat %>%
#   filter(hostGroup == "ungulates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 97 records
# fullDat_unggroup <- fullDat_unggroup[complete.cases(fullDat_unggroup),] # only 61 records
# dim(fullDat_unggroup)
# 
# fullDat_unggroup$combIUCN <- factor(fullDat_unggroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_unggroup_c <- fullDat_unggroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE))
# write_csv(fullDat_unggroup_c, "./Data/JPEK/allDat_totrich_unggroup.csv")
# 
# # Quick check correlation:
# ggpairs(data = subset(fullDat_unggroup_c, select = -c(hostName, parRich)))

# ### PARASITE RICHNESS IS RESPONSE, FULL MODELS----
# # ...carnivores----
# full_totrich_carngroup_mod <- brm(
#   data = fullDat_carngroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*groupSizeCar +
#                  logHostSpeciesRange_c +
#                  logHostMass_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))
# 
# # quick checks
# summary(full_totrich_carngroup_mod)
# plot(full_totrich_carngroup_mod)
# pp_check(full_totrich_carngroup_mod, nsamples = 500)
# 
# # add information criteria
# full_totrich_carngroup_mod <- readRDS("./FINAL/full/full_totrich_carngroup_mod.RDS")
# full_totrich_carngroup_mod <- add_ic(full_totrich_carngroup_mod, ic = "loo", reloo = TRUE)
# full_totrich_carngroup_mod <- add_ic(full_totrich_carngroup_mod, ic = "kfold")
# saveRDS(full_totrich_carngroup_mod, "./FINAL/full/full_totrich_carngroup_mod.RDS")
# 
# # model fits, predictions and residuals
# full_totrich_carngroup_mod <- readRDS("./FINAL/full/full_totrich_carngroup_mod.RDS")
# full_totrich_carngroup_fitted <- fitted(full_totrich_carngroup_mod)
# full_totrich_carngroup_predict <- predict(full_totrich_carngroup_mod)
# full_totrich_carngroup_resid <- residuals(full_totrich_carngroup_mod, type = "pearson")
# 
# # save outputs
# saveRDS(full_totrich_carngroup_fitted, "./FINAL/full/full_totrich_carngroup_fitted.RDS")
# saveRDS(full_totrich_carngroup_predict, "./FINAL/full/full_totrich_carngroup_predict.RDS")
# saveRDS(full_totrich_carngroup_resid, "./FINAL/full/full_totrich_carngroup_resid.RDS")
# 
# # marginal effects
# full_totrich_carngroup_marginal <- marginal_effects(full_totrich_carngroup_mod, method = "fitted")
# saveRDS(full_totrich_carngroup_marginal, "./FINAL/full/full_totrich_carngroup_marginal.RDS")

# specifically for plotting missing marginal effects
full_totrich_carngroup_mod <- readRDS("./FINAL/full/full_totrich_carngroup_mod.RDS")
full_totrich_carngroup_marginal_mass <- marginal_effects(full_totrich_carngroup_mod, effects = "logHostMass_c", method = "fitted")
saveRDS(full_totrich_carngroup_marginal_mass, "./FINAL/full/full_totrich_carngroup_marginal_mass.RDS")

full_totrich_carngroup_marginal_range <- marginal_effects(full_totrich_carngroup_mod, effects = "logHostSpeciesRange_c", method = "fitted")
saveRDS(full_totrich_carngroup_marginal_range, "./FINAL/full/full_totrich_carngroup_marginal_range.RDS")

# # specifically for plotting groupsize-by-threat effects
# full_totrich_carngroup_mod <- readRDS("./FINAL/full/full_totrich_carngroup_mod.RDS")
# full_totrich_carngroup_marginal_groupsize <- marginal_effects(full_totrich_carngroup_mod, effects = "groupSizeCar:combIUCN", method = "fitted")
# saveRDS(full_totrich_carngroup_marginal_groupsize, "./FINAL/full/full_totrich_carngroup_marginal_groupsize.RDS")

# # ...primates
# full_totrich_primgroup_mod <- brm(
#   data = fullDat_primgroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logGroupSizePriUng_c +
#                  logHostSpeciesRange_c +
#                  logHostMass_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(full_totrich_primgroup_mod)
# plot(full_totrich_primgroup_mod)
# pp_check(full_totrich_primgroup_mod, nsamples = 500)
# 
# # add information criteria
# full_totrich_primgroup_mod <- add_ic(full_totrich_primgroup_mod, ic = "loo", reloo = TRUE)
# full_totrich_primgroup_mod <- add_ic(full_totrich_primgroup_mod, ic = "kfold")
# saveRDS(full_totrich_primgroup_mod, "./FINAL/full/full_totrich_primgroup_mod.RDS")

# # model fits, predictions and residuals
# full_totrich_primgroup_mod <- readRDS("./FINAL/full/full_totrich_primgroup_mod.RDS")
# full_totrich_primgroup_fitted <- fitted(full_totrich_primgroup_mod)
# full_totrich_primgroup_predict <- predict(full_totrich_primgroup_mod)
# full_totrich_primgroup_resid <- residuals(full_totrich_primgroup_mod, type = "pearson")
# 
# # save outputs
# saveRDS(full_totrich_primgroup_fitted, "./FINAL/full/full_totrich_primgroup_fitted.RDS")
# saveRDS(full_totrich_primgroup_predict, "./FINAL/full/full_totrich_primgroup_predict.RDS")
# saveRDS(full_totrich_primgroup_resid, "./FINAL/full/full_totrich_primgroup_resid.RDS")
# 
# # marginal effects
# full_totrich_primgroup_marginal <- marginal_effects(full_totrich_primgroup_mod, method = "fitted")
# saveRDS(full_totrich_primgroup_marginal, "./FINAL/full/full_totrich_primgroup_marginal.RDS")
# 
# # specifically for plotting groupsize-by-threat effects
# full_totrich_primgroup_marginal_groupsize <- marginal_effects(full_totrich_primgroup_mod, effects = "logGroupSizePriUng_c:combIUCN", method = "fitted")
# saveRDS(full_totrich_primgroup_marginal_groupsize, "./FINAL/full/full_totrich_primgroup_marginal_groupsize.RDS")

# # ...ungulates
# full_totrich_unggroup_mod <- brm(
#   data = fullDat_unggroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logGroupSizePriUng_c +
#                  logHostSpeciesRange_c +
#                  logHostMass_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(full_totrich_unggroup_mod)
# plot(full_totrich_unggroup_mod)
# pp_check(full_totrich_unggroup_mod, nsamples = 500)
# 
# # add information criteria
# full_totrich_unggroup_mod <- readRDS("./FINAL/full/full_totrich_unggroup_mod.RDS")
# full_totrich_unggroup_mod <- add_ic(full_totrich_unggroup_mod, ic = "loo", reloo = TRUE)
# full_totrich_unggroup_mod <- add_ic(full_totrich_unggroup_mod, ic = "kfold")
# saveRDS(full_totrich_unggroup_mod, "./FINAL/full/full_totrich_unggroup_mod.RDS")

# # model fits, predictions and residuals
# full_totrich_unggroup_mod <- readRDS("./FINAL/full/full_totrich_unggroup_mod.RDS")
# full_totrich_unggroup_fitted <- fitted(full_totrich_unggroup_mod)
# full_totrich_unggroup_predict <- predict(full_totrich_unggroup_mod)
# full_totrich_unggroup_resid <- residuals(full_totrich_unggroup_mod, type = "pearson")
# 
# # save outputs
# saveRDS(full_totrich_unggroup_fitted, "./FINAL/full/full_totrich_unggroup_fitted.RDS")
# saveRDS(full_totrich_unggroup_predict, "./FINAL/full/full_totrich_unggroup_predict.RDS")
# saveRDS(full_totrich_unggroup_resid, "./FINAL/full/full_totrich_unggroup_resid.RDS")
# 
# # marginal effects
# full_totrich_unggroup_marginal <- marginal_effects(full_totrich_unggroup_mod, method = "fitted")
# saveRDS(full_totrich_unggroup_marginal, "./FINAL/full/full_totrich_unggroup_marginal.RDS")
# 
# # specifically for plotting groupsize-by-threat effects
# full_totrich_unggroup_marginal_groupsize <- marginal_effects(full_totrich_unggroup_mod, effects = "logGroupSizePriUng_c:combIUCN", method = "fitted")
# saveRDS(full_totrich_unggroup_marginal_groupsize, "./FINAL/full/full_totrich_unggroup_marginal_groupsize.RDS")


### PARASITE RICHNESS IS RESPONSE, SIMPLE MODELS ----
# # ...carnivores
# simple_totrich_carngroup_mod <- brm(
#   data = fullDat_carngroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simple_totrich_carngroup_mod)
# plot(simple_totrich_carngroup_mod)
# pp_check(simple_totrich_carngroup_mod, nsamples = 500)
# 
# # add information criteria
# simple_totrich_carngroup_mod <- add_ic(simple_totrich_carngroup_mod, ic = "loo", reloo = TRUE)
# simple_totrich_carngroup_mod <- add_ic(simple_totrich_carngroup_mod, ic = "kfold")
# saveRDS(simple_totrich_carngroup_mod, "./FINAL/simple/simple_totrich_carngroup_mod.RDS")

# # model fits, predictions and residuals
# simple_totrich_carngroup_mod <- readRDS("./FINAL/simple/simple_totrich_carngroup_mod.RDS")
# simple_totrich_carngroup_fitted <- fitted(simple_totrich_carngroup_mod)
# simple_totrich_carngroup_predict <- predict(simple_totrich_carngroup_mod)
# simple_totrich_carngroup_resid <- residuals(simple_totrich_carngroup_mod, type = "pearson")
# 
# # save outputs
# saveRDS(simple_totrich_carngroup_fitted, "./FINAL/simple/simple_totrich_carngroup_fitted.RDS")
# saveRDS(simple_totrich_carngroup_predict, "./FINAL/simple/simple_totrich_carngroup_predict.RDS")
# saveRDS(simple_totrich_carngroup_resid, "./FINAL/simple/simple_totrich_carngroup_resid.RDS")
# 
# # marginal effects
# simple_totrich_carngroup_marginal <- marginal_effects(simple_totrich_carngroup_mod, method = "fitted")
# saveRDS(simple_totrich_carngroup_marginal, "./FINAL/simple/simple_totrich_carngroup_marginal.RDS")
# 
# #...primates
# simple_totrich_primgroup_mod <- brm(
#   data = fullDat_primgroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simple_totrich_primgroup_mod)
# plot(simple_totrich_primgroup_mod)
# pp_check(simple_totrich_primgroup_mod, nsamples = 500)
# 
# # add information criteria
# simple_totrich_primgroup_mod <- add_ic(simple_totrich_primgroup_mod, ic = "loo", reloo = TRUE)
# simple_totrich_primgroup_mod <- add_ic(simple_totrich_primgroup_mod, ic = "kfold")
# saveRDS(simple_totrich_primgroup_mod, "./FINAL/simple/simple_totrich_primgroup_mod.RDS")
# 
# # model fits, predictions and residuals
# simple_totrich_primgroup_mod <- readRDS("./FINAL/simple/simple_totrich_primgroup_mod.RDS")
# simple_totrich_primgroup_fitted <- fitted(simple_totrich_primgroup_mod)
# simple_totrich_primgroup_predict <- predict(simple_totrich_primgroup_mod)
# simple_totrich_primgroup_resid <- residuals(simple_totrich_primgroup_mod, type = "pearson")
# 
# # save outputs
# saveRDS(simple_totrich_primgroup_fitted, "./FINAL/simple/simple_totrich_primgroup_fitted.RDS")
# saveRDS(simple_totrich_primgroup_predict, "./FINAL/simple/simple_totrich_primgroup_predict.RDS")
# saveRDS(simple_totrich_primgroup_resid, "./FINAL/simple/simple_totrich_primgroup_resid.RDS")
# 
# # marginal effects
# simple_totrich_primgroup_marginal <- marginal_effects(simple_totrich_primgroup_mod, method = "fitted")
# saveRDS(simple_totrich_primgroup_marginal, "./FINAL/simple/simple_totrich_primgroup_marginal.RDS")
# 
# #...ungulates
# simple_totrich_unggroup_mod <- brm(
#   data = fullDat_unggroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simple_totrich_unggroup_mod)
# plot(simple_totrich_unggroup_mod)
# pp_check(simple_totrich_unggroup_mod, nsamples = 500)
# 
# # add information criteria
# simple_totrich_unggroup_mod <- add_ic(simple_totrich_unggroup_mod, ic = "loo", reloo = TRUE)
# simple_totrich_unggroup_mod <- add_ic(simple_totrich_unggroup_mod, ic = "kfold")
# saveRDS(simple_totrich_unggroup_mod, "./FINAL/simple/simple_totrich_unggroup_mod.RDS")

# # model fits, predictions and residuals
# simple_totrich_unggroup_mod <- readRDS("./FINAL/simple/simple_totrich_unggroup_mod.RDS")
# simple_totrich_unggroup_fitted <- fitted(simple_totrich_unggroup_mod)
# simple_totrich_unggroup_predict <- predict(simple_totrich_unggroup_mod)
# simple_totrich_unggroup_resid <- residuals(simple_totrich_unggroup_mod, type = "pearson")
# 
# # save outputs
# saveRDS(simple_totrich_unggroup_fitted, "./FINAL/simple/simple_totrich_unggroup_fitted.RDS")
# saveRDS(simple_totrich_unggroup_predict, "./FINAL/simple/simple_totrich_unggroup_predict.RDS")
# saveRDS(simple_totrich_unggroup_resid, "./FINAL/simple/simple_totrich_unggroup_resid.RDS")
# 
# # marginal effects
# simple_totrich_unggroup_marginal <- marginal_effects(simple_totrich_unggroup_mod, method = "fitted")
# saveRDS(simple_totrich_unggroup_marginal, "./FINAL/simple/simple_totrich_unggroup_marginal.RDS")