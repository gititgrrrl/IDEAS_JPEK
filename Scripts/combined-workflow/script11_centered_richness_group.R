### PARASITE RICHNESS FOR MODELS INCLUDING GROUP SIZE
# # CENTER ALL THE DATA FOR MODELING ----
# # > I realized that if we want readers to be able to relate marginal effects plots to the model outputs I needed to center all predictors (I should have done this anyway as good stats practice--totally forgot). I also re-leveled the categorical variables, so reference level is non-threatened (and for carnivores, non-group)
# 
# rm(list=ls())
# 
# ### LOAD PACKAGES ----
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
#
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
#
# # PARASITE RICHNESS IS RESPONSE, COMPLEX MODELS ----
# # ...carnivores
# fullDat_carngroup <- allDat %>%
#   filter(hostGroup == "carnivores") %>%
#   select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() # 148 records
# fullDat_carngroup <- fullDat_carngroup[complete.cases(fullDat_carngroup),] # only 74 records
# 
# fullDat_carngroup$groupSizeCar <- factor(fullDat_carngroup$groupSizeCar, levels = c("non_group", "group"))
# fullDat_carngroup$combIUCN <- factor(fullDat_carngroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_carngroup_c <- fullDat_carngroup %>%
#   mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# allIntBrm_carngroup_c <- brm(
#   data = fullDat_carngroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*groupSizeCar +
#                  combIUCN*logHostMass_c +
#                  combIUCN*hostMaxLifespan_c +
#                  combIUCN*absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# summary(allIntBrm_carngroup_c)
# saveRDS(allIntBrm_carngroup_c, "./Data/JPEK/allInt/allInt_brm_carngroup_c_all.RDS")

# marginal effects
allIntBrm_carngroup_c <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_all.RDS")
allInt_carngroup_c_me <- plot(marginal_effects(allIntBrm_carngroup_c), method = "fitted", plot = FALSE)
saveRDS(allInt_carngroup_c_me, "./Data/JPEK/allInt/allInt_brm_carngroup_c_me.RDS")
# specifically for plotting groupsize-by-threat effects
allInt_carngroup_c_me_groupsize <- marginal_effects(allIntBrm_carngroup_c, effects = "groupSizeCar:combIUCN", method = "fitted")
saveRDS(allInt_carngroup_c_me_groupsize, "./Data/JPEK/allInt/allInt_brm_carngroup_c_me_groupsize.RDS")

# #...primates
# fullDat_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 142 records
# fullDat_primgroup <- fullDat_primgroup[complete.cases(fullDat_primgroup),] # only 73 records
# 
# fullDat_primgroup$combIUCN <- factor(fullDat_primgroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_primgroup_c <- fullDat_primgroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE), 
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# fullBrm_primgroup_c <- brm(
#   data = fullDat_primgroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*logGroupSizePriUng_c +
#                  logHostMass_c +
#                  hostMaxLifespan_c +
#                  absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# saveRDS(fullBrm_primgroup_c, "./Data/JPEK/full/full_brm_primgroup_c_all.RDS") # model
# 
# marginal effects
fullBrm_primgroup_c <- readRDS("./Data/JPEK/full/full_brm_primgroup_c_all.RDS")
full_primgroup_c_me <- plot(marginal_effects(fullBrm_primgroup_c), method = "fitted", plot = FALSE)
saveRDS(full_primgroup_c_me, "./Data/JPEK/full/full_brm_primgroup_c_me.RDS")
# specifically for plotting groupsize-by-threat effects
full_primgroup_c_me_groupsize <- marginal_effects(fullBrm_primgroup_c, effects = "logGroupSizePriUng_c:combIUCN", method = "fitted")
saveRDS(full_primgroup_c_me_groupsize, "./Data/JPEK/full/full_brm_primgroup_c_me_groupsize.RDS")
# specifically for plotting range-by-threat effects
full_primgroup_c_me_speciesrange <- marginal_effects(fullBrm_primgroup_c, effects = "logHostSpeciesRange_c:combIUCN", method = "fitted")
saveRDS(full_primgroup_c_me_speciesrange, "./Data/JPEK/full/full_brm_primgroup_c_me_speciesrange.RDS")
#
# #...ungulates
# fullDat_unggroup <- allDat %>%
#   filter(hostGroup == "ungulates") %>%
#   select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() # 102 records
# fullDat_unggroup <- fullDat_unggroup[complete.cases(fullDat_unggroup),] # only 60 records
# 
# fullDat_unggroup$combIUCN <- factor(fullDat_unggroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_unggroup_c <- fullDat_unggroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# # ...full model (some threat interactions omitted)
# fullBrm_unggroup_c <- brm(
#   data = fullDat_unggroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*logGroupSizePriUng_c +
#                  logHostMass_c +
#                  hostMaxLifespan_c +
#                  absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
#
# saveRDS(fullBrm_unggroup_c, "./Data/JPEK/full/full_brm_unggroup_c_all.RDS")

# marginal effects
fullBrm_unggroup_c <- readRDS("./Data/JPEK/full/full_brm_unggroup_c_all.RDS")
full_unggroup_c_me <- plot(marginal_effects(fullBrm_unggroup_c), method = "fitted", plot = FALSE)
saveRDS(full_unggroup_c_me, "./Data/JPEK/full/full_brm_unggroup_c_me.RDS")
# specifically for plotting groupsize-by-threat effects
full_unggroup_c_me_groupsize <- marginal_effects(fullBrm_unggroup_c, effects = "logGroupSizePriUng_c:combIUCN", method = "fitted")
saveRDS(full_unggroup_c_me_groupsize, "./Data/JPEK/full/full_brm_unggroup_c_me_groupsize.RDS")
# specifically for plotting range-by-threat effects
full_unggroup_c_me_speciesrange <- marginal_effects(fullBrm_unggroup_c, effects = "logHostSpeciesRange_c:combIUCN", method = "fitted")
saveRDS(full_unggroup_c_me_speciesrange, "./Data/JPEK/full/full_brm_unggroup_c_me_speciesrange.RDS")
#
# # PARASITE RICHNESS IS RESPONSE, SIMPLE MODELS ----
# ...carnivores
# simpBrm_carngroup_c <- brm(
#   data = fullDat_carngroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# summary(simpBrm_carngroup_c)
# 
# saveRDS(simpBrm_carngroup_c, "./Data/JPEK/simple/simp_brm_carngroup_c_all.RDS") # model
# 
# marginal effects
simpBrm_carngroup_c <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_all.RDS") # model
simp_carngroup_c_me <- plot(marginal_effects(simpBrm_carngroup_c), method = "fitted", plot = FALSE)
saveRDS(simp_carngroup_c_me, "./Data/JPEK/simple/simp_brm_carngroup_c_me.RDS")
#
# #...primates
# simpBrm_primgroup_c <- brm(
#   data = fullDat_primgroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# summary(simpBrm_primgroup_c)
# saveRDS(simpBrm_primgroup_c, "./Data/JPEK/simple/simp_brm_primgroup_c_all.RDS") # model
#
# marginal effects
simpBrm_primgroup_c <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_all.RDS")
simp_primgroup_c_me <- plot(marginal_effects(simpBrm_primgroup_c), method = "fitted", plot = FALSE)
saveRDS(simp_primgroup_c_me, "./Data/JPEK/simple/simp_brm_primgroup_c_me.RDS")
#
# #...ungulates
# simpBrm_unggroup_c <- brm(
#   data = fullDat_unggroup_c,
#   family = poisson,
#   formula = bf(parRich | trunc(lb = 1) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# summary(simpBrm_unggroup_c)
# saveRDS(simpBrm_unggroup_c, "./Data/JPEK/simple/simp_brm_unggroup_c_all.RDS") # model
#  
# marginal effects
simpBrm_unggroup_c <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_all.RDS")
simp_unggroup_c_me <- plot(marginal_effects(simpBrm_unggroup_c), method = "fitted", plot = FALSE)
saveRDS(simp_unggroup_c_me, "./Data/JPEK/simple/simp_brm_unggroup_c_me.RDS")
