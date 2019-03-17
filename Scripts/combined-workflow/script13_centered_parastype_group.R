### PARASITE TYPE FOR MODELS INCLUDING GROUP SIZE
# # CENTER ALL THE DATA FOR MODELING ----
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
# # PARASITE TYPE IS RESPONSE, COMPLEX MODELS ----
# # ...carnivores
# allDat <- read_csv("Data/JPEK/script4.csv")
# 
# fullDat_parastype_carngroup <- allDat %>%
#   filter(hostGroup == "carnivores") %>%
#   select(hostName, logNumHostCitations, combIUCN, groupSizeCar, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() %>%
#   mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
#          parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 148
# 
# fullDat_parastype_carngroup <- fullDat_parastype_carngroup[complete.cases(fullDat_parastype_carngroup),] # only 74 records
# 
# fullDat_parastype_carngroup$groupSizeCar <- factor(fullDat_parastype_carngroup$groupSizeCar, levels = c("non_group", "group"))
# fullDat_parastype_carngroup$combIUCN <- factor(fullDat_parastype_carngroup$combIUCN, levels = c("not_threatened", "threatened"))
# 
# fullDat_parastype_carngroup_c <- fullDat_parastype_carngroup %>%
#   mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# # ...model with ALL threat interactions
# allIntBrm_carngroup_c_parastype <- brm(
#   data = fullDat_parastype_carngroup_c,
#   family = binomial,
#   formula = bf(parRich_micro|trials(parRich_alltypes) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*groupSizeCar +
#                  combIUCN*logHostMass_c +
#                  combIUCN*hostMaxLifespan_c +
#                  combIUCN*absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # # quick checks
# summary(allIntBrm_carngroup_c_parastype)
# saveRDS(allIntBrm_carngroup_c_parastype, "./Data/JPEK/allInt/allInt_brm_carngroup_c_parastype.RDS")
# 
# # marginal effects
# allIntMarg_carngroup_c_parastype <- plot(marginal_effects(allIntBrm_carngroup_c_parastype), method = "fitted", plot = FALSE)
# saveRDS(allIntMarg_carngroup_c_parastype, "./Data/JPEK/allInt/allInt_brm_carngroup_c_parastype_me.RDS")
# 
# #...ungulates
# fullDat_parastype_unggroup <- allDat %>%
#   filter(hostGroup == "ungulates") %>%
#   select(hostName, logNumHostCitations, combIUCN, groupSizePriUng, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() %>%
#   mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
#          parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 102 records
# 
# fullDat_parastype_unggroup <- fullDat_parastype_unggroup[complete.cases(fullDat_parastype_unggroup),] # only 60 records
# 
# fullDat_parastype_unggroup$combIUCN <- factor(fullDat_parastype_unggroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_parastype_unggroup_c <- fullDat_parastype_unggroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# # ...full model (some threat interactions omitted)
# fullBrm_unggroup_c_parastype <- brm(
#   data = fullDat_parastype_unggroup_c,
#   family = binomial,
#   formula = bf(parRich_micro|trials(parRich_alltypes) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*logGroupSizePriUng_c +
#                  logHostMass_c +
#                  hostMaxLifespan_c +
#                  absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(fullBrm_unggroup_c_parastype)
# saveRDS(fullBrm_unggroup_c_parastype, "./Data/JPEK/full/full_brm_unggroup_c_parastype.RDS") # model
# 
# # marginal effects
# fullMarg_unggroup_c_parastype <- plot(marginal_effects(fullBrm_unggroup_c_parastype), method = "fitted", plot = FALSE)
# saveRDS(fullMarg_unggroup_c_parastype, "./Data/JPEK/full/full_brm_unggroup_c_parastype_me.RDS")
# 
# #...primates
# fullDat_parastype_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, logNumHostCitations, combIUCN, groupSizePriUng, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() %>%
#   mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
#          parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 142 records
# 
# fullDat_parastype_primgroup <- fullDat_parastype_primgroup[complete.cases(fullDat_parastype_primgroup),] # only 73 records
# 
# fullDat_parastype_primgroup$combIUCN <- factor(fullDat_parastype_primgroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_parastype_primgroup_c <- fullDat_parastype_primgroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# # ...model with ALL threat interactions
# allIntBrm_primgroup_c_parastype <- brm(
#   data = fullDat_parastype_primgroup_c,
#   family = binomial,
#   formula = bf(parRich_micro|trials(parRich_alltypes) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*logGroupSizePriUng_c +
#                  combIUCN*logHostMass_c +
#                  combIUCN*hostMaxLifespan_c +
#                  combIUCN*absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(allIntBrm_primgroup_c_parastype)
# saveRDS(allIntBrm_primgroup_c_parastype, "./Data/JPEK/allInt/allInt_brm_primgroup_c_parastype.RDS")
# 
# # marginal effects
# allIntMarg_primgroup_c_parastype <- plot(marginal_effects(allIntBrm_primgroup_c_parastype), method = "fitted", plot = FALSE)
# saveRDS(allIntMarg_primgroup_c_parastype, "./Data/JPEK/allInt/allInt_brm_primgroup_c_parastype_me.RDS")
# 
# # PARASITE TYPE IS RESPONSE, SIMPLE MODELS ----
# # ...carnivores
# simpBrm_carngroup_c_parastype <- brm(
#   data = fullDat_parastype_carngroup_c,
#   family = binomial,
#   formula = bf(parRich_micro|trials(parRich_alltypes) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_carngroup_c_parastype)
# saveRDS(simpBrm_carngroup_c_parastype, "./Data/JPEK/simple/simp_brm_carngroup_c_parastype.RDS")
# 
# # marginal effects
# simpMarg_carngroup_c_parastype <- plot(marginal_effects(simpBrm_carngroup_c_parastype), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_carngroup_c_parastype, "./Data/JPEK/simple/simp_brm_carngroup_c_parastype_me.RDS")
# 
# #...ungulates
# simpBrm_unggroup_c_parastype <- brm(
#   data = fullDat_parastype_unggroup_c,
#   family = binomial,
#   formula = bf(parRich_micro|trials(parRich_alltypes) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_unggroup_c_parastype)
# saveRDS(simpBrm_unggroup_c_parastype, "./Data/JPEK/simple/simp_brm_unggroup_c_parastype.RDS")
# 
# # marginal effects
# simpMarg_unggroup_c_parastype <- plot(marginal_effects(simpBrm_unggroup_c_parastype), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_unggroup_c_parastype, "./Data/JPEK/simple/simp_brm_unggroup_c_parastype_me.RDS")
# 
# #...primates
# simpBrm_primgroup_c_parastype <- brm(
#   data = fullDat_parastype_primgroup_c,
#   family = binomial,
#   formula = bf(parRich_micro|trials(parRich_alltypes) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_primgroup_c_parastype)
# saveRDS(simpBrm_primgroup_c_parastype, "./Data/JPEK/simple/simp_brm_primgroup_c_parastype.RDS")
# 
# # marginal effects
# simpMarg_primgroup_c_parastype <- plot(marginal_effects(simpBrm_primgroup_c_parastype), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_primgroup_c_parastype, "./Data/JPEK/simple/simp_brm_primgroup_c_parastype_me.RDS")
