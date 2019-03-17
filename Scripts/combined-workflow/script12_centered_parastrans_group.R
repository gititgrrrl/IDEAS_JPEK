### PARASITE TRANSMISSION FOR MODELS INCLUDING GROUP SIZE
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
# # PARASITE TRANSMISSION IS RESPONSE, COMPLEX MODELS ----
# # ...carnivores
# allDat <- read_csv("Data/JPEK/script4.csv")
# 
# fullDat_parastrans_carngroup <- allDat %>%
#   filter(hostGroup == "carnivores") %>%
#   select(hostName, parRichCloseOnly, parRichTransKnown, groupSizeCar, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() %>%
#   filter(parRichTransKnown > 0) # 139
# 
# fullDat_parastrans_carngroup <- fullDat_parastrans_carngroup[complete.cases(fullDat_parastrans_carngroup),] # only 71 records
# 
# fullDat_parastrans_carngroup$groupSizeCar <- factor(fullDat_parastrans_carngroup$groupSizeCar, levels = c("non_group", "group"))
# fullDat_parastrans_carngroup$combIUCN <- factor(fullDat_parastrans_carngroup$combIUCN, levels = c("not_threatened", "threatened"))
# 
# fullDat_parastrans_carngroup_c <- fullDat_parastrans_carngroup %>%
#   mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# # ...full model (some threat interactions omitted)
# fullBrm_carngroup_c_parastrans <- brm(
#   data = fullDat_parastrans_carngroup_c,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations_c +
#                  combIUCN*logHostSpeciesRange_c +
#                  combIUCN*groupSizeCar +
#                  logHostMass_c +
#                  hostMaxLifespan_c +
#                  absHostMeanLat_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# summary(fullBrm_carngroup_c_parastrans)
# saveRDS(fullBrm_carngroup_c_parastrans, "./Data/JPEK/full/full_brm_carngroup_c_parastrans.RDS") # model
# 
# # marginal effects
# fullMarg_carngroup_c_parastrans <- plot(marginal_effects(fullBrm_carngroup_c_parastrans), method = "fitted", plot = FALSE)
# saveRDS(fullMarg_carngroup_c_parastrans, "./Data/JPEK/full/full_brm_carngroup_c_parastrans_me.RDS")
#
# #...ungulates
# fullDat_parastrans_unggroup <- allDat %>%
#   filter(hostGroup == "ungulates") %>%
#   select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() %>%
#   filter(parRichTransKnown > 0) # 93
# 
# fullDat_parastrans_unggroup <- fullDat_parastrans_unggroup[complete.cases(fullDat_parastrans_unggroup),] # only 58 records
# 
# fullDat_parastrans_unggroup$combIUCN <- factor(fullDat_parastrans_unggroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_parastrans_unggroup_c <- fullDat_parastrans_unggroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# 
# # ...full model (some threat interactions omitted)
# fullBrm_unggroup_c_parastrans <- brm(
#   data = fullDat_parastrans_unggroup_c,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
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
# summary(fullBrm_unggroup_c_parastrans)
# saveRDS(fullBrm_unggroup_c_parastrans, "./Data/JPEK/full/full_brm_unggroup_c_parastrans.RDS") # model
# 
# # marginal effects
# fullMarg_unggroup_c_parastrans <- plot(marginal_effects(fullBrm_unggroup_c_parastrans), method = "fitted", plot = FALSE)
# saveRDS(fullMarg_unggroup_c_parastrans, "./Data/JPEK/full/full_brm_unggroup_c_parastrans_me.RDS")
#
# #...primates
# fullDat_parastrans_primgroup <- allDat %>%
#   filter(hostGroup == "primates") %>%
#   select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
#   distinct() %>%
#   filter(parRichTransKnown > 0) # 130
# 
# fullDat_parastrans_primgroup <- fullDat_parastrans_primgroup[complete.cases(fullDat_parastrans_primgroup),] # only 69 records
# 
# fullDat_parastrans_primgroup$combIUCN <- factor(fullDat_parastrans_primgroup$combIUCN, levels = c("not_threatened", "threatened"))
# fullDat_parastrans_primgroup_c <- fullDat_parastrans_primgroup %>%
#   mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
#          logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
#          logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
#          logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
#          hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
#          absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# # ...model with ALL threat interactions
# allIntBrm_primgroup_c_parastrans <- brm(
#   data = fullDat_parastrans_primgroup_c,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
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
# summary(allIntBrm_primgroup_c_parastrans)
# saveRDS(allIntBrm_primgroup_c_parastrans, "./Data/JPEK/allInt/allInt_brm_primgroup_c_parastrans.RDS")
# 
# # marginal effects
# allIntMarg_primgroup_c_parastrans <- plot(marginal_effects(allIntBrm_primgroup_c_parastrans), method = "fitted", plot = FALSE)
# saveRDS(allIntMarg_primgroup_c_parastrans, "./Data/JPEK/allInt/allInt_brm_primgroup_c_parastrans_me.RDS")
# 
# # PARASITE RICHNESS IS RESPONSE, SIMPLE MODELS ----
# ...carnivores
# simpBrm_carngroup_c_parastrans <- brm(
#   data = fullDat_parastrans_carngroup_c,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_carngroup_c_parastrans)
# saveRDS(simpBrm_carngroup_c_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_c_parastrans.RDS")
# 
# # marginal effects
# simpMarg_carngroup_c_parastrans <- plot(marginal_effects(simpBrm_carngroup_c_parastrans), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_carngroup_c_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_c_parastrans_me.RDS")
#
# #...ungulates
# simpBrm_unggroup_c_parastrans <- brm(
#   data = fullDat_parastrans_unggroup_c,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_unggroup_c_parastrans)
# saveRDS(simpBrm_unggroup_c_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_c_parastrans.RDS")
# 
# # marginal effects
# simpMarg_unggroup_c_parastrans <- plot(marginal_effects(simpBrm_unggroup_c_parastrans), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_unggroup_c_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_c_parastrans_me.RDS")
#
# #...primates
# simpBrm_primgroup_c_parastrans <- brm(
#   data = fullDat_parastrans_primgroup_c,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations_c),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_primgroup_c_parastrans)
# saveRDS(simpBrm_primgroup_c_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_c_parastrans.RDS")
# 
# # marginal effects
# simpMarg_primgroup_c_parastrans <- plot(marginal_effects(simpBrm_primgroup_c_parastrans), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_primgroup_c_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_c_parastrans_me.RDS")
