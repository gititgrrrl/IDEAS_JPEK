### PARASITE TYPE FOR MODELS INCLUDING GROUP SIZE

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

### CARNIVORES WITH PARASITE TYPE AS RESPONSE ----
# Create data
allDat <- read_csv("Data/JPEK/script4.csv")

fullDat_parastype_carngroup <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, logNumHostCitations, combIUCN, groupSizeCar, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 393

fullDat_parastype_carngroup <- fullDat_parastype_carngroup[complete.cases(fullDat_parastype_carngroup),] # only 74 records

# ...model with ALL threat interactions
allIntBrm_carngroup_parastype <- brm(
  data = fullDat_parastype_carngroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*groupSizeCar +
                 combIUCN*logHostMass +
                 combIUCN*hostMaxLifespan +
                 combIUCN*absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(allIntBrm_carngroup_parastype)
# plot(allIntBrm_carngroup_parastype)
# pp_check(allIntBrm_carngroup_parastype, nsamples = 500)

# add information criteria
allIntBrm_carngroup_parastype <- add_ic(allIntBrm_carngroup_parastype, ic = "loo", reloo = TRUE)
allIntBrm_carngroup_parastype <- add_ic(allIntBrm_carngroup_parastype, ic = "kfold")
saveRDS(allIntBrm_carngroup_parastype, "./Data/JPEK/allInt/allInt_brm_carngroup_parastype.RDS")

# marginal effects
allIntMarg_carngroup_parastype <- plot(marginal_effects(allIntBrm_carngroup_parastype), method = "fitted", plot = FALSE)
saveRDS(allIntMarg_carngroup_parastype, "./Data/JPEK/allInt/allInt_brm_carngroup_parastype_me.RDS")

# ...full model (some threat interactions omitted)
fullBrm_carngroup_parastype <- brm(
  data = fullDat_parastype_carngroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*groupSizeCar +
                 logHostMass +
                 hostMaxLifespan +
                 absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(fullBrm_carngroup_parastype)
# plot(fullBrm_carngroup_parastype)
# pp_check(fullBrm_carngroup_parastype, nsamples = 500)

# add information criteria
fullBrm_carngroup_parastype <- add_ic(fullBrm_carngroup_parastype, ic = "loo", reloo = TRUE)
fullBrm_carngroup_parastype <- add_ic(fullBrm_carngroup_parastype, ic = "kfold")
saveRDS(fullBrm_carngroup_parastype, "./Data/JPEK/full/full_brm_carngroup_parastype.RDS")

# # model fits, predictions and residuals
# fullMu_carngroup_parastype <- fitted(fullBrm_carngroup_parastype)
# fullPredict_carngroup_parastype <- predict(fullBrm_carngroup_parastype)
# fullResid_carngroup_parastype <- residuals(fullBrm_carngroup_parastype, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_carngroup_parastype, "./Data/JPEK/full/full_brm_carngroup_parastype_mu.RDS")
# saveRDS(fullPredict_carngroup_parastype, "./Data/JPEK/full/full_brm_carngroup_parastype_predict.RDS")
# saveRDS(fullResid_carngroup_parastype, "./Data/JPEK/full/full_brm_carngroup_parastype_resid.RDS")
# 
# # marginal effects
# fullMarg_carngroup_parastype <- plot(marginal_effects(fullBrm_carngroup_parastype), method = "fitted", plot = FALSE)
# saveRDS(fullMarg_carngroup_parastype, "./Data/JPEK/full/full_brm_carngroup_parastype_me.RDS")

### UNGULATES WITH PARASITE TYPE AS RESPONSE ----
# Create data
fullDat_parastype_unggroup <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, logNumHostCitations, combIUCN, groupSizePriUng, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich))

fullDat_parastype_unggroup <- fullDat_parastype_unggroup[complete.cases(fullDat_parastype_unggroup),] # only 60 records

# ...model with ALL threat interactions
allIntBrm_unggroup_parastype <- brm(
  data = fullDat_parastype_unggroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 combIUCN*logHostMass +
                 combIUCN*hostMaxLifespan +
                 combIUCN*absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(allIntBrm_unggroup_parastype)
# plot(allIntBrm_unggroup_parastype)
# pp_check(allIntBrm_unggroup_parastype, nsamples = 500)

# add information criteria
allIntBrm_unggroup_parastype <- add_ic(allIntBrm_unggroup_parastype, ic = "loo", reloo = TRUE)
allIntBrm_unggroup_parastype <- add_ic(allIntBrm_unggroup_parastype, ic = "kfold")
saveRDS(allIntBrm_unggroup_parastype, "./Data/JPEK/allInt/allInt_brm_unggroup_parastype.RDS")

# ...full model (some threat interactions omitted)
fullBrm_unggroup_parastype <- brm(
  data = fullDat_parastype_unggroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 logHostMass +
                 hostMaxLifespan +
                 absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(fullBrm_unggroup_parastype)
# plot(fullBrm_unggroup_parastype)
# pp_check(fullBrm_unggroup_parastype, nsamples = 500)

# add information criteria
fullBrm_unggroup_parastype <- add_ic(fullBrm_unggroup_parastype, ic = "loo", reloo = TRUE)
fullBrm_unggroup_parastype <- add_ic(fullBrm_unggroup_parastype, ic = "kfold")
saveRDS(fullBrm_unggroup_parastype, "./Data/JPEK/full/full_brm_unggroup_parastype.RDS")

# # model fits, predictions and residuals
# fullMu_unggroup_parastype <- fitted(fullBrm_unggroup_parastype)
# fullPredict_unggroup_parastype <- predict(fullBrm_unggroup_parastype)
# fullResid_unggroup_parastype <- residuals(fullBrm_unggroup_parastype, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_unggroup_parastype, "./Data/JPEK/full/full_brm_unggroup_parastype_mu.RDS")
# saveRDS(fullPredict_unggroup_parastype, "./Data/JPEK/full/full_brm_unggroup_parastype_predict.RDS")
# saveRDS(fullResid_unggroup_parastype, "./Data/JPEK/full/full_brm_unggroup_parastype_resid.RDS")
# 
# marginal effects
fullMarg_unggroup_parastype <- plot(marginal_effects(fullBrm_unggroup_parastype), method = "fitted", plot = FALSE)
saveRDS(fullMarg_unggroup_parastype, "./Data/JPEK/full/full_brm_unggroup_parastype_me.RDS")

### PRIMATES WITH PARASITE TYPE AS RESPONSE ----
# Create data
fullDat_parastype_primgroup <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, logNumHostCitations, combIUCN, groupSizePriUng, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich))

fullDat_parastype_primgroup <- fullDat_parastype_primgroup[complete.cases(fullDat_parastype_primgroup),] # only 60 records

# ...model with ALL threat interactions
allIntBrm_primgroup_parastype <- brm(
  data = fullDat_parastype_primgroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 combIUCN*logHostMass +
                 combIUCN*hostMaxLifespan +
                 combIUCN*absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(allIntBrm_primgroup_parastype)
# plot(allIntBrm_primgroup_parastype)
# pp_check(allIntBrm_primgroup_parastype, nsamples = 500)

# add information criteria
allIntBrm_primgroup_parastype <- add_ic(allIntBrm_primgroup_parastype, ic = "loo", reloo = TRUE)
allIntBrm_primgroup_parastype <- add_ic(allIntBrm_primgroup_parastype, ic = "kfold")
saveRDS(allIntBrm_primgroup_parastype, "./Data/JPEK/allInt/allInt_brm_primgroup_parastype.RDS")

# marginal effects
allIntMarg_primgroup_parastype <- plot(marginal_effects(allIntBrm_primgroup_parastype), method = "fitted", plot = FALSE)
saveRDS(allIntMarg_primgroup_parastype, "./Data/JPEK/allInt/allInt_brm_primgroup_parastype_me.RDS")

# ...full model (some threat interactions omitted)
fullBrm_primgroup_parastype <- brm(
  data = fullDat_parastype_primgroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 logHostMass +
                 hostMaxLifespan +
                 absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(fullBrm_primgroup_parastype)
plot(fullBrm_primgroup_parastype)
pp_check(fullBrm_primgroup_parastype, nsamples = 500)
saveRDS(fullBrm_primgroup_parastype, "./Data/JPEK/full/full_brm_primgroup_parastype.RDS") # model

# add information criteria
fullBrm_primgroup_parastype <- add_ic(fullBrm_primgroup_parastype, ic = "loo", reloo = TRUE)
fullBrm_primgroup_parastype <- add_ic(fullBrm_primgroup_parastype, ic = "kfold")
saveRDS(fullBrm_primgroup_parastype, "./Data/JPEK/full/full_brm_primgroup_parastype.RDS")

# # model fits, predictions and residuals
# fullMu_primgroup_parastype <- fitted(fullBrm_primgroup_parastype)
# fullPredict_primgroup_parastype <- predict(fullBrm_primgroup_parastype)
# fullResid_primgroup_parastype <- residuals(fullBrm_primgroup_parastype, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_primgroup_parastype, "./Data/JPEK/full/full_brm_primgroup_parastype_mu.RDS")
# saveRDS(fullPredict_primgroup_parastype, "./Data/JPEK/full/full_brm_primgroup_parastype_predict.RDS")
# saveRDS(fullResid_primgroup_parastype, "./Data/JPEK/full/full_brm_primgroup_parastype_resid.RDS")
# 
# # marginal effects
# fullMarg_primgroup_parastype <- plot(marginal_effects(fullBrm_primgroup_parastype), method = "fitted", plot = FALSE)
# saveRDS(fullMarg_primgroup_parastype, "./Data/JPEK/full/full_brm_primgroup_parastype_me.RDS")



## SIMPLE MODEL, CARNIVORES WITH PARASITE TYPE AS RESPONSE ----
# ...simple model
simpBrm_carngroup_parastype <- brm(
  data = fullDat_parastype_carngroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(simpBrm_carngroup_parastype)
# plot(simpBrm_carngroup_parastype)
# pp_check(simpBrm_carngroup_parastype, nsamples = 500)

# add information criteria
simpBrm_carngroup_parastype <- add_ic(simpBrm_carngroup_parastype, ic = "loo", reloo = TRUE)
simpBrm_carngroup_parastype <- add_ic(simpBrm_carngroup_parastype, ic = "kfold")
saveRDS(simpBrm_carngroup_parastype, "./Data/JPEK/simple/simp_brm_carngroup_parastype.RDS")

# # model fits, predictions and residuals
# simpMu_carngroup_parastype <- fitted(simpBrm_carngroup_parastype)
# simpPredict_carngroup_parastype <- predict(simpBrm_carngroup_parastype)
# simpResid_carngroup_parastype <- residuals(simpBrm_carngroup_parastype, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_carngroup_parastype, "./Data/JPEK/simple/simp_brm_carngroup_parastype_mu.RDS")
# saveRDS(simpPredict_carngroup_parastype, "./Data/JPEK/simple/simp_brm_carngroup_parastype_predict.RDS")
# saveRDS(simpResid_carngroup_parastype, "./Data/JPEK/simple/simp_brm_carngroup_parastype_resid.RDS")

# marginal effects
simpMarg_carngroup_parastype <- plot(marginal_effects(simpBrm_carngroup_parastype), method = "fitted", plot = FALSE)
saveRDS(simpMarg_carngroup_parastype, "./Data/JPEK/simple/simp_brm_carngroup_parastype_me.RDS")

### UNGULATES WITH PARASITE TYPE AS RESPONSE ----
# ...simple model
simpBrm_unggroup_parastype <- brm(
  data = fullDat_parastype_unggroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(simpBrm_unggroup_parastype)
# plot(simpBrm_unggroup_parastype)
# pp_check(simpBrm_unggroup_parastype, nsamples = 500)

# add information criteria
simpBrm_unggroup_parastype <- add_ic(simpBrm_unggroup_parastype, ic = "loo", reloo = TRUE)
simpBrm_unggroup_parastype <- add_ic(simpBrm_unggroup_parastype, ic = "kfold")
saveRDS(simpBrm_unggroup_parastype, "./Data/JPEK/simple/simp_brm_unggroup_parastype.RDS")

# # model fits, predictions and residuals
# simpMu_unggroup_parastype <- fitted(simpBrm_unggroup_parastype)
# simpPredict_unggroup_parastype <- predict(simpBrm_unggroup_parastype)
# simpResid_unggroup_parastype <- residuals(simpBrm_unggroup_parastype, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_unggroup_parastype, "./Data/JPEK/simple/simp_brm_unggroup_parastype_mu.RDS")
# saveRDS(simpPredict_unggroup_parastype, "./Data/JPEK/simple/simp_brm_unggroup_parastype_predict.RDS")
# saveRDS(simpResid_unggroup_parastype, "./Data/JPEK/simple/simp_brm_unggroup_parastype_resid.RDS")

# marginal effects
simpMarg_unggroup_parastype <- plot(marginal_effects(simpBrm_unggroup_parastype), method = "fitted", plot = FALSE)
saveRDS(simpMarg_unggroup_parastype, "./Data/JPEK/simple/simp_brm_unggroup_parastype_me.RDS")

### PRIMATES WITH PARASITE TYPE AS RESPONSE ----
# ...simple model
simpBrm_primgroup_parastype <- brm(
  data = fullDat_parastype_primgroup,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
# summary(simpBrm_primgroup_parastype)
# plot(simpBrm_primgroup_parastype)
# pp_check(simpBrm_primgroup_parastype, nsamples = 500)

# add information criteria
simpBrm_primgroup_parastype <- add_ic(simpBrm_primgroup_parastype, ic = "loo", reloo = TRUE)
simpBrm_primgroup_parastype <- add_ic(simpBrm_primgroup_parastype, ic = "kfold")
saveRDS(simpBrm_primgroup_parastype, "./Data/JPEK/simple/simp_brm_primgroup_parastype.RDS")

# # model fits, predictions and residuals
# simpMu_primgroup_parastype <- fitted(simpBrm_primgroup_parastype)
# simpPredict_primgroup_parastype <- predict(simpBrm_primgroup_parastype)
# simpResid_primgroup_parastype <- residuals(simpBrm_primgroup_parastype, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_primgroup_parastype, "./Data/JPEK/simple/simp_brm_primgroup_parastype_mu.RDS")
# saveRDS(simpPredict_primgroup_parastype, "./Data/JPEK/simple/simp_brm_primgroup_parastype_predict.RDS")
# saveRDS(simpResid_primgroup_parastype, "./Data/JPEK/simple/simp_brm_primgroup_parastype_resid.RDS")
# 
# marginal effects
simpMarg_primgroup_parastype <- plot(marginal_effects(simpBrm_primgroup_parastype), method = "fitted", plot = FALSE)
saveRDS(simpMarg_primgroup_parastype, "./Data/JPEK/simple/simp_brm_primgroup_parastype_me.RDS")