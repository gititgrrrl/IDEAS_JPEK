### PARASITE TRANSMISSION FOR MODELS INCLUDING GROUP SIZE

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

### CARNIVORES WITH PARASITE TRANSMISSION AS RESPONSE ----
# Create data
allDat <- read_csv("Data/JPEK/script4.csv")

fullDat_parastrans_carngroup <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizeCar, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

fullDat_parastrans_carngroup <- fullDat_parastrans_carngroup[complete.cases(fullDat_parastrans_carngroup),] # only 71 records

# ...model with ALL threat interactions
allIntBrm_carngroup_parastrans <- brm(
  data = fullDat_parastrans_carngroup,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*groupSizeCar +
                 combIUCN*logHostMass +
                 combIUCN*hostMaxLifespan +
                 combIUCN*absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(allIntBrm_carngroup_parastrans)
plot(allIntBrm_carngroup_parastrans)
pp_check(allIntBrm_carngroup_parastrans, nsamples = 500)

# add information criteria
allIntBrm_carngroup_parastrans <- add_ic(allIntBrm_carngroup_parastrans, ic = "loo", reloo = TRUE)
allIntBrm_carngroup_parastrans <- add_ic(allIntBrm_carngroup_parastrans, ic = "kfold")
saveRDS(allIntBrm_carngroup_parastrans, "./Data/JPEK/allInt/allInt_brm_carngroup_parastrans.RDS")

# ...full model (some threat interactions omitted)
fullBrm_carngroup_parastrans <- brm(
  data = fullDat_parastrans_carngroup,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*groupSizeCar +
                 logHostMass +
                 hostMaxLifespan +
                 absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(fullBrm_carngroup_parastrans)
plot(fullBrm_carngroup_parastrans)
pp_check(fullBrm_carngroup_parastrans, nsamples = 500)
saveRDS(fullBrm_carngroup_parastrans, "./Data/JPEK/full/full_brm_carngroup_parastrans.RDS") # model

# add information criteria
fullBrm_carngroup_parastrans <- add_ic(fullBrm_carngroup_parastrans, ic = "loo", reloo = TRUE)
fullBrm_carngroup_parastrans <- add_ic(fullBrm_carngroup_parastrans, ic = "kfold")
saveRDS(fullBrm_carngroup_parastrans, "./Data/JPEK/full/full_brm_carngroup_parastrans.RDS")

# # model fits, predictions and residuals
# fullMu_carngroup_parastrans <- fitted(fullBrm_carngroup_parastrans)
# fullPredict_carngroup_parastrans <- predict(fullBrm_carngroup_parastrans)
# fullResid_carngroup_parastrans <- residuals(fullBrm_carngroup_parastrans, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_carngroup_parastrans, "./Data/JPEK/full/full_brm_carngroup_parastrans_mu.RDS")
# saveRDS(fullPredict_carngroup_parastrans, "./Data/JPEK/full/full_brm_carngroup_parastrans_predict.RDS")
# saveRDS(fullResid_carngroup_parastrans, "./Data/JPEK/full/full_brm_carngroup_parastrans_resid.RDS")
# 
# marginal effects
fullMarg_carngroup_parastrans <- plot(marginal_effects(fullBrm_carngroup_parastrans), method = "fitted", plot = FALSE)
saveRDS(fullMarg_carngroup_parastrans, "./Data/JPEK/full/full_brm_carngroup_parastrans_me.RDS")

### UNGULATES WITH PARASITE TRANSMISSION AS RESPONSE ----
# Create data
fullDat_parastrans_unggroup <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

fullDat_parastrans_unggroup <- fullDat_parastrans_unggroup[complete.cases(fullDat_parastrans_unggroup),] # only 60 records

# ...model with ALL threat interactions
allIntBrm_unggroup_parastrans <- brm(
  data = fullDat_parastrans_unggroup,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 combIUCN*logHostMass +
                 combIUCN*hostMaxLifespan +
                 combIUCN*absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(allIntBrm_unggroup_parastrans)
plot(allIntBrm_unggroup_parastrans)
pp_check(allIntBrm_unggroup_parastrans, nsamples = 500)

# add information criteria
allIntBrm_unggroup_parastrans <- add_ic(allIntBrm_unggroup_parastrans, ic = "loo", reloo = TRUE)
allIntBrm_unggroup_parastrans <- add_ic(allIntBrm_unggroup_parastrans, ic = "kfold")
saveRDS(allIntBrm_unggroup_parastrans, "./Data/JPEK/allInt/allInt_brm_unggroup_parastrans.RDS")

# ...full model (some threat interactions omitted)
fullBrm_unggroup_parastrans <- brm(
  data = fullDat_parastrans_unggroup,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 logHostMass +
                 hostMaxLifespan +
                 absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(fullBrm_unggroup_parastrans)
plot(fullBrm_unggroup_parastrans)
pp_check(fullBrm_unggroup_parastrans, nsamples = 500)
saveRDS(fullBrm_unggroup_parastrans, "./Data/JPEK/full/full_brm_unggroup_parastrans.RDS") # model

# add information criteria
fullBrm_unggroup_parastrans <- add_ic(fullBrm_unggroup_parastrans, ic = "loo", reloo = TRUE)
fullBrm_unggroup_parastrans <- add_ic(fullBrm_unggroup_parastrans, ic = "kfold")
saveRDS(fullBrm_unggroup_parastrans, "./Data/JPEK/full/full_brm_unggroup_parastrans.RDS")

# # model fits, predictions and residuals
# fullMu_unggroup_parastrans <- fitted(fullBrm_unggroup_parastrans)
# fullPredict_unggroup_parastrans <- predict(fullBrm_unggroup_parastrans)
# fullResid_unggroup_parastrans <- residuals(fullBrm_unggroup_parastrans, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_unggroup_parastrans, "./Data/JPEK/full/full_brm_unggroup_parastrans_mu.RDS")
# saveRDS(fullPredict_unggroup_parastrans, "./Data/JPEK/full/full_brm_unggroup_parastrans_predict.RDS")
# saveRDS(fullResid_unggroup_parastrans, "./Data/JPEK/full/full_brm_unggroup_parastrans_resid.RDS")
# 
# marginal effects
fullMarg_unggroup_parastrans <- plot(marginal_effects(fullBrm_unggroup_parastrans), method = "fitted", plot = FALSE)
saveRDS(fullMarg_unggroup_parastrans, "./Data/JPEK/full/full_brm_unggroup_parastrans_me.RDS")

### PRIMATES WITH PARASITE TRANSMISSION AS RESPONSE ----
# Create data
fullDat_parastrans_primgroup <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

fullDat_parastrans_primgroup <- fullDat_parastrans_primgroup[complete.cases(fullDat_parastrans_primgroup),] # only 60 records

# ...model with ALL threat interactions
allIntBrm_primgroup_parastrans <- brm(
  data = fullDat_parastrans_primgroup,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 combIUCN*logHostMass +
                 combIUCN*hostMaxLifespan +
                 combIUCN*absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(allIntBrm_primgroup_parastrans)
plot(allIntBrm_primgroup_parastrans)
pp_check(allIntBrm_primgroup_parastrans, nsamples = 500)

# add information criteria
allIntBrm_primgroup_parastrans <- add_ic(allIntBrm_primgroup_parastrans, ic = "loo", reloo = TRUE)
allIntBrm_primgroup_parastrans <- add_ic(allIntBrm_primgroup_parastrans, ic = "kfold")
saveRDS(allIntBrm_primgroup_parastrans, "./Data/JPEK/allInt/allInt_brm_primgroup_parastrans.RDS")

# marginal effects
allIntMarg_primgroup_parastrans <- plot(marginal_effects(allIntBrm_primgroup_parastrans), method = "fitted", plot = FALSE)
saveRDS(allIntMarg_primgroup_parastrans, "./Data/JPEK/allInt/allInt_brm_primgroup_parastrans_me.RDS")

# ...full model (some threat interactions omitted)
fullBrm_primgroup_parastrans <- brm(
  data = fullDat_parastrans_primgroup,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations +
                 combIUCN*logHostSpeciesRange +
                 combIUCN*logGroupSizePriUng +
                 logHostMass +
                 hostMaxLifespan +
                 absHostMeanLat),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(fullBrm_primgroup_parastrans)
plot(fullBrm_primgroup_parastrans)
pp_check(fullBrm_primgroup_parastrans, nsamples = 500)
saveRDS(fullBrm_primgroup_parastrans, "./Data/JPEK/full/full_brm_primgroup_parastrans.RDS") # model

# add information criteria
fullBrm_primgroup_parastrans <- add_ic(fullBrm_primgroup_parastrans, ic = "loo", reloo = TRUE)
fullBrm_primgroup_parastrans <- add_ic(fullBrm_primgroup_parastrans, ic = "kfold")
saveRDS(fullBrm_primgroup_parastrans, "./Data/JPEK/full/full_brm_primgroup_parastrans.RDS")

# # model fits, predictions and residuals
# fullMu_primgroup_parastrans <- fitted(fullBrm_primgroup_parastrans)
# fullPredict_primgroup_parastrans <- predict(fullBrm_primgroup_parastrans)
# fullResid_primgroup_parastrans <- residuals(fullBrm_primgroup_parastrans, type = "pearson")
# 
# # save outputs
# saveRDS(fullMu_primgroup_parastrans, "./Data/JPEK/full/full_brm_primgroup_parastrans_mu.RDS")
# saveRDS(fullPredict_primgroup_parastrans, "./Data/JPEK/full/full_brm_primgroup_parastrans_predict.RDS")
# saveRDS(fullResid_primgroup_parastrans, "./Data/JPEK/full/full_brm_primgroup_parastrans_resid.RDS")
# 
# # marginal effects
# fullMarg_primgroup_parastrans <- plot(marginal_effects(fullBrm_primgroup_parastrans), method = "fitted", plot = FALSE)
# saveRDS(fullMarg_primgroup_parastrans, "./Data/JPEK/full/full_brm_primgroup_parastrans_me.RDS")
# 
# 
# 
# ## SIMPLE MODEL, CARNIVORES WITH PARASITE TRANSMISSION AS RESPONSE ----
# # ...simple model
# simpBrm_carngroup_parastrans <- brm(
#   data = fullDat_parastrans_carngroup,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_carngroup_parastrans)
# plot(simpBrm_carngroup_parastrans)
# pp_check(simpBrm_carngroup_parastrans, nsamples = 500)
# 
# # add information criteria
# simpBrm_carngroup_parastrans <- add_ic(simpBrm_carngroup_parastrans, ic = "loo", reloo = TRUE)
# simpBrm_carngroup_parastrans <- add_ic(simpBrm_carngroup_parastrans, ic = "kfold")
# saveRDS(simpBrm_carngroup_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_parastrans.RDS")
# 
# # model fits, predictions and residuals
# simpMu_carngroup_parastrans <- fitted(simpBrm_carngroup_parastrans)
# simpPredict_carngroup_parastrans <- predict(simpBrm_carngroup_parastrans)
# simpResid_carngroup_parastrans <- residuals(simpBrm_carngroup_parastrans, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_carngroup_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_parastrans_mu.RDS")
# saveRDS(simpPredict_carngroup_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_parastrans_predict.RDS")
# saveRDS(simpResid_carngroup_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_parastrans_resid.RDS")
# 
# # marginal effects
# simpMarg_carngroup_parastrans <- plot(marginal_effects(simpBrm_carngroup_parastrans), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_carngroup_parastrans, "./Data/JPEK/simple/simp_brm_carngroup_parastrans_me.RDS")
# 
# ### UNGULATES WITH PARASITE TRANSMISSION AS RESPONSE ----
# # ...simple model
# simpBrm_unggroup_parastrans <- brm(
#   data = fullDat_parastrans_unggroup,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_unggroup_parastrans)
# plot(simpBrm_unggroup_parastrans)
# pp_check(simpBrm_unggroup_parastrans, nsamples = 500)
# 
# # add information criteria
# simpBrm_unggroup_parastrans <- add_ic(simpBrm_unggroup_parastrans, ic = "loo", reloo = TRUE)
# simpBrm_unggroup_parastrans <- add_ic(simpBrm_unggroup_parastrans, ic = "kfold")
# saveRDS(simpBrm_unggroup_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_parastrans.RDS")
# 
# # model fits, predictions and residuals
# simpMu_unggroup_parastrans <- fitted(simpBrm_unggroup_parastrans)
# simpPredict_unggroup_parastrans <- predict(simpBrm_unggroup_parastrans)
# simpResid_unggroup_parastrans <- residuals(simpBrm_unggroup_parastrans, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_unggroup_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_parastrans_mu.RDS")
# saveRDS(simpPredict_unggroup_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_parastrans_predict.RDS")
# saveRDS(simpResid_unggroup_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_parastrans_resid.RDS")
# 
# # marginal effects
# simpMarg_unggroup_parastrans <- plot(marginal_effects(simpBrm_unggroup_parastrans), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_unggroup_parastrans, "./Data/JPEK/simple/simp_brm_unggroup_parastrans_me.RDS")
# 
# ### PRIMATES WITH PARASITE TRANSMISSION AS RESPONSE ----
# # ...simple model
# simpBrm_primgroup_parastrans <- brm(
#   data = fullDat_parastrans_primgroup,
#   family = binomial,
#   formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
#                  combIUCN*logNumHostCitations),
#   iter =4000, warmup = 2000, chains = 4, cores = 4,
#   control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
# 
# # quick checks
# summary(simpBrm_primgroup_parastrans)
# plot(simpBrm_primgroup_parastrans)
# pp_check(simpBrm_primgroup_parastrans, nsamples = 500)
# 
# # add information criteria
# simpBrm_primgroup_parastrans <- add_ic(simpBrm_primgroup_parastrans, ic = "loo", reloo = TRUE)
# simpBrm_primgroup_parastrans <- add_ic(simpBrm_primgroup_parastrans, ic = "kfold")
# saveRDS(simpBrm_primgroup_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_parastrans.RDS")
# 
# # model fits, predictions and residuals
# simpMu_primgroup_parastrans <- fitted(simpBrm_primgroup_parastrans)
# simpPredict_primgroup_parastrans <- predict(simpBrm_primgroup_parastrans)
# simpResid_primgroup_parastrans <- residuals(simpBrm_primgroup_parastrans, type = "pearson")
# 
# # save outputs
# saveRDS(simpMu_primgroup_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_parastrans_mu.RDS")
# saveRDS(simpPredict_primgroup_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_parastrans_predict.RDS")
# saveRDS(simpResid_primgroup_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_parastrans_resid.RDS")
# 
# # marginal effects
# simpMarg_primgroup_parastrans <- plot(marginal_effects(simpBrm_primgroup_parastrans), method = "fitted", plot = FALSE)
# saveRDS(simpMarg_primgroup_parastrans, "./Data/JPEK/simple/simp_brm_primgroup_parastrans_me.RDS")