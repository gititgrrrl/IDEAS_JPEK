### PARASITE TRANSMISSION FOR MODELS INCLUDING GROUP SIZE
rm(list=ls())
allDat <- read_csv("Data/JPEK/script4.csv")

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

## FORMAT THE DATA----
#...carnivores----
fullDat_parastrans_carngroup <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizeCar, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139
fullDat_parastrans_carngroup <- fullDat_parastrans_carngroup[complete.cases(fullDat_parastrans_carngroup),] # only 75 records
dim(fullDat_parastrans_carngroup)
  
fullDat_parastrans_carngroup$groupSizeCar <- factor(fullDat_parastrans_carngroup$groupSizeCar, levels = c("non_group", "group"))
fullDat_parastrans_carngroup$combIUCN <- factor(fullDat_parastrans_carngroup$combIUCN, levels = c("not_threatened", "threatened"))
fullDat_parastrans_carngroup_c <- fullDat_parastrans_carngroup %>%
  mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE))
write_csv(fullDat_parastrans_carngroup_c, "./Data/JPEK/allDat_parastrans_carngroup.csv")

#...primates----
fullDat_parastrans_primgroup <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 130
fullDat_parastrans_primgroup <- fullDat_parastrans_primgroup[complete.cases(fullDat_parastrans_primgroup),] # only 87 records
dim(fullDat_parastrans_primgroup)

fullDat_parastrans_primgroup$combIUCN <- factor(fullDat_parastrans_primgroup$combIUCN, levels = c("not_threatened", "threatened"))
fullDat_parastrans_primgroup_c <- fullDat_parastrans_primgroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE))
write_csv(fullDat_parastrans_primgroup_c, "./Data/JPEK/allDat_parastrans_primgroup.csv")

#...ungulates----
fullDat_parastrans_unggroup <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 93
fullDat_parastrans_unggroup <- fullDat_parastrans_unggroup[complete.cases(fullDat_parastrans_unggroup),] # only 60 records
dim(fullDat_parastrans_unggroup)

fullDat_parastrans_unggroup$combIUCN <- factor(fullDat_parastrans_unggroup$combIUCN, levels = c("not_threatened", "threatened"))
fullDat_parastrans_unggroup_c <- fullDat_parastrans_unggroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE))
write_csv(fullDat_parastrans_unggroup_c, "./Data/JPEK/allDat_parastrans_unggroup.csv")

### RICHNESS OF CLOSELY TRANSMITTED PARASITES IS RESPONSE, FULL MODELS----
# ...carnivores----
full_parastrans_carngroup_mod <- brm(
  data = fullDat_parastrans_carngroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*groupSizeCar +
                 logHostSpeciesRange_c +
                 logHostMass_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(full_parastrans_carngroup_mod)
plot(full_parastrans_carngroup_mod)
pp_check(full_parastrans_carngroup_mod, nsamples = 500)

# add information criteria
full_parastrans_carngroup_mod <- add_ic(full_parastrans_carngroup_mod, ic = "loo", reloo = TRUE)
full_parastrans_carngroup_mod <- add_ic(full_parastrans_carngroup_mod, ic = "kfold")
saveRDS(full_parastrans_carngroup_mod, "./FINAL/full/full_parastrans_carngroup_mod.RDS")

# model fits, predictions and residuals
full_parastrans_carngroup_mod <- readRDS("./FINAL/full/full_parastrans_carngroup_mod.RDS")
full_parastrans_carngroup_fitted <- fitted(full_parastrans_carngroup_mod)
full_parastrans_carngroup_predict <- predict(full_parastrans_carngroup_mod)
full_parastrans_carngroup_resid <- residuals(full_parastrans_carngroup_mod, type = "pearson")

# save outputs
saveRDS(full_parastrans_carngroup_fitted, "./FINAL/full/full_parastrans_carngroup_fitted.RDS")
saveRDS(full_parastrans_carngroup_predict, "./FINAL/full/full_parastrans_carngroup_predict.RDS")
saveRDS(full_parastrans_carngroup_resid, "./FINAL/full/full_parastrans_carngroup_resid.RDS")

# marginal effects
full_parastrans_carngroup_marginal <- plot(marginal_effects(full_parastrans_carngroup_mod), method = "fitted", plot = FALSE)
saveRDS(full_parastrans_carngroup_marginal, "./FINAL/full/full_parastrans_carngroup_marginal.RDS")

# specifically for plotting groupsize-by-threat effects
full_parastrans_carngroup_marginal_groupsize <- marginal_effects(full_parastrans_carngroup_mod, effects = "groupSizeCar:combIUCN", method = "fitted")
saveRDS(full_parastrans_carngroup_marginal_groupsize, "./FINAL/full/full_parastrans_carngroup_marginal_groupsize.RDS")

#...primates
full_parastrans_primgroup_mod <- brm(
  data = fullDat_parastrans_primgroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logGroupSizePriUng_c +
                 logHostSpeciesRange_c +
                 logHostMass_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(full_parastrans_primgroup_mod)
plot(full_parastrans_primgroup_mod)
pp_check(full_parastrans_primgroup_mod, nsamples = 500)

# add information criteria
full_parastrans_primgroup_mod <- add_ic(full_parastrans_primgroup_mod, ic = "loo", reloo = TRUE)
full_parastrans_primgroup_mod <- add_ic(full_parastrans_primgroup_mod, ic = "kfold")
saveRDS(full_parastrans_primgroup_mod, "./FINAL/full/full_parastrans_primgroup_mod.RDS")

# model fits, predictions and residuals
full_parastrans_primgroup_mod <- readRDS("./FINAL/full/full_parastrans_primgroup_mod.RDS")
full_parastrans_primgroup_fitted <- fitted(full_parastrans_primgroup_mod)
full_parastrans_primgroup_predict <- predict(full_parastrans_primgroup_mod)
full_parastrans_primgroup_resid <- residuals(full_parastrans_primgroup_mod, type = "pearson")

# save outputs
saveRDS(full_parastrans_primgroup_fitted, "./FINAL/full/full_parastrans_primgroup_fitted.RDS")
saveRDS(full_parastrans_primgroup_predict, "./FINAL/full/full_parastrans_primgroup_predict.RDS")
saveRDS(full_parastrans_primgroup_resid, "./FINAL/full/full_parastrans_primgroup_resid.RDS")

# marginal effects
full_parastrans_primgroup_marginal <- plot(marginal_effects(full_parastrans_primgroup_mod), method = "fitted", plot = FALSE)
saveRDS(full_parastrans_primgroup_marginal, "./FINAL/full/full_parastrans_primgroup_marginal.RDS")

# specifically for plotting groupsize-by-threat effects
full_parastrans_primgroup_marginal_groupsize <- marginal_effects(full_parastrans_primgroup_mod, effects = "logGroupSizePriUng_c:combIUCN", method = "fitted")
saveRDS(full_parastrans_primgroup_marginal_groupsize, "./FINAL/full/full_parastrans_primgroup_marginal_groupsize.RDS")

#...ungulates
full_parastrans_unggroup_mod <- brm(
  data = fullDat_parastrans_unggroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logGroupSizePriUng_c +
                 logHostSpeciesRange_c +
                 logHostMass_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(full_parastrans_unggroup_mod)
plot(full_parastrans_unggroup_mod)
pp_check(full_parastrans_unggroup_mod, nsamples = 500)

# add information criteria
full_parastrans_unggroup_mod <- add_ic(full_parastrans_unggroup_mod, ic = "loo", reloo = TRUE)
full_parastrans_unggroup_mod <- add_ic(full_parastrans_unggroup_mod, ic = "kfold")
saveRDS(full_parastrans_unggroup_mod, "./FINAL/full/full_parastrans_unggroup_mod.RDS")

# model fits, predictions and residuals
full_parastrans_unggroup_mod <- readRDS("./FINAL/full/full_parastrans_unggroup_mod.RDS")
full_parastrans_unggroup_fitted <- fitted(full_parastrans_unggroup_mod)
full_parastrans_unggroup_predict <- predict(full_parastrans_unggroup_mod)
full_parastrans_unggroup_resid <- residuals(full_parastrans_unggroup_mod, type = "pearson")

# save outputs
saveRDS(full_parastrans_unggroup_fitted, "./FINAL/full/full_parastrans_unggroup_fitted.RDS")
saveRDS(full_parastrans_unggroup_predict, "./FINAL/full/full_parastrans_unggroup_predict.RDS")
saveRDS(full_parastrans_unggroup_resid, "./FINAL/full/full_parastrans_unggroup_resid.RDS")

# marginal effects
full_parastrans_unggroup_marginal <- plot(marginal_effects(full_parastrans_unggroup_mod), method = "fitted", plot = FALSE)
saveRDS(full_parastrans_unggroup_marginal, "./FINAL/full/full_parastrans_unggroup_marginal.RDS")

# specifically for plotting groupsize-by-threat effects
full_parastrans_unggroup_marginal_groupsize <- marginal_effects(full_parastrans_unggroup_mod, effects = "logGroupSizePriUng_c:combIUCN", method = "fitted")
saveRDS(full_parastrans_unggroup_marginal_groupsize, "./FINAL/full/full_parastrans_unggroup_marginal_groupsize.RDS")

### RICHNESS OF CLOSELY TRANSMITTED PARASITES IS RESPONSE, SIMPLE MODELS----
# ...carnivores----
simple_parastrans_carngroup_mod <- brm(
  data = fullDat_parastrans_carngroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(simple_parastrans_carngroup_mod)
plot(simple_parastrans_carngroup_mod)
pp_check(simple_parastrans_carngroup_mod, nsamples = 500)

# add information criteria
simple_parastrans_carngroup_mod <- add_ic(simple_parastrans_carngroup_mod, ic = "loo", reloo = TRUE)
simple_parastrans_carngroup_mod <- add_ic(simple_parastrans_carngroup_mod, ic = "kfold")
saveRDS(simple_parastrans_carngroup_mod, "./FINAL/simple/simple_parastrans_carngroup_mod.RDS")

# model fits, predictions and residuals
simple_parastrans_carngroup_mod <- readRDS("./FINAL/simple/simple_parastrans_carngroup_mod.RDS")
simple_parastrans_carngroup_fitted <- fitted(simple_parastrans_carngroup_mod)
simple_parastrans_carngroup_predict <- predict(simple_parastrans_carngroup_mod)
simple_parastrans_carngroup_resid <- residuals(simple_parastrans_carngroup_mod, type = "pearson")

# save outputs
saveRDS(simple_parastrans_carngroup_fitted, "./FINAL/simple/simple_parastrans_carngroup_fitted.RDS")
saveRDS(simple_parastrans_carngroup_predict, "./FINAL/simple/simple_parastrans_carngroup_predict.RDS")
saveRDS(simple_parastrans_carngroup_resid, "./FINAL/simple/simple_parastrans_carngroup_resid.RDS")

# marginal effects
simple_parastrans_carngroup_marginal <- plot(marginal_effects(simple_parastrans_carngroup_mod), method = "fitted", plot = FALSE)
saveRDS(simple_parastrans_carngroup_marginal, "./FINAL/simple/simple_parastrans_carngroup_marginal.RDS")

#...primates
simple_parastrans_primgroup_mod <- brm(
  data = fullDat_parastrans_primgroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(simple_parastrans_primgroup_mod)
plot(simple_parastrans_primgroup_mod)
pp_check(simple_parastrans_primgroup_mod, nsamples = 500)

# add information criteria
simple_parastrans_primgroup_mod <- add_ic(simple_parastrans_primgroup_mod, ic = "loo", reloo = TRUE)
simple_parastrans_primgroup_mod <- add_ic(simple_parastrans_primgroup_mod, ic = "kfold")
saveRDS(simple_parastrans_primgroup_mod, "./FINAL/simple/simple_parastrans_primgroup_mod.RDS")

# model fits, predictions and residuals
simple_parastrans_primgroup_mod <- readRDS("./FINAL/simple/simple_parastrans_primgroup_mod.RDS")
simple_parastrans_primgroup_fitted <- fitted(simple_parastrans_primgroup_mod)
simple_parastrans_primgroup_predict <- predict(simple_parastrans_primgroup_mod)
simple_parastrans_primgroup_resid <- residuals(simple_parastrans_primgroup_mod, type = "pearson")

# save outputs
saveRDS(simple_parastrans_primgroup_fitted, "./FINAL/simple/simple_parastrans_primgroup_fitted.RDS")
saveRDS(simple_parastrans_primgroup_predict, "./FINAL/simple/simple_parastrans_primgroup_predict.RDS")
saveRDS(simple_parastrans_primgroup_resid, "./FINAL/simple/simple_parastrans_primgroup_resid.RDS")

# marginal effects
simple_parastrans_primgroup_marginal <- plot(marginal_effects(simple_parastrans_primgroup_mod), method = "fitted", plot = FALSE)
saveRDS(simple_parastrans_primgroup_marginal, "./FINAL/simple/simple_parastrans_primgroup_marginal.RDS")

#...ungulates
simple_parastrans_unggroup_mod <- brm(
  data = fullDat_parastrans_unggroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(simple_parastrans_unggroup_mod)
plot(simple_parastrans_unggroup_mod)
pp_check(simple_parastrans_unggroup_mod, nsamples = 500)

# add information criteria
simple_parastrans_unggroup_mod <- add_ic(simple_parastrans_unggroup_mod, ic = "loo", reloo = TRUE)
simple_parastrans_unggroup_mod <- add_ic(simple_parastrans_unggroup_mod, ic = "kfold")
saveRDS(simple_parastrans_unggroup_mod, "./FINAL/simple/simple_parastrans_unggroup_mod.RDS")

# model fits, predictions and residuals
simple_parastrans_unggroup_mod <- readRDS("./FINAL/simple/simple_parastrans_unggroup_mod.RDS")
simple_parastrans_unggroup_fitted <- fitted(simple_parastrans_unggroup_mod)
simple_parastrans_unggroup_predict <- predict(simple_parastrans_unggroup_mod)
simple_parastrans_unggroup_resid <- residuals(simple_parastrans_unggroup_mod, type = "pearson")

# save outputs
saveRDS(simple_parastrans_unggroup_fitted, "./FINAL/simple/simple_parastrans_unggroup_fitted.RDS")
saveRDS(simple_parastrans_unggroup_predict, "./FINAL/simple/simple_parastrans_unggroup_predict.RDS")
saveRDS(simple_parastrans_unggroup_resid, "./FINAL/simple/simple_parastrans_unggroup_resid.RDS")

# marginal effects
simple_parastrans_unggroup_marginal <- plot(marginal_effects(simple_parastrans_unggroup_mod), method = "fitted", plot = FALSE)
saveRDS(simple_parastrans_unggroup_marginal, "./FINAL/simple/simple_parastrans_unggroup_marginal.RDS")