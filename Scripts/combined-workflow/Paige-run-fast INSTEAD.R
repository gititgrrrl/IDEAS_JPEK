# Additional code for Paige to run on fast computers -- DON'T WORRY ABOUT GETTING LOO & KFOLD ON THE PARAS TRANS AND PARAS TYPE MODELS. I JUST REALIZED BECAUSE THOSE ARE NOT TRUNCATED POISSONS, I CAN RUN THEM ON MY OWN COMPUTER! (IT'S THE TRUNCATED POISSON THAT TAKES TOO MUCH MEMORY FOR ME TO RUN). SO PLEASE JUST RUN THE UNCOMMENTED CODE IN THIS SCRIPT.


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


# Get residuals for full and simple_full data models
# RICHNESS
simple_brm_all_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_fulldat.RDS") 
full_brm_all <- readRDS("./Data/JPEK/full/full_brm_all.RDS") 

simple_brm_all_resid_fulldat <- residuals(simple_brm_all_fulldat, type = "pearson")
full_brm_all_resid <- residuals(full_brm_all, type = "pearson")

saveRDS(simple_brm_all_resid_fulldat, "./Data/JPEK/simple/simple_brm_all_resid_fulldat.RDS")
saveRDS(full_brm_all_resid, "./Data/JPEK/full/full_brm_all_resid.RDS")

### -- GROUP SIZE SUB-ANALYSES ----
# # --- subset data for full model ALL PARASITES for CARNIVORES ----
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

### -- SIMPLE MODELS ----
# --- subset data for full model ALL PARASITES for CARNIVORES ----
simpBrm_carngroup <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all.RDS") # model
# add information criteria
simpBrm_carngroup <- add_ic(simpBrm_carngroup, ic = "loo", reloo = TRUE)
simpBrm_carngroup <- add_ic(simpBrm_carngroup, ic = "kfold")

# model fits, predictions and residuals
simpMu_carngroup <- fitted(simpBrm_carngroup)
simpPredict_carngroup <- predict(simpBrm_carngroup)
simpResid_carngroup <- residuals(simpBrm_carngroup, type = "pearson")

# save outputs
saveRDS(simpMu_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all_mu.RDS")
saveRDS(simpPredict_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all_predict.RDS")
saveRDS(simpResid_carngroup, "./Data/JPEK/simple/simp_brm_carngroup_all_resid.RDS")

# marginal effects
simp_carngroup_me <- plot(marginal_effects(simpBrm_carngroup), method = "fitted", plot = FALSE)
saveRDS(simp_carngroup_me, "./Data/JPEK/simple/simp_brm_carngroup_me.RDS")

# # # --- subset data for full model ALL PARASITES for UNGULATES ----
simpBrm_unggroup <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all.RDS") # model
# add information criteria
simpBrm_unggroup <- add_ic(simpBrm_unggroup, ic = "loo", reloo = TRUE)
simpBrm_unggroup <- add_ic(simpBrm_unggroup, ic = "kfold")

# model fits, predictions and residuals
simpMu_unggroup <- fitted(simpBrm_unggroup)
simpPredict_unggroup <- predict(simpBrm_unggroup)
simpResid_unggroup <- residuals(simpBrm_unggroup, type = "pearson")

# save outputs
saveRDS(simpMu_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all_mu.RDS")
saveRDS(simpPredict_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all_predict.RDS")
saveRDS(simpResid_unggroup, "./Data/JPEK/simple/simp_brm_unggroup_all_resid.RDS")

# marginal effects
simp_unggroup_me <- plot(marginal_effects(simpBrm_unggroup), method = "fitted", plot = FALSE)
saveRDS(simp_unggroup_me, "./Data/JPEK/simple/simp_brm_unggroup_me.RDS")

# # --- subset data for full model ALL PARASITES for PRIMATES ----
simpBrm_primgroup <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all.RDS") # model
# add information criteria
simpBrm_primgroup <- add_ic(simpBrm_primgroup, ic = "loo", reloo = TRUE)
simpBrm_primgroup <- add_ic(simpBrm_primgroup, ic = "kfold")

# model fits, predictions and residuals
simpMu_primgroup <- fitted(simpBrm_primgroup)
simpPredict_primgroup <- predict(simpBrm_primgroup)
simpResid_primgroup <- residuals(simpBrm_primgroup, type = "pearson")

# save outputs
saveRDS(simpMu_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all_mu.RDS")
saveRDS(simpPredict_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all_predict.RDS")
saveRDS(simpResid_primgroup, "./Data/JPEK/simple/simp_brm_primgroup_all_resid.RDS")

# marginal effects
simp_primgroup_me <- plot(marginal_effects(simpBrm_primgroup), method = "fitted", plot = FALSE)
saveRDS(simp_primgroup_me, "./Data/JPEK/simple/simp_brm_primgroup_me.RDS")