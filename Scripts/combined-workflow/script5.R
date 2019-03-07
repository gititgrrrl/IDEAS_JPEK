# fifth script in workflow 
# "SIMPLE MODELS": richness ~ hostGroup + threat status (yes no) + citations  
# broken down by transmission type (close non close) & parasite type (micro macroparasite)

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
allDat <- read_csv("Data/JPEK/script4.csv") # columns 28 & 29 = # host and parasite citations 

# --- subset data for simple model ALL PARASITE TYPES --- 

simpleDat <- allDat %>%
  select(hostName, parRich, numHostCitations, combIUCN, hostGroup) %>%
  distinct()
simpleDat <- data.frame(simpleDat[complete.cases(simpleDat),])                 # 408 records

# --- run simple model ALL PARASITE TYPES---

# MOD 1 - ZERO TRUNCATED POISSON. FULL MODEL ----
# > What we need is a one-inflated, zero-truncated Poisson!
# > NOTES ON MODEL 1

simpleBrm <- brm(
  data = simpleDat, 
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 numHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# add information criteria
simpleBrm <- add_ic(simpleBrm, ic = "loo", reloo = TRUE)  
simpleBrm <- add_ic(simpleBrm, ic = "kfold")

# model fits and predictions
simpleMu <- fitted(simpleBrm)
simplePredict <- predict(simpleBrm)

# --- subset data for simple model CLOSE PARASITES --- 

simpleDat_close <- allDat %>%
  select(hostName, parRich_close, numHostCitations, combIUCN, hostGroup) %>%
  filter(parRich_close > 0) %>%
  distinct()
simpleDat_close <- data.frame(simpleDat_close[complete.cases(simpleDat_close),])     # 300

# --- run simple model CLOSE PARASITES ---

simpleBrm_close <- brm(
  data = simpleDat_close, 
  family = poisson,
  formula = bf(parRich_close | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 numHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 8,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# add information criteria
simpleBrm_close <- add_ic(simpleBrm_close, ic = "loo", reloo = TRUE)  
simpleBrm_close <- add_ic(simpleBrm_close, ic = "kfold")

# model fits and predictions
simpleMu_close <- fitted(simpleBrm_close)
simplePredict_close <- predict(simpleBrm_close)

# --- subset data for simple model NON CLOSE PARASITES --- 

simpleDat_nonclose <- allDat %>%
  select(hostName, parRich, parRich_close, numHostCitations, combIUCN, hostGroup) %>%
  mutate(parRich_nonclose=(parRich-parRich_close)) %>% select(-parRich, -parRich_close)  %>%
  filter(parRich_nonclose > 0) %>%
  distinct()
simpleDat_nonclose <- simpleDat_nonclose[complete.cases(simpleDat_nonclose),]     # 367

# --- run simple model NON CLOSE PARASITES ---

simpleBrm_nonclose <- brm(
  data = simpleDat_nonclose, 
  family = poisson,
  formula = bf(parRich_nonclose | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 numHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 8,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# add information criteria
# simpleBrm_nonclose <- add_ic(simpleBrm_nonclose, ic = "loo", reloo = TRUE)  
# simpleBrm_nonclose <- add_ic(simpleBrm_nonclose, ic = "kfold")

# model fits and predictions
simpleMu_nonclose <- fitted(simpleBrm_nonclose)
simplePredict_nonclose <- predict(simpleBrm_nonclose)

# --- subset data for simple model MICROSPARASITES (bacteria, viruses, protozoa) --- 

simpleDat_micro <- allDat %>%
  select(hostName, numHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich) %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich)) %>% 
  select(-bacteriaRich, -virusRich, -protozoaRich)%>%
  filter(parRich_micro> 0) %>%
  distinct()

simpleDat_micro <- simpleDat_micro[complete.cases(simpleDat_micro),]    # 312 

# --- run simple model MICROSPARASITES---

simpleBrm_micro <- brm(
  data = simpleDat_micro, 
  family = poisson,
  formula = bf(parRich_micro | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 numHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 8,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# add information criteria
# simpleBrm_micro <- add_ic(simpleBrm_micro, ic = "loo", reloo = TRUE)  
# simpleBrm_micro <- add_ic(simpleBrm_micro, ic = "kfold")

# model fits and predictions
simpleMu_micro <- fitted(simpleBrm_micro)
simplePredict_micro <- predict(simpleBrm_micro)

# --- subset data for simple model MACROSPARASITES --- 

simpleDat_macro <- allDat %>%
  select(hostName, numHostCitations, combIUCN, hostGroup, helminthRich, arthropodRich) %>%
  mutate(parRich_macro=(helminthRich + arthropodRich)) %>% 
  select(-helminthRich, -arthropodRich)%>%
  filter(parRich_macro> 0) %>%
  distinct()
simpleDat_macro <- simpleDat_macro[complete.cases(simpleDat_macro),]    # 345

# --- run simple model MACROSPARASITES---

simpleBrm_macro <- brm(
  data = simpleDat_macro, 
  family = poisson,
  formula = bf(parRich_macro | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 numHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 8,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# add information criteria
# simpleBrm_macro <- add_ic(simpleBrm_macro, ic = "loo", reloo = TRUE)  
# simpleBrm_macro <- add_ic(simpleBrm_macro, ic = "kfold")

# model fits and predictions
simpleMu_macro <- fitted(simpleBrm_macro)
simplePredict_macro <- predict(simpleBrm_macro)

# --- 

saveRDS(simpleBrm, "./Data/JPEK/simple/simple_brm_all.RDS") # model of parasite spp richness total
saveRDS(simpleMu, "./Data/JPEK/simple/simple_brm_all_mu.RDS")
saveRDS(simplePredict, "./Data/JPEK/simple/simple_brm_all_predict.RDS")

saveRDS(simpleBrm_close, "./Data/JPEK/simple/simple_brm_close.RDS") #   model of closely transmitted parasite spp richness 
saveRDS(simpleMu_close, "./Data/JPEK/simple/simple_brm_close_mu.RDS")
saveRDS(simplePredict_close, "./Data/JPEK/simple/simple_brm_close_predict.RDS")

saveRDS(simpleBrm_nonclose, "./Data/JPEK/simple/simple_brm_nonclose.RDS") #   model of NON closely transmitted parasite spp richness 
saveRDS(simpleMu_nonclose, "./Data/JPEK/simple/simple_brm_nonclose_mu.RDS")
saveRDS(simplePredict_nonclose, "./Data/JPEK/simple/simple_brm_nonclose_predict.RDS")

saveRDS(simpleBrm_micro, "./Data/JPEK/simple/simple_brm_micro.RDS") #   model of micro spp richness 
saveRDS(simpleMu_micro, "./Data/JPEK/simple/simple_brm_micro_mu.RDS")
saveRDS(simplePredict_micro, "./Data/JPEK/simple/simple_brm_micro_predict.RDS")

saveRDS(simpleBrm_macro, "./Data/JPEK/simple/simple_brm_macro.RDS") #   model of macro spp richness
aveRDS(simpleMu_macro, "./Data/JPEK/simple/simple_brm_macro_mu.RDS")
saveRDS(simplePredict_macro, "./Data/JPEK/simple/simple_brm_macro_predict.RDS")
