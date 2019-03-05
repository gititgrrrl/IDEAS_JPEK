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
  select(hostName, parRich, numHostCitations, combIUCN, hostGroup) 
simpleDat <- simpleDat[complete.cases(simpleDat),]     # 18691 records

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
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 
simpleBrm <- add_ic(simpleBrm, ic = "loo", reloo = TRUE)  
simpleBrm <- add_ic(simpleBrm, ic = "kfold")

# --- subset data for simple model CLOSE PARASITES --- 

# --- run simple model CLOSE PARASITES ---

# --- subset data for simple model NON CLOSE PARASITES --- 

# --- run simple model NON CLOSE PARASITES---

# --- subset data for simple model MICROSPARASITES --- 

# --- run simple model MICROSPARASITES---

# --- subset data for simple model MACROSPARASITES --- 

# --- run simple model MACROSPARASITES---

# --- 

saveRDS(simpleBrm, "./Data/JPEK/simple_brm_all.RDS") # model of parasite spp richness total
saveRDS(simpleBrm_close, "./Data/JPEK/simple_brm_close.RDS") #   model of closely transmitted parasite spp richness 
saveRDS(simpleBrm_nonclose, "./Data/JPEK/simple_brm_nonclose.RDS") #   model of NON closely transmitted parasite spp richness 
saveRDS(simpleBrm_micro, "./Data/JPEK/simple_brm_micro.RDS") #   model of micro spp richness 
saveRDS(simpleBrm_macro, "./Data/JPEK/simple_brm_macro.RDS") #   model of macro spp richness

