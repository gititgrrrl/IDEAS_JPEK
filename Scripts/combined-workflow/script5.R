# fifth script in workflow 
# "SIMPLE MODELS": richness ~ hostGroup + threat status (yes no) + citations  
# broken down by transmission type (close non close) & parasite type (micro macroparasite)

# if running code first time:
rm(list=ls())

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

# --- subset data for simple model ALL PARASITES --- 

simpleDat <- allDat %>%
  select(hostName, parRich, logNumHostCitations, combIUCN, hostGroup) %>% # <<<<< (Ellen) grab the logNumHostCitations (not numHostCitations)
  distinct()
simpleDat <- simpleDat[complete.cases(simpleDat),]              # 362 hosts

### --- Simple model on simple data --- ###

simpleBrm <- brm(
  data = simpleDat, 
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN * hostGroup + 
                 combIUCN * logNumHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# quick checks
summary(simpleBrm)
plot(simpleBrm)
pp_check(simpleBrm, nsamples = 500)

# add information criteria
simpleBrm <- add_ic(simpleBrm, ic = "loo", reloo = TRUE)
simpleBrm <- add_ic(simpleBrm, ic = "kfold")

# model fits and predictions
simpleMu <- fitted(simpleBrm)
simplePredict <- predict(simpleBrm)

# --- subset data for simple model PARASITE TRANSMISSION --- 

# (Ellen comment--I think it should be this... use the parRichCloseOnly in a binomial regression, log citations interact with IUCN group
simpleDat_parastrans <- allDat %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, logNumHostCitations, combIUCN, hostGroup) %>%
  distinct() %>%
  filter(parRichTransKnown > 0)

simpleDat_parastrans <- simpleDat_parastrans[complete.cases(simpleDat_parastrans),]     # 333

# --- run simple model PARASITE TRANSMISSION  ---

# BINOMIAL regression--only citations and its interaction with IUCNstatus is significant
simpleBrm_parastrans <- brm(
  data = simpleDat_parastrans,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN * hostGroup + 
                 combIUCN * logNumHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# quick checks
summary(simpleBrm_parastrans)
plot(simpleBrm_parastrans)
pp_check(simpleBrm_parastrans, nsamples = 500)

# add information criteria
simpleBrm_parastrans <- add_ic(simpleBrm_parastrans, ic = "loo", reloo = TRUE)
simpleBrm_parastrans <- add_ic(simpleBrm_parastrans, ic = "kfold")

# model fits and predictions
simpleMu_parastrans <- fitted(simpleBrm_parastrans)
simplePredict_parastrans <- predict(simpleBrm_parastrans)

# simpleDat_close <- allDat %>%
#   select(hostName, parRich_close, numHostCitations, combIUCN, hostGroup) %>%
#   filter(parRich_close > 0) %>%
#   distinct()
# simpleDat_close <- data.frame(simpleDat_close[complete.cases(simpleDat_close),])     # 251
# 
# # --- run simple model CLOSE PARASITES ---
# 
# simpleBrm_close <- brm(
#   data = simpleDat_close, 
#   family = poisson,
#   formula = bf(parRich_close | trunc(lb = 1) ~
#                  combIUCN * hostGroup + 
#                  numHostCitations),
#   iter = 4000, warmup = 2000, chains = 4, cores = 8,
#   control = list(adapt_delta = .8, max_treedepth = 10)) 

# # --- subset data for simple model NON CLOSE PARASITES --- 
# 
# simpleDat_nonclose <- allDat %>%
#   select(hostName, parRich, parRich_close, numHostCitations, combIUCN, hostGroup) %>%
#   mutate(parRich_nonclose=(parRich-parRich_close)) %>% select(-parRich, -parRich_close)  %>%
#   filter(parRich_nonclose > 0) %>%
#   distinct()
# simpleDat_nonclose <- simpleDat_nonclose[complete.cases(simpleDat_nonclose),]     # 334
# 
# # --- run simple model NON CLOSE PARASITES ---
# 
# simpleBrm_nonclose <- brm(
#   data = simpleDat_nonclose, 
#   family = poisson,
#   formula = bf(parRich_nonclose | trunc(lb = 1) ~
#                  combIUCN * hostGroup + 
#                  numHostCitations),
#   iter = 4000, warmup = 2000, chains = 4, cores = 8,
#   control = list(adapt_delta = .8, max_treedepth = 10)) 
# 
# # add information criteria
# # simpleBrm_nonclose <- add_ic(simpleBrm_nonclose, ic = "loo", reloo = TRUE)  
# # simpleBrm_nonclose <- add_ic(simpleBrm_nonclose, ic = "kfold")
# 
# # model fits and predictions
# simpleMu_nonclose <- fitted(simpleBrm_nonclose)
# simplePredict_nonclose <- predict(simpleBrm_nonclose)

# --- subset data for simple model PARASITE TYPE

simpleDat_parastype <- allDat %>%
  select(hostName, logNumHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich))

simpleDat_parastype <- simpleDat_parastype[complete.cases(simpleDat_parastype),]    # 362 

# --- run simple model PARASITE TYPE ---

# BINOMIAL regression--only citations and its interaction with IUCNstatus is significant
simpleBrm_parastype <- brm(
  data = simpleDat_parastype,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN * hostGroup + 
                 combIUCN * logNumHostCitations),
  iter = 4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

# quick checks
summary(simpleBrm_parastype)
plot(simpleBrm_parastype)
pp_check(simpleBrm_parastype, nsamples = 500)

# add information criteria
simpleBrm_parastype <- add_ic(simpleBrm_parastype, ic = "loo", reloo = TRUE)
simpleBrm_parastype <- add_ic(simpleBrm_parastype, ic = "kfold")

# model fits and predictions
simpleMu_parastype <- fitted(simpleBrm_parastype)
simplePredict_parastype <- predict(simpleBrm_parastype)

# simpleDat_micro <- allDat %>%
#   select(hostName, numHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich) %>%
#   mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich)) %>% 
#   select(-bacteriaRich, -virusRich, -protozoaRich)%>%
#   filter(parRich_micro> 0) %>%
#   distinct()
# 
# simpleDat_micro <- simpleDat_micro[complete.cases(simpleDat_micro),]    # 298 
# 
# # --- run simple model MICROSPARASITES---
# 
# simpleBrm_micro <- brm(
#   data = simpleDat_micro, 
#   family = poisson,
#   formula = bf(parRich_micro | trunc(lb = 1) ~
#                  combIUCN * hostGroup + 
#                  numHostCitations),
#   iter = 4000, warmup = 2000, chains = 4, cores = 8,
#   control = list(adapt_delta = .8, max_treedepth = 10)) 
# 
# # add information criteria
# # simpleBrm_micro <- add_ic(simpleBrm_micro, ic = "loo", reloo = TRUE)  
# # simpleBrm_micro <- add_ic(simpleBrm_micro, ic = "kfold")
# 
# # model fits and predictions
# simpleMu_micro <- fitted(simpleBrm_micro)
# simplePredict_micro <- predict(simpleBrm_micro)
# 
# # --- subset data for simple model MACROSPARASITES --- 
# 
# simpleDat_macro <- allDat %>%
#   select(hostName, numHostCitations, combIUCN, hostGroup, helminthRich, arthropodRich) %>%
#   mutate(parRich_macro=(helminthRich + arthropodRich)) %>% 
#   select(-helminthRich, -arthropodRich)%>%
#   filter(parRich_macro> 0) %>%
#   distinct()
# simpleDat_macro <- simpleDat_macro[complete.cases(simpleDat_macro),]    # 291
# 
# # --- run simple model MACROSPARASITES---
# 
# simpleBrm_macro <- brm(
#   data = simpleDat_macro, 
#   family = poisson,
#   formula = bf(parRich_macro | trunc(lb = 1) ~
#                  combIUCN * hostGroup + 
#                  numHostCitations),
#   iter = 4000, warmup = 2000, chains = 4, cores = 8,
#   control = list(adapt_delta = .8, max_treedepth = 10)) 
# 
# # add information criteria
# # simpleBrm_macro <- add_ic(simpleBrm_macro, ic = "loo", reloo = TRUE)  
# # simpleBrm_macro <- add_ic(simpleBrm_macro, ic = "kfold")
# 
# # model fits and predictions
# simpleMu_macro <- fitted(simpleBrm_macro)
# simplePredict_macro <- predict(simpleBrm_macro)
# 
# # --- 
# 
saveRDS(simpleBrm, "./Data/JPEK/simple/simple_brm_all.RDS") # model of parasite spp richness total
saveRDS(simpleMu, "./Data/JPEK/simple/simple_brm_all_mu.RDS")
saveRDS(simplePredict, "./Data/JPEK/simple/simple_brm_all_predict.RDS")

saveRDS(simpleBrm_parastrans, "./Data/JPEK/simple/simple_brm_parastrans.RDS") # model of parasite transmission mode
saveRDS(simpleMu_parastrans, "./Data/JPEK/simple/simple_brm_parastrans_mu.RDS")
saveRDS(simplePredict_parastrans, "./Data/JPEK/simple/simple_brm_parastrans_predict.RDS")

saveRDS(simpleBrm_parastype, "./Data/JPEK/simple/simple_brm_parastype.RDS") # model of parasite type
saveRDS(simpleMu_parastype, "./Data/JPEK/simple/simple_brm_parastype_mu.RDS")
saveRDS(simplePredict_parastype, "./Data/JPEK/simple/simple_brm_parastype_predict.RDS")
# 
# saveRDS(simpleBrm_close, "./Data/JPEK/simple/simple_brm_close.RDS") #   model of closely transmitted parasite spp richness 
# saveRDS(simpleMu_close, "./Data/JPEK/simple/simple_brm_close_mu.RDS")
# saveRDS(simplePredict_close, "./Data/JPEK/simple/simple_brm_close_predict.RDS")
# 
# saveRDS(simpleBrm_nonclose, "./Data/JPEK/simple/simple_brm_nonclose.RDS") #   model of NON closely transmitted parasite spp richness 
# saveRDS(simpleMu_nonclose, "./Data/JPEK/simple/simple_brm_nonclose_mu.RDS")
# saveRDS(simplePredict_nonclose, "./Data/JPEK/simple/simple_brm_nonclose_predict.RDS")
# 
# saveRDS(simpleBrm_micro, "./Data/JPEK/simple/simple_brm_micro.RDS") #   model of micro spp richness 
# saveRDS(simpleMu_micro, "./Data/JPEK/simple/simple_brm_micro_mu.RDS")
# saveRDS(simplePredict_micro, "./Data/JPEK/simple/simple_brm_micro_predict.RDS")
# 
# saveRDS(simpleBrm_macro, "./Data/JPEK/simple/simple_brm_macro.RDS") #   model of macro spp richness
# saveRDS(simpleMu_macro, "./Data/JPEK/simple/simple_brm_macro_mu.RDS")
# saveRDS(simplePredict_macro, "./Data/JPEK/simple/simple_brm_macro_predict.RDS")
# 
simple_me <- plot(marginal_effects(simpleBrm), method = "fitted", plot = FALSE)
simple_me_parastrans <- plot(marginal_effects(simpleBrm_parastrans), method = "fitted", plot = FALSE)
simple_me_parastype <- plot(marginal_effects(simpleBrm_parastype), method = "fitted", plot = FALSE)

# simple_me_close <- plot(marginal_effects(simpleBrm_close), method = "fitted", plot = FALSE)
# simple_me_nonclose <- plot(marginal_effects(simpleBrm_nonclose), method = "fitted", plot = FALSE)
# simple_me_mic <- plot(marginal_effects(simpleBrm_micro), method = "fitted", plot = FALSE)
# simple_me_mac <- plot(marginal_effects(simpleBrm_macro), method = "fitted", plot = FALSE)
# 
# saveRDS(simple_me, "./Data/JPEK/simple/simple_brm_me.RDS")
# saveRDS(simple_me_close, "./Data/JPEK/simple/simple_brm_me_close.RDS")
# saveRDS(simple_me_nonclose, "./Data/JPEK/simple/simple_brm_me_nonclose.RDS")
# saveRDS(simple_me_mic, "./Data/JPEK/simple/simple_brm_me_micro.RDS")
# saveRDS(simple_me_mac, "./Data/JPEK/simple/simple_brm_me_macro.RDS")

saveRDS(simple_me, "./Data/JPEK/simple/simple_brm_me.RDS")
saveRDS(simple_me_parastrans, "./Data/JPEK/simple/simple_brm_parastrans_me.RDS")
saveRDS(simple_me_parastype, "./Data/JPEK/simple/simple_brm_parastype_me.RDS")