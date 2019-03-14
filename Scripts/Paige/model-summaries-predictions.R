# Model observations and predictions

library(magrittr)
library(cowplot)
library(tidyverse)
library(rstan)
library(brms)
library(broom)
library(tidybayes)
library(purrr)

### --- DATA --- ###
allDat <- read_csv("./Data/JPEK/script4.csv") # columns 28 & 29 = # host and parasite citations 

# SIMPLE RICHNESS
simpleDat <- allDat %>%
  select(hostName, parRich, logNumHostCitations, combIUCN, hostGroup) %>% # <<<<< (Ellen) grab the logNumHostCitations (not numHostCitations)
  distinct()
simpleDat <- simpleDat[complete.cases(simpleDat),]              # 362 hosts


# FULL RICHNESS
fullDat <- allDat %>%
  select(hostName, parRich, logNumHostCitations, combIUCN, 
         hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() # 392 records
fullDat <- fullDat[complete.cases(fullDat),] # only 232 records (before, it was 255 records--not sure why it's fewer now <<<<< ?????)

# FULL PARAS TRANS
fullDat_parastrans <- allDat %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, logNumHostCitations, combIUCN, hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 362

fullDat_parastrans <- fullDat_parastrans[complete.cases(fullDat_parastrans),]   # 220

# FULL PARA TYPE
fullDat_parastype <- allDat %>%
  select(hostName, logNumHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 393

fullDat_parastype <- fullDat_parastype[complete.cases(fullDat_parastype),] # 232

### --- MODELS --- ###
# Total richness
simple_mod_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_fulldat.RDS")
full_mod <- readRDS("./Data/JPEK/full/full_brm_all.RDS")

# Parasite transmission
simple_mod_parastrans_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastrans_fulldat.RDS")
full_mod_parastrans <- readRDS("./Data/JPEK/full/full_brm_parastrans.RDS")

# Parasite type
simple_mod_parastype_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastype_fulldat.RDS")
full_mod_parastype <- readRDS("./Data/JPEK/full/full_brm_parastype.RDS")

### --- MODEL PARAMETERS --- ###

simpleMu <- readRDS("./Data/JPEK/simple/simple_brm_all_mu.RDS")
simpleMu_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_mu_fulldat.RDS")
fullMu <- readRDS("./Data/JPEK/full/full_brm_all_mu.RDS")

### -- MODEL SUMMARIES --- ###
# EXAMINE COEFFICIENTS ----
# NOTES:
# > For total richness response, model coefficients should be exponentiated to get the estimated richness
# > For parasite transmission and type responses, model coefficients are log-odds. Exponentiate the coefficient to get odds. To get probability, use the equation: exp(coef) / (1 + exp(coef)). The function 'plogis' calculates the latter, so you just use plogis(coef)

# Response = total richness
fixef(simple_mod_fulldat)
fixef(full_mod)

# Response = % close parasite transmission
# NOTES:
# > The overall average % close in the real data is 15%. Mean parasite richness of close only is almost 2
# INTERESTING FINDINGS:
# > The more citations, the smaller the % close parasites--Therefore, as a host species is studied more, then the relative number of non-close parasites discovered increases at a much faster rate than close parasites. 
fixef(simple_mod_parastrans_fulldat)
fixef(full_mod_parastrans) # <<<<<<<<<<<<<<<<<<< ELLEN, PICK UP FROM HERE--WHAT IS MARGINAL EFFECTS ACTUALLY CALCULATING? COULD BE THE EXPECTED NUMBER OF CLOSE PARASITES GIVEN AVERAGE LOG CITATIONS, BUT THEN IT'S A BIT TOO LOW

# EXAMINE COEFFICIENT PLOTS ----
# NOTE: Only look at coefficients that both simple and full models estimate.
# UPSHOT: For all three response variables, the full vs. simple (using same data) coefficients are pretty different, especially for primates. But the confidence intervals for full model coefficients are quite large, so often still overlap the simple model coefficients.
FuncPlotCoef <- function(mod_list, terms_vec = NULL, exclude_terms = FALSE) {
  # Coefficient plots, comparing across models
  #
  # Args:
  #   mod_list:  List of (named) fitted regression models
  #   terms_vec: The model terms to plot. To plot ALL terms, leave as NULL
  #   exclude_terms: TRUE if parameters in 'terms_vec' should be excluded from the plot, rather than included
  # Returns:
  #   For multiple models, plot of coefficient estimates with 95% credible intervals. Each term is in a separate facet.
  # 
  df <- mod_list %>% purrr::map(function(.) { 
    broom::tidy(., conf.int = TRUE, par_type = "non-varying") }) %>% 
    dplyr::bind_rows(.id = "model") 
  if(exclude_terms) {
    df <- df %>% 
      filter(!term %in% terms_vec)
  } else {
    if(!is.null(terms_vec)) {
      df <- df %>%
        filter(term %in% terms_vec)
    }
  }
  p <- df %>% 
    ggplot(aes(term, estimate, ymin = lower, ymax = upper, color = model)) +
    geom_pointrange(position = position_dodge(width = 0.3)) + coord_flip() +
    labs(y = "Coefficient estimate (95% CI)", x = "Model parameter") +
    geom_hline(yintercept = 0,  linetype = "dotted") +
    theme_bw() +
    theme(legend.position="top")
  return(p)
}

# Response = total richness
FuncPlotCoef(mod_list = list(simple_fulldat = simple_mod_fulldat, full = full_mod), terms_vec = c("Intercept", "combIUCNthreatened", "hostGroupprimates", "hostGroupungulates", "logNumHostCitations", "combIUCNthreatened:hostGroupprimates", "combIUCNthreatened:hostGroupungulates", "combIUCNthreatened:logNumHostCitations"))

# Response = % close parasite transmission
FuncPlotCoef(mod_list = list(simple_parastrans_fulldat = simple_mod_parastrans_fulldat, full_parastrans = full_mod_parastrans), terms_vec = c("Intercept", "combIUCNthreatened", "hostGroupprimates", "hostGroupungulates", "logNumHostCitations", "combIUCNthreatened:hostGroupprimates", "combIUCNthreatened:hostGroupungulates", "combIUCNthreatened:logNumHostCitations"))

# Response = % micro parasite
FuncPlotCoef(mod_list = list(simple_parastype_fulldat = simple_mod_parastype_fulldat, full_parastype = full_mod_parastype), terms_vec = c("Intercept", "combIUCNthreatened", "hostGroupprimates", "hostGroupungulates", "logNumHostCitations", "combIUCNthreatened:hostGroupprimates", "combIUCNthreatened:hostGroupungulates", "combIUCNthreatened:logNumHostCitations"))

### -- MODEL PREDICTIONS --- ###

simplePred <- readRDS("./Data/JPEK/simple/simple_brm_all_predict.RDS")
simplePred_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_predict_fulldat.RDS")
fullPred <- readRDS("./Data/JPEK/full/full_brm_all_predict.RDS")

### --- MODEL OBSERVATIONS PREDICTIONS --- ###

