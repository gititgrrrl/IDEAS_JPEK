### SCRIPT FOR MODEL SUMMARIES INCL. MARGINAL EFFECTS

### QUESTIONS/COMMENTS ----
# > Should have run an interaction between IUCN status and ALL predictors! That means addin ght interaction for mass, lifespan, mean latitude

rm(list=ls())

### LOAD PACKAGES ----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "rstan", "brms", "broom", "tidybayes", "purrr", "glmmTMB")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run chains on multiple cores

### READ IN DATA ----
allDat <- read_csv("Data/JPEK/script4.csv")
# ...carnivores
dat_rich_carn <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() # 148 records
dat_rich_carn <- dat_rich_carn[complete.cases(dat_rich_carn),] # only 74 records
# ...ungulates
dat_rich_ung <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() # 102 records
dat_rich_ung <- dat_rich_ung[complete.cases(dat_rich_ung),] # only 60 records
# ...primates
dat_rich_prim <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() # 142 records
dat_rich_prim <- dat_rich_prim[complete.cases(dat_rich_prim),] # only 73 records

### READ IN MODELS ----
# Parasite richness
# ...full
mod_rich_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_all.RDS")
mod_rich_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_all.RDS")
mod_rich_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_all.RDS") 
# ...simple
mod_rich_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all.RDS")
mod_rich_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all.RDS") 
mod_rich_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all.RDS")

### READ IN FITTED ----
# Parasite richness
# ...full
fitted_rich_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_all_mu.RDS")
fitted_rich_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_all_mu.RDS")
fitted_rich_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_all_mu.RDS")
# ...simple
fitted_rich_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all_mu.RDS") 
fitted_rich_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all_mu.RDS") 
fitted_rich_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all_mu.RDS") 

### READ IN PREDICTED ----
# Parasite richness
# ...full
predict_rich_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_all_predict.RDS")
predict_rich_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_all_predict.RDS")
predict_rich_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_all_predict.RDS")
# ...simple
predict_rich_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all_predict.RDS")
predict_rich_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all_predict.RDS")
predict_rich_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all_predict.RDS")

### READ IN RESIDUALS ----
# Parasite richness
# ...full
resid_rich_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_all_resid.RDS")
resid_rich_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_all_resid.RDS")
resid_rich_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_all_resid.RDS")
# ...simple
resid_rich_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all_resid.RDS")
resid_rich_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all_resid.RDS")
resid_rich_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all_resid.RDS")

### Marginal Effects ----
# Parasite richness
# ...full
marginal_rich_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_me.RDS")
marginal_rich_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_me.RDS")
marginal_rich_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_me.RDS")
# ...simple
marginal_rich_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_me.RDS")
marginal_rich_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_me.RDS")
marginal_rich_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_me.RDS")

# Parasite transmission
# ...allInt
marginal_parastrans_allInt_prim <- readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_parastrans_me.RDS")
# ...full
marginal_parastrans_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_parastrans_me.RDS")
# ...simple
marginal_parastrans_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_parastrans_me.RDS")
marginal_parastrans_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_parastrans_me.RDS")
marginal_parastrans_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_parastrans_me.RDS")

# Parasite type
# ...allInt
marginal_parastype_allInt_carn <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_parastype_me.RDS")
marginal_parastype_allInt_prim <- readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_parastype_me.RDS")
# ...full
marginal_parastype_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_parastype_me.RDS")
# ...simple
marginal_parastype_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_parastype_me.RDS")
marginal_parastype_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_parastype_me.RDS")
marginal_parastype_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_parastype_me.RDS")

# # FULL PARAS TRANS
# fullDat_parastrans <- allDat %>%
#   select(hostName, parRichCloseOnly, parRichTransKnown, logNumHostCitations, combIUCN, hostGroup, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() %>%
#   filter(parRichTransKnown > 0) # 362
# 
# fullDat_parastrans <- fullDat_parastrans[complete.cases(fullDat_parastrans),]   # 220
# 
# # FULL PARA TYPE
# fullDat_parastype <- allDat %>%
#   select(hostName, logNumHostCitations, combIUCN, hostGroup, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
#   distinct() %>%
#   mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
#          parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 393
# 
# fullDat_parastype <- fullDat_parastype[complete.cases(fullDat_parastype),] # 232

### EXAMINE COEFFICIENTS ----
# NOTES:
# > For total richness response, model coefficients should be exponentiated to get the estimated richness
# > For parasite transmission and type responses, model coefficients are log-odds. Exponentiate the coefficient to get odds. To get probability, use the equation: exp(coef) / (1 + exp(coef)). The function 'plogis' calculates the latter, so you just use plogis(coef)
# 
# RESULTS ----
# ...carnivores
# > *(+ CITATIONS) PR INCREASES as host CITATIONS increases. The interaction between host citations & IUCN status is not significant.
# > *(INTERACTION IUCN:SPECIES RANGE) There is a significant interaction between SPECIES RANGE and IUCN STATUS for parasite richness. Therefore, we can't say <<<<<< 
# > *(INTERACTION IUCN:GROUP SIZE) There is a significant interaction between GROUP SIZE and IUCN STATUS for parasite richness. <<<<<Note that for carnivores the data only distinguish group sizes of 1,2, and >2.
# > *(- BODY MASS) larger species (larger BODY MASS) have LOWER PR
# > *(+ LIFESPAN) species that have longer LIFESPAN have HIGHER PR
# > NS(LATITUDE) absolute mean LATITUDE is NOT IMPORTANT

# ...ungulates
# > *(INTERACTION IUCN:CITATIONS) There is a significant interaction between HOST CITATIONS and IUCN STATUS for parasite richness. Therefore, we can't say <<<<<< 
# > *(INTERACTION IUCN:SPECIES RANGE) There is a significant interaction between SPECIES RANGE and IUCN STATUS for parasite richness. Therefore, we can't say <<<<<< 
# > *(- GROUP SIZE) As GROUP SIZE increase, PR DECLINES
# > NS(BODY MASS) BODY MASS is NOT IMPORTANT
# > *(- LIFESPAN) species that have longer LIFESPAN have LOWER PR
# > *(- LATITUDE) species FARTHER FROM EQUATOR have LOWER PR

# ...primates
# > *(+ CITATIONS) Host species with more CITATIONS have HIGHER PR
# > *(INTERACTION IUCN:SPECIES RANGE) There is a significant interaction between SPECIES RANGE and IUCN STATUS for parasite richness. Therefore, we can't say <<<<<< 
# > NS(GROUP SIZE) GROUP SIZE is NOT IMPORTANT
# > *(+ BODY MASS) species with larger BODY MASS have HIGHER PR
# > *(+ LIFESPAN) species that have longer LIFESPAN have HIGHER PR
# > NS(LATITUDE) absolute mean LATITUDE is NOT IMPORTANT

# ...simple models
# > For CARNIVORES, PR richness did not differ by threat status
# > For UNGULATES, threatened species had lower PR
# > For PRIMATES, 
# modlist <- list(carnivores = mod_rich_full_carn, ungulates = mod_rich_full_ung, primates = mod_rich_full_prim)
modlist <- list(carnivores = mod_rich_simple_carn, ungulates = mod_rich_simple_ung, primates = mod_rich_simple_prim)
# Parasite richness
lapply(modlist_simple, fixef)

# GENERATE COEFFICIENT PLOTS ----
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

### MARGINAL EFFECTS ----
marg_list <- marginal_parastype_full_ung # <<<<<<< USER SET
group_nam <- "ung" # <<<<<<< USER SET

for(i in 1:length(marg_list)) {
pdf(paste0("./Results/model-summaries-group/marginal_parastype_full/marginal_parastype_full_", group_nam, "_", names(marg_list)[i], ".pdf")) # <<<<<<< USER SET
plot(marg_list[[i]] +
  labs(y = "Parasite type") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 10))
dev.off()
}
