### For all best models, calculate 95% CI for threatened species

rm(list=ls())

### LOAD PACKAGES----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "plyr", "brms")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})

# Data frames of best models only----
# Parasite richness
modRichBest_list <- list(
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_all.RDS"),
  primates = readRDS("./Data/JPEK/full/full_brm_primgroup_c_all.RDS"),
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_all.RDS"))

# Parasite transmission
modParasTransBest_list <- list( 
  carnivores = readRDS("./Data/JPEK/full/full_brm_carngroup_c_parastrans.RDS"),
  primates = readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_c_parastrans.RDS"),
  ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_parastrans.RDS"))

# Parasite type
modParasTypeBest_list <- list(
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_ParasType.RDS"),
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_ParasType.RDS"),
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_ParasType.RDS"))
