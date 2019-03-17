### FIGURES FOR DRAFT PAPER

rm(list=ls())

### LOAD PACKAGES ----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
}

best_complex <- data.frame(
  
### FIG.1 ----
# Effects of threat status on response (richness, % direct transmission, % micro-parasites) for simple vs. full models, separately for each host group
# READ IN DATA
margRich <- margParasTrans <- margParasType <- list()
 # >>> JUST GRAB THE DATA
### Marginal Effects ----
# Parasite richness
margRich$simple <- 
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

# CREATE FIGURE
fig1_R
plot(marg_list[[i]] +
       labs(y = "Parasite transmission") +
       scale_color_brewer(palette = "Dark2") +
       scale_fill_brewer(palette = "Dark2") +
       theme_bw(base_size = 10))
  
  