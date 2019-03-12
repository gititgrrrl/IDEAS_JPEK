# Model diagnostics
library(magrittr)
library(cowplot)
library(tidyverse)
library(rstan)
library(brms)
library(broom)
library(tidybayes)
library(purrr)

# Richness
simpleBrm <- readRDS("./Data/JPEK/simple/simple_brm_all.RDS") 
simpleBrm_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_fulldat.RDS") 
fullBrm <- readRDS("./Data/JPEK/full/full_brm_all.RDS") 

# Parasite transmission 
simpleBrm_parastrans <- readRDS("./Data/JPEK/simple/simple_brm_parastrans.RDS") 
simpleBrm_parastrans_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastrans_fulldat.RDS") 
fullBrm_parastrans <- readRDS("./Data/JPEK/full/full_brm_parastrans.RDS") 

# Parasite type
simpleBrm_parastype <- readRDS("./Data/JPEK/simple/simple_brm_parastype.RDS") 
simpleBrm_parastype_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastype_fulldat.RDS") 
fullBrm_parastype <- readRDS("./Data/JPEK/full/full_brm_parastype.RDS") 

########### RICHNESS #############

### --- fit of simple model of richness --- ###

summary(simpleBrm)

### --- fit of full model of richness --- ###

### --- fit of simple model of richness (on full dataset) --- ###

########### PROP CLOSE #############

### --- fit of simple model of prop close --- ###

### --- fit of full model of prop close --- ###

### --- fit of simple model of prop close (on full dataset) --- ###

########### PROP MICRO #############

### --- fit of simple model of prop micro --- ###

### --- fit of full model of prop micro --- ###

### --- fit of simple model of prop micro (on full dataset) --- ###