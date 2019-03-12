# Model diagnostics

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

### --- summary and fit simple model of richness --- ###
summary(simpleBrm)

for(g in c("carnivores", "ungulates", "primates")) {
  pdf(paste0("./Results/model-diagnostics/richness_dens_plot_simple_", g, ".pdf")) 
  plot(pp_check(simpleBrm, type = "dens_overlay", nsamples = 300, 
                newdata = subset(simpleDat, hostGroup == g)) +
         ggtitle(paste0(g, " density plot")) +
         theme(plot.title = element_text(size = 12, face = "bold")))
  dev.off()
}
pdf("./Results/model-diagnostics/richness_dens_plot_simple.pdf") 
pp_check(simpleBrm, type = "violin_grouped", nsamples = 300, group = "combIUCN")
dev.off()

### --- summary and fit of simple model of richness (on full dataset) --- ###
summary(simpleBrm_fulldat)
for(g in c("carnivores", "ungulates", "primates")) {
  pdf(paste0("./Results/model-diagnostics/richness_dens_plot_simple_fulldat_", g, ".pdf")) 
  plot(pp_check(simpleBrm_fulldat, type = "dens_overlay", nsamples = 300, 
                newdata = subset(fullDat, hostGroup == g)) +
         ggtitle(paste0(g, " density plot")) +
         theme(plot.title = element_text(size = 12, face = "bold")))
  dev.off()
}
pdf("./Results/model-diagnostics/richness_dens_plot_simple_fulldat_iucn.pdf")
pp_check(simpleBrm_fulldat, type = "violin_grouped", nsamples = 300, group = "combIUCN")
dev.off()


### --- summary and fit of full model of richness --- ###
summary(fullBrm)
for(g in c("carnivores", "ungulates", "primates")) {
  pdf(paste0("./Results/model-diagnostics/richness_dens_plot_full_", g, ".pdf")) 
  plot(pp_check(fullBrm, type = "dens_overlay", nsamples = 300, 
                newdata = subset(fullDat, hostGroup == g)) +
         ggtitle(paste0(g, " density plot")) +
         theme(plot.title = element_text(size = 12, face = "bold")))
  dev.off()
}
pdf("./Results/model-diagnostics/richness_dens_plot_full_iucn.pdf")
pp_check(fullBrm, type = "violin_grouped", nsamples = 300, group = "combIUCN")
dev.off()

