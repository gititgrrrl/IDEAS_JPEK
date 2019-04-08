library(dplyr)
library(magrittr)
library(tidyverse)
library(ggplot2)

### READ IN RESIDUALS FROM SIMPLE MODELS <<<<<<<<<<<<<<<<<<<<<<<<<
filenames <- list.files(path = "./FINAL/simple/", pattern = "_resid.RDS$") 
for(i in 1:length(filenames)){
  outname <- unlist(strsplit(x = filenames[i], split = '.RDS'))[1]
  dat <- readRDS(paste0("./FINAL/simple/", filenames[i]))
  assign(x = outname, value = dat)
}
### END NEW CODE <<<<<<<<<<<<<<<<<<<<<<<<<


# read in parasite richness data
carnRICH <- read.csv("./Data/JPEK/allDat_totrich_carngroup.csv") %>% dplyr::rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Carnivores")
unggRICH <- read.csv("./Data/JPEK/allDat_totrich_unggroup.csv") %>% dplyr::rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Ungulates")
primRICH <- read.csv("./Data/JPEK/allDat_totrich_primgroup.csv") %>% dplyr::rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Primates")


### ADD COLUMN OF DATA TO ORIGINAL DATA <<<<<<<<<<<<<<<<<<<<<<<<<
# Here I'm adding a new column called "Residuals" to each of your datasets
carnRICH <- cbind(carnRICH, Residuals = simple_totrich_carngroup_resid[, "Estimate"])
unggRICH <- cbind(unggRICH, Residuals = simple_totrich_unggroup_resid[, "Estimate"])
primRICH <- cbind(primRICH, Residuals = simple_totrich_primgroup_resid[, "Estimate"])
### END NEW CODE <<<<<<<<<<<<<<<<<<<<<<<<<

combRICH <- full_join(unggRICH, primRICH) %>% full_join(carnRICH) # <<<<<<<<< NOW YOU WILL SEE A NEW COLUMN CALLED 'RESIDUALS' HERE, AND THAT IS THE COLUMN YOU SHOULD USE FOR Y-AXIS INSTEAD OF 'TOTRICH' IN PLOTS

# read in parasite transmission mode data
carnTM <- read.csv("./Data/JPEK/allDat_parastrans_carngroup.csv")
unggTM <- read.csv("./Data/JPEK/allDat_parastrans_unggroup.csv")
primTM <- read.csv("./Data/JPEK/allDat_parastrans_primgroup.csv")

### ADD COLUMN OF DATA TO ORIGINAL DATA <<<<<<<<<<<<<<<<<<<<<<<<<
# Here I'm adding a new column called "Residuals" to each of your datasets
carnTM <- cbind(carnTM, Residuals = simple_parastrans_carngroup_resid[, "Estimate"])
unggTM <- cbind(unggTM, Residuals = simple_parastrans_unggroup_resid[, "Estimate"])
primTM <- cbind(primTM, Residuals = simple_parastrans_primgroup_resid[, "Estimate"])
### END NEW CODE <<<<<<<<<<<<<<<<<<<<<<<<<

# read in parasite type data
carnTYPE <- read.csv("./Data/JPEK/allDat_parastype_carngroup.csv")
unggTYPE <- read.csv("./Data/JPEK/allDat_parastype_unggroup.csv")
primTYPE <- read.csv("./Data/JPEK/allDat_parastype_primgroup.csv")

### ADD COLUMN OF DATA TO ORIGINAL DATA <<<<<<<<<<<<<<<<<<<<<<<<<
# Here I'm adding a new column called "Residuals" to each of your datasets
carnTYPE <- cbind(carnTYPE, Residuals = simple_parastype_carngroup_resid[, "Estimate"])
unggTYPE <- cbind(unggTYPE, Residuals = simple_parastype_unggroup_resid[, "Estimate"])
primTYPE <- cbind(primTYPE, Residuals = simple_parastype_primgroup_resid[, "Estimate"])
### END NEW CODE <<<<<<<<<<<<<<<<<<<<<<<<<

# Basic violin plots

carnRICH_violin <- ggplot(carnRICH, aes(x=Threat_status, y=parRich, fill=Threat_status)) +
  geom_violin(trim = F) +
  geom_boxplot(width=0.1) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill = 'grey20')
carnRICH_violin + 
  scale_fill_brewer(palette="Dark2")

unggRICH_violin <- ggplot(unggRICH, aes(x=Threat_status, y=parRich, fill=Threat_status)) +
  geom_violin(trim = F) +
  geom_boxplot(width=0.1) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill = 'grey20')
unggRICH_violin + 
  scale_fill_brewer(palette="Dark2")

primRICH_violin <- ggplot(primRICH, aes(x=Threat_status, y=parRich, fill=Threat_status)) +
  geom_violin(trim = F) +
  geom_boxplot(width=0.1) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill = 'grey20')
primRICH_violin + 
  scale_fill_brewer(palette="Dark2")

combRICH_violin <- ggplot(combRICH, aes(x=Host_Order, y=parRich, fill=Threat_status)) +
  geom_violin(trim = F, position=position_dodge(1)) +
  geom_boxplot(width=0.05, position=position_dodge(1)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, position = position_dodge(1), aes(group = Threat_status), fill = "white") # <<<< THIS IS THE ONLY LINE I CHANGED (I ALSO CHANGED THE FILL COLOR JUST SO IT'S EASIER TO SEE, BUT YOU CAN CHANGE IT BACK IF YOU WANT)
combRICH_violin + 
  scale_fill_brewer(palette="Dark2")
