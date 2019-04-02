library(dplyr)
library(magrittr)
library(tidyverse)
library(ggplot2)

# read in parasite richness data
carnRICH <- read.csv("./Data/JPEK/allDat_totrich_carngroup.csv") %>% rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Carnivores")
unggRICH <- read.csv("./Data/JPEK/allDat_totrich_unggroup.csv") %>% rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Ungulates")
primRICH <- read.csv("./Data/JPEK/allDat_totrich_primgroup.csv") %>% rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Primates")

combRICH <- full_join(unggRICH, primRICH) %>% full_join(carnRICH)

# read in parasite transmission mode data
carnTM <- read.csv("./Data/JPEK/allDat_parastrans_carngroup.csv")
unggTM <- read.csv("./Data/JPEK/allDat_parastrans_unggroup.csv")
primTM <- read.csv("./Data/JPEK/allDat_parastrans_primgroup.csv")

# read in parasite type data
carnTYPE <- read.csv("./Data/JPEK/allDat_parastype_carngroup.csv")
unggTYPE <- read.csv("./Data/JPEK/allDat_parastype_unggroup.csv")
primTYPE <- read.csv("./Data/JPEK/allDat_parastype_primgroup.csv")

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
  stat_summary(fun.y=mean, geom="point", shape=23, size=2, fill = 'grey20')
combRICH_violin + 
  scale_fill_brewer(palette="Dark2")
