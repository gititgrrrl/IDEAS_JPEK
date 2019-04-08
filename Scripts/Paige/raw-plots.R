library(tidyverse)
library(magrittr)
library(ggplot2) 
library(cowplot)

# richness
carnRICH <- read.csv("./Data/JPEK/allDat_totrich_carngroup.csv") %>% dplyr::rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Carnivores")
unggRICH <- read.csv("./Data/JPEK/allDat_totrich_unggroup.csv") %>% dplyr::rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Ungulates")
primRICH <- read.csv("./Data/JPEK/allDat_totrich_primgroup.csv") %>% dplyr::rename(Threat_status = combIUCN) %>% mutate(Host_Order = "Primates")

# mode 
carnTM <- read.csv("./Data/JPEK/allDat_parastrans_carngroup.csv")%>% mutate(Host_Order = "Carnivores")
unggTM <- read.csv("./Data/JPEK/allDat_parastrans_unggroup.csv")%>% mutate(Host_Order = "Ungulates")
primTM <- read.csv("./Data/JPEK/allDat_parastrans_primgroup.csv") %>% mutate(Host_Order = "Primates")

# type
carnTYPE <- read.csv("./Data/JPEK/allDat_parastype_carngroup.csv")%>% mutate(Host_Order = "Carnivores")
unggTYPE <- read.csv("./Data/JPEK/allDat_parastype_unggroup.csv")%>% mutate(Host_Order = "Ungulates")
primTYPE <- read.csv("./Data/JPEK/allDat_parastype_primgroup.csv") %>% mutate(Host_Order = "Primates")

# fulls
richy <- bind_rows(carnRICH, unggRICH, primRICH) 
trans <- bind_rows(carnTM, primTM, unggTM) %>% mutate(propClose=parRichCloseOnly/parRichTransKnown)
parTypes <- bind_rows(carnTYPE, primTYPE, unggTYPE) %>% mutate(propMicro=parRich_micro/parRich_alltypes)


levels(richy$Threat_status) <- c("NT", "T")
levels(trans$combIUCN) <- c("NT", "T")
levels(parTypes$combIUCN) <- c("NT", "T")

richy %>%
  ggplot(aes(x=Threat_status, y=parRich)) + 
  geom_boxplot()+ 
  ylab("Parasite richness \n (uncorrected)") + xlab("") + 
  facet_grid(.~Host_Order) -> pRich

trans %>% 
  ggplot(aes(x=combIUCN, y=propClose)) + 
  geom_boxplot() + 
  ylab("Proportion close-contact \n (uncorrected)") +xlab("") + 
  facet_grid(.~Host_Order) -> pClose

parTypes %>% 
  ggplot(aes(x=combIUCN, y=propMicro)) + 
  geom_boxplot() + 
  ylab("Proportion microparasite \n (uncorrected)") + xlab("") + 
  facet_grid(.~Host_Order) -> pType

plot_grid(pRich, pClose, pType ,ncol=1)
