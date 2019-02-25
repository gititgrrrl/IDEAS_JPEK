library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

GMPD_pairs <- read_csv("GMPD_pairs.csv")
GMPD_threat <- read_csv("GMPD_threat.csv")

par_rich <- GMPD_threat %>%
  group_by(hostName, parasiteName, IUCN.Status.1.2) %>%
  distinct() %>% # to get distinct host parasite pairs
  group_by(hostName) %>%
  tally() # to get parasite richness for each host


parRichness_threat <- full_join(par_rich, GMPD_threat, by = "hostName")

IUCN.levels <- c("NE", "DD", "LC", "NT", "VU", "EN", "CR", "EW")

parRichness_threat$IUCN.levels <- factor(merge$IUCN.Status.1.2, levels = IUCN.levels)


violin_plots <- parRichness_threat %>% filter(IUCN.levels != "DD", !is.na(IUCN.levels)) %>%
ggplot(aes(x = factor(IUCN.levels), y = log(n), color = IUCN.levels)) +
  geom_violin() + facet_grid(Group ~ .) +
  ggtitle("Parasite richness across IUCN levels in GMPD") +
  labs (y = "log(Parasite Richness)", x = "IUCN Classifications") +
  theme(legend.position="none")


table <- parRichness_threat %>%
  select(hostName, Group, IUCN.levels) %>%
    distinct() %>%
  count(Group, IUCN.levels)


grid.table(table)
colnames(table)[which(names(table) == "n")] <- "parasite richness"


print(violin_plots)


