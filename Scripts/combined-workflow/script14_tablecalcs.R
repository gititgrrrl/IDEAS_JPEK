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

### Lists of best models only----
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
<<<<<<< HEAD
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_ParasType.RDS"),
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_ParasType.RDS"),
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_ParasType.RDS"))

### Relevel threat status ----
# Parasite richness
#...carnivores
allDat <- read_csv("Data/JPEK/script4.csv")

fullDat_carngroup <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() # 148 records
fullDat_carngroup <- fullDat_carngroup[complete.cases(fullDat_carngroup),] # only 74 records

fullDat_carngroup$groupSizeCar <- factor(fullDat_carngroup$groupSizeCar, levels = c("non_group", "group"))
fullDat_carngroup$combIUCN <- factor(fullDat_carngroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_carngroup_c <- fullDat_carngroup %>%
  mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

REV_allIntBrm_carngroup_c <- brm(
  data = fullDat_carngroup_c,
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*groupSizeCar +
                 combIUCN*logHostMass_c +
                 combIUCN*hostMaxLifespan_c +
                 combIUCN*absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
summary(REV_allIntBrm_carngroup_c)
saveRDS(REV_allIntBrm_carngroup_c, "./Data/JPEK/allInt/REV_allInt_brm_carngroup_c_all.RDS")

#...primates
fullDat_primgroup <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() # 142 records
fullDat_primgroup <- fullDat_primgroup[complete.cases(fullDat_primgroup),] # only 73 records

fullDat_primgroup$combIUCN <- factor(fullDat_primgroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_primgroup_c <- fullDat_primgroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

REV_fullBrm_primgroup_c <- brm(
  data = fullDat_primgroup_c,
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*logGroupSizePriUng_c +
                 logHostMass_c +
                 hostMaxLifespan_c +
                 absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
summary(REV_fullBrm_primgroup_c)
saveRDS(REV_fullBrm_primgroup_c, "./Data/JPEK/full/REV_full_brm_primgroup_c_all.RDS") # model

#...ungulates
fullDat_unggroup <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() # 102 records
fullDat_unggroup <- fullDat_unggroup[complete.cases(fullDat_unggroup),] # only 60 records

fullDat_unggroup$combIUCN <- factor(fullDat_unggroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_unggroup_c <- fullDat_unggroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

# ...full model (some threat interactions omitted)
REV_fullBrm_unggroup_c <- brm(
  data = fullDat_unggroup_c,
  family = poisson,
  formula = bf(parRich | trunc(lb = 1) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*logGroupSizePriUng_c +
                 logHostMass_c +
                 hostMaxLifespan_c +
                 absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
summary(REV_fullBrm_unggroup_c)
saveRDS(REV_fullBrm_unggroup_c, "./Data/JPEK/full/REV_full_brm_unggroup_c_all.RDS")

# Parasite transmission
#...carnivores
fullDat_parastrans_carngroup <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizeCar, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

fullDat_parastrans_carngroup <- fullDat_parastrans_carngroup[complete.cases(fullDat_parastrans_carngroup),] # only 71 records

fullDat_parastrans_carngroup$groupSizeCar <- factor(fullDat_parastrans_carngroup$groupSizeCar, levels = c("non_group", "group"))
fullDat_parastrans_carngroup$combIUCN <- factor(fullDat_parastrans_carngroup$combIUCN, levels = c("threatened", "not_threatened"))

fullDat_parastrans_carngroup_c <- fullDat_parastrans_carngroup %>%
  mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

# ...full model (some threat interactions omitted)
REV_fullBrm_carngroup_c_parastrans <- brm(
  data = fullDat_parastrans_carngroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*groupSizeCar +
                 logHostMass_c +
                 hostMaxLifespan_c +
                 absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
summary(REV_fullBrm_carngroup_c_parastrans)
saveRDS(REV_fullBrm_carngroup_c_parastrans, "./Data/JPEK/full/REV_full_brm_carngroup_c_parastrans.RDS") # model

#...primates
fullDat_parastrans_primgroup <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 130

fullDat_parastrans_primgroup <- fullDat_parastrans_primgroup[complete.cases(fullDat_parastrans_primgroup),] # only 69 records

fullDat_parastrans_primgroup$combIUCN <- factor(fullDat_parastrans_primgroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_parastrans_primgroup_c <- fullDat_parastrans_primgroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))
# ...model with ALL threat interactions
REV_allIntBrm_primgroup_c_parastrans <- brm(
  data = fullDat_parastrans_primgroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*logGroupSizePriUng_c +
                 combIUCN*logHostMass_c +
                 combIUCN*hostMaxLifespan_c +
                 combIUCN*absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(REV_allIntBrm_primgroup_c_parastrans)
saveRDS(REV_allIntBrm_primgroup_c_parastrans, "./Data/JPEK/allInt/REV_allInt_brm_primgroup_c_parastrans.RDS")

#...ungulates
fullDat_parastrans_unggroup <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 93

fullDat_parastrans_unggroup <- fullDat_parastrans_unggroup[complete.cases(fullDat_parastrans_unggroup),] # only 58 records

fullDat_parastrans_unggroup$combIUCN <- factor(fullDat_parastrans_unggroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_parastrans_unggroup_c <- fullDat_parastrans_unggroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

REV_simpBrm_unggroup_c_parastrans <- brm(
  data = fullDat_parastrans_unggroup_c,
  family = binomial,
  formula = bf(parRichCloseOnly|trials(parRichTransKnown) ~
                 combIUCN*logNumHostCitations_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged
summary(REV_simpBrm_unggroup_c_parastrans)
saveRDS(REV_simpBrm_unggroup_c_parastrans, "./Data/JPEK/simple/REV_simp_brm_unggroup_c_parastrans.RDS")

# Parasite type
#...carnivores
fullDat_parastype_carngroup <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, logNumHostCitations, combIUCN, groupSizeCar, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 148

fullDat_parastype_carngroup <- fullDat_parastype_carngroup[complete.cases(fullDat_parastype_carngroup),] # only 74 records

fullDat_parastype_carngroup$groupSizeCar <- factor(fullDat_parastype_carngroup$groupSizeCar, levels = c("non_group", "group"))
fullDat_parastype_carngroup$combIUCN <- factor(fullDat_parastype_carngroup$combIUCN, levels = c("threatened", "not_threatened"))

fullDat_parastype_carngroup_c <- fullDat_parastype_carngroup %>%
  mutate(logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

# ...model with ALL threat interactions
REV_allIntBrm_carngroup_c_parastype <- brm(
  data = fullDat_parastype_carngroup_c,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*groupSizeCar +
                 combIUCN*logHostMass_c +
                 combIUCN*hostMaxLifespan_c +
                 combIUCN*absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# # quick checks
summary(REV_allIntBrm_carngroup_c_parastype)
saveRDS(REV_allIntBrm_carngroup_c_parastype, "./Data/JPEK/allInt/REV_allInt_brm_carngroup_c_parastype.RDS")

#...primates
fullDat_parastype_primgroup <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, logNumHostCitations, combIUCN, groupSizePriUng, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 142 records

fullDat_parastype_primgroup <- fullDat_parastype_primgroup[complete.cases(fullDat_parastype_primgroup),] # only 73 records

fullDat_parastype_primgroup$combIUCN <- factor(fullDat_parastype_primgroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_parastype_primgroup_c <- fullDat_parastype_primgroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

REV_simpBrm_primgroup_c_parastype <- brm(
  data = fullDat_parastype_primgroup_c,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(REV_simpBrm_primgroup_c_parastype)
saveRDS(REV_simpBrm_primgroup_c_parastype, "./Data/JPEK/simple/REV_simp_brm_primgroup_c_parastype.RDS")

#...ungulates
fullDat_parastype_unggroup <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, logNumHostCitations, combIUCN, groupSizePriUng, bacteriaRich, virusRich, protozoaRich, fungusRich, prionRich, helminthRich, arthropodRich, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  mutate(parRich_micro=(bacteriaRich + virusRich + protozoaRich + fungusRich + prionRich),
         parRich_alltypes =(parRich_micro + helminthRich + arthropodRich)) # 102 records

fullDat_parastype_unggroup <- fullDat_parastype_unggroup[complete.cases(fullDat_parastype_unggroup),] # only 60 records

fullDat_parastype_unggroup$combIUCN <- factor(fullDat_parastype_unggroup$combIUCN, levels = c("threatened", "not_threatened"))
fullDat_parastype_unggroup_c <- fullDat_parastype_unggroup %>%
  mutate(logGroupSizePriUng_c = logGroupSizePriUng - mean(logGroupSizePriUng, na.rm = TRUE),
         logNumHostCitations_c = logNumHostCitations - mean(logNumHostCitations, na.rm = TRUE),
         logHostSpeciesRange_c = logHostSpeciesRange - mean(logHostSpeciesRange, na.rm = TRUE),
         logHostMass_c = logHostMass - mean(logHostMass, na.rm = TRUE),
         hostMaxLifespan_c = hostMaxLifespan - mean(hostMaxLifespan, na.rm = TRUE),
         absHostMeanLat_c = absHostMeanLat - mean(absHostMeanLat, na.rm = TRUE))

# ...full model (some threat interactions omitted)
REV_fullBrm_unggroup_c_parastype <- brm(
  data = fullDat_parastype_unggroup_c,
  family = binomial,
  formula = bf(parRich_micro|trials(parRich_alltypes) ~
                 combIUCN*logNumHostCitations_c +
                 combIUCN*logHostSpeciesRange_c +
                 combIUCN*logGroupSizePriUng_c +
                 logHostMass_c +
                 hostMaxLifespan_c +
                 absHostMeanLat_c),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  # may wnat to increase max_tree depth although all converged

# quick checks
summary(REV_fullBrm_unggroup_c_parastype)
saveRDS(REV_fullBrm_unggroup_c_parastype, "./Data/JPEK/full/REV_full_brm_unggroup_c_parastype.RDS") # model

###Read in Excel and convert to RDS
richness_summary <- read_csv("Data/JPEK/richness_summary.csv")
saveRDS(richness_summary, "./Results/FINAL-FIGS/richness_summary.RDS")
paratrans_summary <- read_csv("Data/JPEK/paratrans_summary.csv")
saveRDS(paratrans_summary, "./Results/FINAL-FIGS/paratrans_summary.RDS")
paratype_summary <- read_csv("Data/JPEK/paratype_summary.csv")
saveRDS(paratype_summary, "./Results/FINAL-FIGS/paratype_summary.RDS")
=======
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_parastype.RDS"),
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_parastype.RDS"),
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_parastype.RDS"))
>>>>>>> 06b5f4b67a41254eecc08b4096d53479d28e3bbc
