### QUESTIONS
# > My estimates of parasit richness didn't match Paige's

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

### IMPORT & FORMAT DATA ----
browseURL("http://esapubs.org/archive/ecol/e090/184/metadata.htm") # link to metadata for Pantheria
download.file(url = "http://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR93_Aug2008.txt", destfile = "Data/LH/LH.txt")  
# download mammals life history data from the internet, save to the location specified by 'destfile'

LH <- read_delim("~/Google Drive/GMPD/GlobalParasites/Data/LH/LH.txt", "\t", escape_double = FALSE, na = "-999", trim_ws = TRUE) # Set NA's to -999

LH_sub <- LH[, c(5, 6, 18, 23, 29, 32, 41, 44, 47)] # create a new data frame with a subset of the columns
names(LH_sub) <- c("HostName", "HostActivCycle", "HostHomeRange", "HostMaxLifespan", "HostGroupSize", "HostTrophic", "HostSpeciesRange", "HostMeanLat", "HostMeanLong")
LH_sub$HostName <- gsub(" ", "_", LH_sub$HostName) # replace space with underscore
LH_sub %<>% 
  mutate(HostActivCycle = factor(HostActivCycle, labels = c("nocturnal", "crepuscular", "diurnal")),
         HostTrophic = factor(HostTrophic, labels=c("herbivore", "omnivore", "carnivore")), 
         HostSpeciesRange = HostSpeciesRange/1000,
         HostMaxLifespan = HostMaxLifespan/12)

Dat <- read_csv("Data/JPEK/dat-residual-corrected.csv") # errors on import, but unrelated to the columns we are interested in
TempDat <- Dat %>%
  rename(IUCN = IUCN.Status.1.2,
         HostName = hostName,
         HostGroup = Group,
         ParName = parasiteName,
         ParTrans_Close = close,
         NumHostCitations = numHostCitations,
         UngPrimMeanGroupSize = comGrpSize, 
         CarnGroupSize = carnGrp) %>%
  mutate(HostMass = Mass.g/1000,
         HostIUCNcomb = ifelse(IUCN %in% c("CR", "EN", "VU", "NT"), "threatened", ifelse(IUCN=="LC", "not_threatened", NA))) %>% # had 1 data deficient & 1 EP--convert those to NA
  select(HostName, HostMass, HostIUCNcomb, HostGroup, ParName, ParType, ParTrans_Close, NumHostCitations, UngPrimMeanGroupSize, CarnGroupSize) %>%
  distinct()
TempDat$CarnGroupSize <- tolower(TempDat$CarnGroupSize)

# NOTE: ~10% of host-parasite combinations don't have parasite transmission information
RichnessByTrans <- TempDat %>%
  select(HostName, ParName, ParTrans_Close) %>%
  group_by(HostName) %>%
  summarize(ParRich = sum(!is.na(ParName)),
            ParRich_Close = sum(ParTrans_Close, na.rm = TRUE))
RichnessByParType <- TempDat %>%
  select(HostName, ParType) %>%
  group_by(HostName, ParType) %>%
  summarize(ParRich = n()) %>%
  distinct() %>%
  group_by(HostName) %>%
  spread(ParType, ParRich) %>%
  select(-Prion, -Fungus) # these categories don't have enough data
colnames(RichnessByParType)[-1] <- paste0("ParRich_", colnames(RichnessByParType)[-1])
RichnessByParType[is.na(RichnessByParType)] <- 0

FinalDat <- TempDat %>%
  select(-ParName, -ParType, -ParTrans_Close) %>%
  distinct() %>%
  left_join(RichnessByTrans, by = "HostName") %>%
  left_join(RichnessByParType, by = "HostName") %>%
  left_join(LH_sub, by = "HostName") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(logParRich = log(ParRich),
         logHostSpeciesRange = log(HostSpeciesRange),
         logHostMass = log(HostMass),
         logNumHostCitations = log(NumHostCitations),
         propParCloseTrans = ParRich_Close/ParRich)
FinalDat$UngPrimMeanGroupSize[is.nan(FinalDat$UngPrimMeanGroupSize)] <- NA
saveRDS(FinalDat, "~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen/FinalDat.RDS")

### INITIAL PLOTS OF RAW DATA ----
# Density plots of response variable
FuncDensityHist(dat_df = FinalDat, x_nam = "ParRich", facet_nam = "HostIUCNcomb", col_nam = "HostGroup", cens = FALSE, cens_xint = NULL, plot_labs = list(x = "log(Parasite Richness)", y = NULL, title = NULL))

# Pairwise plots
# > carnivores--no host group size except a few had size =1. 
ggpairs(data = FinalDat,
        title = "GMPD - Pantheria: HostTraits",
        mapping = aes(color = HostGroup, alpha = 0.1),
        lower = list(continuous = "smooth_loess"),
        columns = c("logHostMass", "HostMaxLifespan", "HostCombGroupSize", "logHostSpeciesRange", "HostMeanLat", "logParRich"),
        upper = list(continuous = wrap("cor", size = 3))) +
  theme_bw(base_size = 8)

# Scatterplots showing the relationship between maximum longevity and parasite richness, with a different facet for each group and different color for each activity level (log response)
(p_top1 <- ggplot(data = FinalDat, aes(x = HostMaxLifespan, fill = HostGroup)) +
    geom_density(alpha = 0.4) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Host Lifespan (yrs)", y = NULL) +
    theme_bw() +
    theme(legend.position="top"))

(p_top2 <- ggplot(data = FinalDat, aes(x = HostMaxLifespan, y = ParRich, color = HostGroup)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = HostGroup), show.legend = FALSE) +
  geom_point(alpha = 0.6, show.legend = FALSE) +
  scale_y_continuous(trans='log', breaks= pretty_breaks()) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Host Lifespan (yrs)", y = "Parasite Richness") +
  theme_bw() + theme(legend.position="top"))

# ... species range (log-log)
(p_bottom1 <- ggplot(data = FinalDat, aes(x = HostSpeciesRange, fill = HostGroup)) +
  geom_density(alpha = 0.4, show.legend = FALSE) +
  scale_x_continuous(trans='log', breaks= pretty_breaks()) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Host Species Range (1000km2)", y = NULL) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  
(p_bottom2 <- ggplot(data = FinalDat, aes(x = HostSpeciesRange, y = ParRich, color = HostGroup)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = HostGroup), show.legend = FALSE) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_x_continuous(trans='log', breaks= pretty_breaks()) +
  scale_y_continuous(trans='log', breaks= pretty_breaks()) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Host Species Range (1000km2)", y = "Parasite Richness") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)))

# Combine the plots on a single page, with a shared legend
legend <- get_legend(p_top1)
plots_top <- plot_grid(p_top1 + theme(legend.position="none"), p_top2, ncol = 2, rel_widths = c(0.7, 1))
plots_bottom <- plot_grid(p_bottom1, p_bottom2, ncol = 2, rel_widths = c(0.7, 1))
(plot_final <- plot_grid(legend, plots_top, plots_bottom, nrow = 3, rel_heights = c(0.1, 1, 1)))
saveRDS(plot_final, "plot_final.RDS")

### ANALYSES ----

### Pantheria has social group size vs. population group size--Paige is getting group data (THERE SEEM TO BE PROBLEMS--SEE RESULTS FOR CARNIVORES--analyze only primates & ungulates)
### Do sampling bias by the other counts as well!!!!
### Predictors to try--mean latitude
### Responses to try--parasite broken down by transmission type, parasite broken down by type
### For parasites, just a subset of transmission--like % close only, or % vector
### Response = virus type
### 
### TONIGHT:
### > Group size 2 ways (Sonia categ on carnivores & only ungu/prim)*%close
### > Web of science citations??
### > Check residual plots--account for it?
### > What are error bars in plots?
### > ??Response subset by virus type, skip prions & fungus
### Final plots on best model

FinalDat <- readRDS("~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen/FinalDat.RDS")

ModDat <- FinalDat %>%
  select(HostName, ParRich, propParCloseTrans, NumHostCitations, HostIUCNcomb, HostGroup, logHostSpeciesRange, logHostMass, HostMaxLifespan, HostMeanLat) # 450 records
ModDat <- ModDat[complete.cases(ModDat),] #only 255 records b/c 120 missing species range ans 149 missing max lifespan. With host group size included, only 156 records left.
# Tried these models--didn't work well or didn't make a difference
# > Negative binomial on mod1

# MOD 1 - ZERO TRUNCATED POISSON. FULL MODEL ----
# > What we need is a one-inflated, zero-truncated Poisson!
# > Interaction btwn species range & group is not important.
# > Model is not capturing the high number of singletons, and its decline is too sharp. Fits primates best out of the 3 groups. Fits threatened better than not_threatened. # > Modify>> primates have smaller residual variance. Threatened species also have smaller residual variance, but hold off on this. A lot of the issues with residual plots may be fixed by sigma ~ Group
mod1_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 NumHostCitations +
                 logHostSpeciesRange*HostGroup +
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans*HostIUCNcomb),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 
summary(mod1_brm)
mod1_brm <- add_ic(mod1_brm, ic = "loo", reloo = TRUE) # 11 observations with a pareto_k > 0.7
plot(mod1_brm$loo, label_points = TRUE) # most problematic are #85**, 36, 76, 115--they have k >1
View(ModDat[c(85, 36, 76, 115),])
saveRDS(mod1_brm, "mod1_brm.RDS")

# MOD 2 - LOG-TRANSFORM for CITATIONS ----
# Also removed the interaction of species range with host group
# Great improvement on the model fit
# HostMaxLifespan had 95%CI of 0 to 0.01, so can probably drop it
# Interaction btwn IUCN and %close transmission maybe not important
mod2_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations) + # <<<<<<
                 logHostSpeciesRange + # <<<<<<
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans*HostIUCNcomb),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

mod2_brm <- add_ic(mod2_brm, ic = "loo", reloo = TRUE) # 2 observations with k > 0.7. Worked fine with reloo
plot(mod2_brm$loo, label_points = TRUE) # most problematic is 115 (k > 1)
mod2_brm <- add_ic(mod2_brm, ic = "kfold")
saveRDS(mod2_brm, "mod2_brm.RDS")

# MOD 3 - REPLACE IUCN*%CLOSE INTERACTION WITH HOST GROUP*%CLOSE INTERACTION ----
mod3_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations) + 
                 logHostSpeciesRange + 
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans*HostGroup), # <<<<<<
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

mod3_brm <- add_ic(mod3_brm, ic = "loo", reloo = TRUE) # 3 observations with k > 0.7. Worked fine with reloo
mod3_brm <- add_ic(mod3_brm, ic = "kfold")
plot(mod3_brm$loo, label_points = TRUE) # most problematic is 115 (k > 1)
saveRDS(mod3_brm, "mod3_brm.RDS")

# **** BEST MODEL: MOD4 - ADD INTERACTION BTWN CITATIONS AND IUCN STATUS ----
# PP-check shows still slight problems with ungulates and a bit of problem with carnivores & non-threatened, but not bad
# Loo-pit shows S-shape
mod4_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations)*HostIUCNcomb + # <<<<<<
                 logHostSpeciesRange + 
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans*HostGroup), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

mod4_brm <- add_ic(mod4_brm, ic = "loo", reloo = TRUE) # 5 observations with k > 0.7. Worked fine with reloo
mod4_brm <- add_ic(mod4_brm, ic = "kfold")
plot(mod4_brm$loo, label_points = TRUE)
saveRDS(mod4_brm, "mod4_brm.RDS")

# Run best model (MOD 4) with glmmtmb to double-check results
mod4_tmb <- glmmTMB(ParRich ~
                      HostIUCNcomb*HostGroup + 
                      log(NumHostCitations)*HostIUCNcomb +
                      logHostSpeciesRange + 
                      logHostMass*HostGroup +
                      HostMaxLifespan +
                      propParCloseTrans*HostGroup,
                    data = ModDat,
                    family=truncated_poisson)
summary(mod4_tmb) # results are similar to brms!

mod_pred = predict(mod4_tmb, se.fit=FALSE, type = "response")

# MOD5 - ADD INTERACTION BTWN CITATIONS AND HOST GROUP ----
mod5_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations)*HostGroup + # <<<<<<
                 logHostSpeciesRange + 
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans*HostGroup), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

mod5_brm <- add_ic(mod5_brm, ic = "loo", reloo = TRUE) # 4 observations with k > 0.7. Worked fine with reloo
mod5_brm <- add_ic(mod5_brm, ic = "kfold")
plot(mod5_brm$loo, label_points = TRUE) # most problematic are 131 and 115 (k > 1)
saveRDS(mod5_brm, "mod5_brm.RDS")

# MOD6 - ADD INTERACTION BTWN LIFESPAN AND HOST GROUP (BASED ON MOD4) ----
mod6_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations)*HostIUCNcomb + 
                 logHostSpeciesRange + 
                 logHostMass*HostGroup +
                 HostMaxLifespan*HostGroup + # <<<<<<
                 propParCloseTrans*HostGroup), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10))  

mod6_brm <- add_ic(mod6_brm, ic = "loo", reloo = TRUE) # 4 observations with k > 0.7. Worked fine with reloo
mod6_brm <- add_ic(mod6_brm, ic = "kfold")
plot(mod6_brm$loo, label_points = TRUE) # most problematic are 115 and 131 (k > 1)
saveRDS(mod6_brm, "mod6_brm.RDS")

# MOD7 - ADD MEAN LATITUDE*GROUP (BASED ON MOD4) ----
mod7_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations)*HostIUCNcomb + 
                 logHostSpeciesRange + 
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans*HostGroup +
                 HostMeanLat*HostGroup), # <<<<<<
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

mod7_brm <- add_ic(mod7_brm, ic = "loo", reloo = TRUE) # 4 observations with k > 0.7. Worked fine with reloo
mod7_brm <- add_ic(mod7_brm, ic = "kfold")
plot(mod7_brm$loo, label_points = TRUE)
saveRDS(mod7_brm, "mod7_brm.RDS")

# MOD8 - REMOVE %CLOSE VIRUS (BASED ON MOD7) ----
# NEED TO THINK ABOUT WHAT THIS IS TELLING US... <<<<<<<<
mod8_brm <- brm(
  data = ModDat, 
  family = poisson,
  formula = bf(ParRich | trunc(lb = 1) ~
                 HostIUCNcomb*HostGroup + 
                 log(NumHostCitations)*HostIUCNcomb +
                 logHostSpeciesRange + 
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 HostMeanLat*HostGroup), 
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 

mod8_brm <- add_ic(mod8_brm, ic = "loo", reloo = TRUE) # 4 observations with k > 0.7. Worked fine with reloo
mod8_brm <- add_ic(mod8_brm, ic = "kfold")
plot(mod8_brm$loo, label_points = TRUE) # problematic are #85, 115, 131 (k >1)
saveRDS(mod8_brm, "mod8_brm.RDS")





### STILL TO DO ----
### SEPARATELY EXAMINE POTENTIAL IMPORTANCE OF GROUPSIZE
# Note that I can't actually compare it to the other models because I'm not using the exact same data records
# > for carnivores...
ModDatCarn <- FinalDat %>%
  select(HostName, ParRich, propParCloseTrans, NumHostCitations, HostIUCNcomb, HostGroup, logHostSpeciesRange, logHostMass, HostMaxLifespan) %>%
  filter(HostGroup == "carnivore") %>%
  left_join(TempDat[,c("HostName", "CarnGroupSize")], by = "HostName")
ModDatCarn <- ModDatCarn[complete.cases(ModDatCarn),]


# > for ungulates and primates...
ModDatUngPrim <- FinalDat %>%
  select(HostName, ParRich, propParCloseTrans, NumHostCitations, HostIUCNcomb, HostGroup, logHostSpeciesRange, logHostMass, HostMaxLifespan) %>%
  filter(HostGroup %in% c("ungulate", "primate")) %>%
  left_join(TempDat[,c("HostName", "UngPrimMeanGroupSize")], by = "HostName")
ModDatUngPrim <- ModDatUngPrim[complete.cases(ModDatUngPrim),]

