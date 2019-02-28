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
names(LH_sub) <- c("Species", "HostActivCycle", "HostHomeRange", "HostMaxLifespan", "HostGroupSize", "HostTrophic", "HostSpeciesRange", "HostMeanLat", "HostMeanLong")
LH_sub$Species <- gsub(" ", "_", LH_sub$Species) # replace space with underscore
LH_sub %<>% 
  mutate(HostActivCycle = factor(HostActivCycle, labels = c("nocturnal", "crepuscular", "diurnal")),
         HostTrophic = factor(HostTrophic, labels=c("herbivore", "omnivore", "carnivore")), 
         HostSpeciesRange = HostSpeciesRange/1000,
         HostMaxLifespan = HostMaxLifespan/12)
GMPD_threat <- read_csv("~/Google Drive/GMPD/GlobalParasites/Data/JPEK/GMPD_threat.csv") # get the GMPD data with threats

# Check summary of counts by group
# > ungulates have the fewest host species and most parasite species in GMPD
(Summary_counts <- GMPD_threat %>%
    group_by(Group) %>%
    summarize(TotRecords = n(), 
              UniqueHosts = length(unique(hostName)),
              UniqueParas = length(unique(parasiteName))))

(Summary_PR <- GMPD_threat %>%
    mutate(HostMass = Mass.g/1000) %>%
    select(-Mass.g) %>%
    rename(IUCN = IUCN.Status.1.2,
           HostName = hostName,
           HostGroup = Group,
           ParName = parasiteName,
           ParTrans_Close = close,
           ParTrans_Nonclose = nonclose,
           ParTrans_Vector = vector,
           ParTrans_Intermed = intermediate) %>%
    mutate(IUCNcomb = ifelse(IUCN %in% c("CR", "EN", "VU", "NT"), "Threatened", ifelse(IUCN=="LC", "NotThreatened", NA))) %>% # had 1 data deficient & 1 EP--convert those to NA
    left_join(LH_sub, by = c("HostName" = "Species")) %>%
    select(-`Order.1.2`, -`Family.1.2`, -`Genus.1.2`, -ParPhylum, -Prevalence, -ParasiteTraitsCitation, -SamplingType, -HostEnvironment, -Citation, -IUCN))

# NOTE: ~10% of host-parasite combinations don't have parasite transmission information
RichnessByTrans <- Summary_PR %>%
  select(HostName, ParName, ParTrans_Close:ParTrans_Intermed) %>%
  group_by(HostName, ParName) %>%
  distinct() %>%
  group_by(HostName) %>%
  summarize(ParRich = sum(!is.na(ParName)),
            ParRich_Close = sum(ParTrans_Close, na.rm = TRUE),
            ParRich_Nonclose = sum(ParTrans_Nonclose, na.rm = TRUE),
            ParRich_Vector = sum(ParTrans_Vector, na.rm = TRUE),
            ParRich_Intermed = sum(ParTrans_Intermed, na.rm = TRUE))
MM_ParType <- model.matrix(~ParType - 1, data = Summary_PR)
RichnessByParType <- cbind(Summary_PR[c("HostName", "ParName")], MM_ParType) %>%
  group_by(HostName, ParName) %>%
  distinct() %>%
  group_by(HostName) %>%
  summarize(ParRich_Arthropod = sum(ParTypeArthropod, na.rm = TRUE),
            ParRich_Bacteria = sum(ParTypeBacteria, na.rm = TRUE),
            ParRich_Fungus = sum(ParTypeFungus, na.rm = TRUE),
            ParRich_Helminth = sum(ParTypeHelminth, na.rm = TRUE),
            ParRich_Prion = sum(ParTypePrion, na.rm = TRUE),
            ParRich_Protozoa = sum(ParTypeProtozoa, na.rm = TRUE),
            ParRich_Virus = sum(ParTypeVirus, na.rm = TRUE))
SumHostsSampled <- Summary_PR %>%
  group_by(HostName) %>%
  summarize(NumHostsSampled = sum(HostsSampled))
Summary2_PR <- Summary_PR[, -grep("Par", colnames(Summary_PR))] %>%
  select(-HostsSampled) %>%
  distinct() %>%
  left_join(SumHostsSampled, by = "HostName") %>%
  left_join(RichnessByTrans, by = "HostName") %>%
  left_join(RichnessByParType, by = "HostName")

# Get Paige's data
FinalDat <- read_csv("Data/JPEK/dat-residual-corrected.csv") %>% # errors on import, but unrelated to the columns we are interested in
  rename(HostName = hostName, NumHostCitations = numHostCitations, HostCombGroupSize = comGrpSize) %>%
  select(HostName, NumHostCitations, HostCombGroupSize) %>%
  distinct() %>%
  right_join(Summary2_PR, by = "HostName") %>%
  mutate(logParRich = log(ParRich),
         logHostSpeciesRange = log(HostSpeciesRange),
         logHostMass = log(HostMass),
         logNumHostCitations = log(NumHostCitations),
         propParCloseTrans = ParRich_Close/ParRich)

saveRDS(FinalDat, "~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen/FinalDat.RDS")

### INITIAL PLOTS OF RAW DATA ----
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

### Pantheria has social group size vs. population group size--Paige is getting group data (THERE SEEM TO BE PROBLEMS--SEE RESULTS FOR CARNIVORES)
### Do sampling bias by the other counts as well!!!!
### Predictors to try--mean latitude
### Responses to try--parasite broken down by transmission type, parasite broken down by type
### For parasites, just a subset of transmission--like % close only, or % vector

FinalDat <- readRDS("~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen/FinalDat.RDS")

SubDat <- FinalDat %>%
  select(HostName, ParRich, propParCloseTrans, NumHostCitations, IUCNcomb, HostGroup, logHostSpeciesRange, logHostMass, HostMaxLifespan) # 450 records
SubDat <- SubDat[complete.cases(SubDat),] #only 255 records b/c 120 missing species range ans 149 missing max lifespan. With host group size included, only 156 records left.

# MOD 1
# Start with a Poisson regression (since it's counts) with full predictors. Then from model fit diagnostics, check:
# > if zero-inflation and/or overdispersion needs to be modeled (I tried neg. binomial, it was bad!)
# > if need to log-transform some predictors, e.g., NumHostCitations
# These ALMOST DON'T MATTER--species range (except for carnivores), %parasite_close
# >>>>>>>>>>> CHANGE TO POISSON
# mod1_tmb <- glmmTMB(ParRich ~ 
#                       IUCNcomb*HostGroup +
#                       NumHostCitations +
#                       SpeciesRange +
#                       HostMass +
#                       HostMaxLifespan +
#                       propParCloseTrans,
#                     data = SubDat)
# summary(mod1_tmb)
# saveRDS(mod1_tmb, "mod1_tmb.RDS")

mod1_brm <- brm(
  data = SubDat, 
  family = poisson,
  formula = bf(ParRich ~ 
                 IUCNcomb*HostGroup + 
                 NumHostCitations +
                 logHostSpeciesRange*HostGroup +
                 logHostMass*HostGroup +
                 HostMaxLifespan +
                 propParCloseTrans),
  iter =4000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = .8, max_treedepth = 10)) 
summary(mod1_brm)
mod1_brm <- add_ic(mod1_brm, ic = c("loo")) # 11 observations with a pareto_k > 0.7
plot(mod1_brm$loo, label_points = TRUE) # most problematic are 85, 127, 255**
SubDat[c(85, 127, 255), "HostName"] # Diceros_bicornis, Lutra_lutra, Vulpes_vulpes 
saveRDS(mod1_brm, "mod1_brm.RDS")
# Model diagnostics:
# > Model doesn't fit ungulate or carnivore data well--doesn't capture some very low parasite richness, and is also off for non-threatened

# MOD 2
# Try the poisson model without the worst outliers
SubDat2 <- SubDat[-c(85, 127, 255),]
mod2_brm <- update(mod1_brm,
                   newdata = SubDat2)
summary(mod2_brm)
mod2_brm <- add_ic(mod2_brm, ic = c("loo"), reloo = TRUE) 
plot(mod2_brm$loo, label_points = TRUE) # most problematic are 64, 170, 243 (k >1)
saveRDS(mod2_brm, "mod2_brm.RDS")


### CHECK THE MODEL ----
mod <- mod1_brm
plot(mod)

FuncPlotCoef(mod_list = list(mod1 = mod1_brm), terms_vec = c("c_WY_Yr", "c_Precip", "Ic_PrecipE2"), exclude_terms = FALSE) 

FuncPlotCoef(mod_list = list(mod1 = mod1_brm), terms_vec = c("c_WY_Yr", "c_Precip", "Ic_PrecipE2"), exclude_terms = FALSE) # Separately plot coefficients that require finer scale of resolution

pp_check(mod, type = "dens_overlay", nsamples = 300)
pp_check(mod, type = "ecdf_overlay", nsamples = 300)
pp_check(mod, type = "violin_grouped", nsamples = 300, group = "HostGroup")
pp_check(mod, type = "violin_grouped", nsamples = 300, group = "IUCNcomb")
pp_check(mod, type = "loo_pit_qq", nsamples = 300) 

marginal_effects(mod)
plot(marginal_effects(mod, effects = "logSpeciesRange:HostGroup"), points=TRUE) # marginal effects of logHostMass, with different colors by HostGroup

# Posterior predictive fit
FuncPredVObs(dat_df = SubDat, mod_list = list(mod1 = mod), resp_colnam = "ParRich", obsID_colnam = "HostName", color_colnam = "IUCNcomb", coloraxt_title = "IUCN status", facet_colnam = "HostGroup")

# LOO
