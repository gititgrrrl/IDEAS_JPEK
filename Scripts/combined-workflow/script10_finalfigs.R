### FIGURES FOR DRAFT PAPER

rm(list=ls())

### LOAD PACKAGES ----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "plyr", "brms")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
}

### CREATE BEST MODELS DATA FRAME ----
# Data frame summarizing best model for each host-response combination
temp <- expand.grid(
  hostGroup = c("carnivores", "primates", "ungulates"),
  response = c("rich", "parastrans", "parastype"))
bestModSummary <- data.frame(
  temp,
  bestMod = c("complete", "reduced", "reduced", 
               "reduced", "complete", "simple",
               "complete", "simple", "reduced"), # the model with lowest AIC
  bestComplex = c("complete", "reduced", "reduced", 
                   "reduced", "complete", "reduced",
                   "complete", "complete", "reduced")) # between the two complex models, which one has lower AIC
  
  
### FIG.1 ----
# Effects of threat status on response (richness, % direct transmission, % micro-parasites) for simple vs. full models, separately for each host group

margRich <- margParasTrans <- margParasType <- list()

# Parasite richness
margRich$simple <- list(  # >>>>>> WAITING FOR CENTERED MARGINALS
  carnivores = readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")])
margRich_simple <- rbind.fill(margRich$simple)
margRich_simple$hostGroup <- rep(names(margRich$simple), each = 2)
margRich_simple$modType <- "simple"

margRich$complex <- list( 
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "groupSizeCar", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/full/full_brm_carngroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/full/full_brm_carngroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")])
margRich_complex <- rbind.fill(margRich$complex)
margRich_complex$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
margRich_complex$modType <- "complex"

margRich_df <- rbind.fill(margRich_simple, margRich_complex)
levels(margRich_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")

# Parasite transmission
margParasTrans$simple <- list(
  carnivores = readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])

margParasTrans_simple <- rbind.fill(margParasTrans$simple)
margParasTrans_simple$hostGroup <- rep(names(margParasTrans$simple), each = 2)
margParasTrans_simple$modType <- "simple"

margParasTrans$complex <- list( 
  carnivores = readRDS("./Data/JPEK/full/full_brm_carngroup_c_parastrans_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "parRichTransKnown", "groupSizeCar", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])

margParasTrans_complex <- rbind.fill(margParasTrans$complex)
margParasTrans_complex$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
margParasTrans_complex$modType <- "complex"

margParasTrans_df <- rbind.fill(margParasTrans_simple, margParasTrans_complex) %>%
  mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
         `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
         `upper__` = round(100*(`upper__`/parRichTransKnown), 1))
levels(margParasTrans_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
    
# Parasite type
margParasType$simple <- list(
  carnivores = readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])

margParasType_simple <- rbind.fill(margParasType$simple)
margParasType_simple$hostGroup <- rep(names(margParasType$simple), each = 2)
margParasType_simple$modType <- "simple"

margParasType$complex <- list( 
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_ParasType_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "parRich_alltypes", "groupSizeCar", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])

margParasType_complex <- rbind.fill(margParasType$complex)
margParasType_complex$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
margParasType_complex$modType <- "complex"

margParasType_df <- rbind.fill(margParasType_simple, margParasType_complex) %>%
  mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
         `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
         `upper__` = round(100*(`upper__`/parRich_alltypes), 1))
levels(margParasType_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")

# Create figure
FuncPlot1 <- function(dat, y_nam, remove_xlab, remove_faclab) {
  p <- ggplot(data = dat) +
    geom_errorbar(data = subset(dat, modType=="simple"), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0, color = "blue", position = position_nudge(x = -0.15)) +
    geom_point(data = subset(dat, modType=="simple"), aes_string(x = "combIUCN", y = "estimate__"), shape = 21, color = "blue", fill = "white", size = 3, position = position_nudge(x = -0.15)) +
    geom_errorbar(data = subset(dat, modType=="complex" & groupSizeCar %in% c(NA, "non_group")), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0, color = "black", position = position_nudge(x = 0.15)) +
    geom_point(data = subset(dat, modType=="complex" & groupSizeCar %in% c(NA, "non_group")), aes_string(x = "combIUCN", y = "estimate__"), shape = 21, color = "black", fill = "white", size = 3, position = position_nudge(x = 0.15)) +
    scale_y_continuous(breaks= pretty_breaks()) +
    labs(y = y_nam) +
    theme_bw(base_size = 10) +
    facet_grid(. ~ hostGroup)
  if(remove_faclab) {
    p <- p + 
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()
            )
  }
  if(remove_xlab) {
    p <- p +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  } else {
    p <- p +
      labs(x = "Threat Status") 
  }
}

margParasTrans_fig <- FuncPlot1(dat = margParasTrans_df, y_nam = "% Direct-Obligate", remove_faclab = TRUE, remove_xlab = TRUE)
margParasType_fig <- FuncPlot1(dat = margParasType_df, y_nam = "% Micro-Parasite", remove_faclab = TRUE, remove_xlab = FALSE)
fig1_final <- plot_grid(margParasTrans_fig, margParasType_fig, ncol = 1, rel_heights = c(0.85, 1))

pdf("./Results/FINAL-FIGS/fig1_final.pdf")
fig1_final
dev.off()
  
### FIG.2 ----
# Marginal effects plots for IUCN:group size interaction FOR RICHNESS
# <<<<<< WAIT FOR CENTERED RESULTS
allInt_carngroup_c_me_groupsize <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_me_groupsize.RDS")
carn_groupMarg_plot <- plot(allInt_carngroup_c_me_groupsize, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(xxx, xxx) + # <<<<<<<<<<<<<
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = "Group Category", y = "Parasite Richness", subtitle = "carnivores") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") 

full_primgroup_c_me_groupsize <- readRDS("./Data/JPEK/full/full_brm_primgroup_c_me_groupsize.RDS")
prim_groupMarg_plot <- plot(full_primgroup_c_me_groupsize, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(xxx, xxx) + # <<<<<<<<<<<<<
  labs(x = "log(Group Size)", subtitle = "primates") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") 
# +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())

fullBrm_unggroup_c_me_groupsize <- readRDS("./Data/JPEK/full/full_brm_unggroup_c_me_groupsize.RDS")
ung_groupMarg_plot <- plot(fullBrm_unggroup_c_me_groupsize, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(xxx, xxx) + # <<<<<<<<<<<<<
  labs(x = "log(Group Size)", subtitle = "ungulates") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") 
# +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
#         
fig2_final <- plot_grid(carn_groupMarg_plot, prim_groupMarg_plot, ung_groupMarg_plot, nrow = 1, rel_widths = c(1, 0.85, 0.85))

pdf("./Results/FINAL-FIGS/fig2_final.pdf")
fig2_final
dev.off()

### FIG.3 ----
# Marginal effects plots for IUCN:species range interaction FOR RICHNESS
# <<<<<< WAIT FOR CENTERED RESULTS
allInt_carngroup_c_me_speciesrange <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_me_speciesrange.RDS")
carn_rangeMarg_plot <- plot(allInt_carngroup_c_me_speciesrange, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(xxx, xxx) + # <<<<<<<<<<<<<
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = "log(Species Range)", y = "Parasite Richness", subtitle = "carnivores") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") 

full_primgroup_c_me_speciesrange <- readRDS("./Data/JPEK/full/full_brm_primgroup_c_me_speciesrange.RDS")
prim_rangeMarg_plot <- plot(full_primgroup_c_me_speciesrange, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(xxx, xxx) + # <<<<<<<<<<<<<
  labs(x = "log(Species Range)", subtitle = "primates") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") 
# +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())

fullBrm_unggroup_c_me_speciesrange <- readRDS("./Data/JPEK/full/full_brm_unggroup_c_me_speciesrange.RDS")
ung_rangeMarg_plot <- plot(fullBrm_unggroup_c_me_speciesrange, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(xxx, xxx) + # <<<<<<<<<<<<<
  labs(x = "log(Species Range)", subtitle = "ungulates") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") 
# +
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
#         
fig3_final <- plot_grid(carn_rangeMarg_plot, prim_rangeMarg_plot, ung_rangeMarg_plot, nrow = 1, rel_widths = c(1, 0.85, 0.85))

pdf("./Results/FINAL-FIGS/fig3_final.pdf")
fig3_final
dev.off()