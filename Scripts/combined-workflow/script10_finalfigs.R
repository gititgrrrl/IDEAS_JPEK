### FIGURES FOR DRAFT PAPER

rm(list=ls())

### LOAD PACKAGES ----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "plyr", "brms")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})

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

### FIG.1 REDO (best models only) ----
# Effects of threat status on response (richness, % direct transmission, % micro-parasites) for the best-supported model, separately for each host group
# Parasite richness <<<<<<<<<<<<<<<<<<< NEED TO FIX FOR FINAL, TO GET CENTERED RESULTS!
margRichBest_list <- list(
  carnivores = readRDS("./Data/JPEK/full/full_brm_carngroup_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "groupSizeCar", "estimate__", "lower__", "upper__")],# <<<<<<< CHANGE TO "./Data/JPEK/allInt/allInt_brm_carngroup_me.RDS"
  primates = readRDS("./Data/JPEK/full/full_brm_primgroup_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")], # <<<<<<<<<< CHANGE TO "./Data/JPEK/full/full_brm_primgroup_c_me.RDS"
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")]) # <<<<<< CHANGE TO "./Data/JPEK/full/full_brm_unggroup_c_me.RDS"
margRichBest_df <- rbind.fill(margRichBest_list)
margRichBest_df$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
levels(margRichBest_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")

# Parasite transmission
margParasTransBest_list <- list( 
  carnivores = readRDS("./Data/JPEK/full/full_brm_carngroup_c_parastrans_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "parRichTransKnown", "groupSizeCar", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])
margParasTransBest_df <- rbind.fill(margParasTransBest_list)
margParasTransBest_df$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
levels(margParasTransBest_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")

margParasTransBest_df %<>%
  mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
         `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
         `upper__` = round(100*(`upper__`/parRichTransKnown), 1))

# Parasite type
margParasTypeBest_list <- list(
  carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_ParasType_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "parRich_alltypes", "groupSizeCar", "estimate__", "lower__", "upper__")],
  primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])
margParasTypeBest_df <- rbind.fill(margParasTypeBest_list)
margParasTypeBest_df$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
levels(margParasTypeBest_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")

margParasTypeBest_df %<>%
  mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
         `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
         `upper__` = round(100*(`upper__`/parRich_alltypes), 1))

# Create figure
# Need to create a different plot for carnivores (with group and non-group), so can't just facet by hostGroup
FuncPlot1 <- function(dat, y_nam, remove_xlab, add_title) {
  maxY <- max(dat$`upper__`, na.rm = TRUE)
  minY <- min(dat$`lower__`, na.rm = TRUE)
  carn_plot <- ggplot(data = subset(dat, hostGroup=="carnivores"), aes_string(x = "combIUCN", y = "estimate__", color = "groupSizeCar")) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes_string(ymin = "lower__", ymax = "upper__"), width = 0, position = position_dodge(width = 0.3)) +
    ylim(minY, maxY) +
    scale_color_grey() +
    labs(y = y_nam) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          axis.title.x=element_blank())
  
  prim_plot <- ggplot(data = subset(dat, hostGroup=="primates"), aes_string(x = "combIUCN", y = "estimate__")) +
    geom_point(size = 3) +
    geom_errorbar(aes_string(ymin = "lower__", ymax = "upper__"), width = 0) +
    ylim(minY, maxY) +
    scale_color_grey() +
    theme_bw(base_size = 10) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank())
  
  ung_plot <- ggplot(data = subset(dat, hostGroup=="ungulates"), aes_string(x = "combIUCN", y = "estimate__")) +
    geom_point(size = 3) +
    geom_errorbar(aes_string(ymin = "lower__", ymax = "upper__"), width = 0) +
    ylim(minY, maxY) +
    scale_color_grey() +
    theme_bw(base_size = 10) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank())
 
  if(add_title) {
    carn_plot <- carn_plot + labs(title = "carnivores")
    prim_plot <- prim_plot + labs(title = "primates")
    ung_plot <- ung_plot + labs(title = "ungulates")
  }
  if(remove_xlab) {
    carn_plot <- carn_plot +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    prim_plot <- prim_plot +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    ung_plot <- ung_plot +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())}

  plot_row <- plot_grid(carn_plot, prim_plot, ung_plot, nrow = 1, rel_widths = c(1, 0.8, 0.8))
  return(plot_row)
}
margRich_fig <- FuncPlot1(dat = margRichBest_df, y_nam = "Parasite Richness", add_title = TRUE, remove_xlab = TRUE)
margParasTrans_fig <- FuncPlot1(dat = margParasTransBest_df, y_nam = "% Direct-Obligate", add_title = FALSE, remove_xlab = TRUE)
margParasType_fig <- FuncPlot1(dat = margParasTypeBest_df, y_nam = "% Micro-Parasite", add_title = FALSE, remove_xlab = FALSE)
fig1_final <- plot_grid(margRich_fig, margParasTrans_fig, margParasType_fig, ncol = 1, rel_heights = c(0.85, 1, 1))

pdf("./Results/FINAL-FIGS/FIG_condeffects_threat.pdf")
fig1_final
dev.off()

# ### FIG.1 ----
# # Effects of threat status on response (richness, % direct transmission, % micro-parasites) for simple vs. full models, separately for each host group
# 
# margRich <- margParasTrans <- margParasType <- list()
# 
# # Parasite richness
# 
# margRich$simple <- list(  # >>>>>> WAITING FOR CENTERED MARGINALS
#   carnivores = readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")],
#   primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")],
#   ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")])
# margRich_simple <- rbind.fill(margRich$simple)
# margRich_simple$hostGroup <- rep(names(margRich$simple), each = 2)
# margRich_simple$modType <- "simple"
# 
# margRich$complex <- list( 
#   carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "groupSizeCar", "estimate__", "lower__", "upper__")],
#   primates = readRDS("./Data/JPEK/full/full_brm_carngroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")],
#   ungulates = readRDS("./Data/JPEK/full/full_brm_carngroup_c_me.RDS")$combIUCN$data[c("combIUCN", "estimate__", "lower__", "upper__")])
# margRich_complex <- rbind.fill(margRich$complex)
# margRich_complex$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
# margRich_complex$modType <- "complex"
# 
# margRich_df <- rbind.fill(margRich_simple, margRich_complex)
# levels(margRich_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
# 
# # Parasite transmission
# margParasTrans$simple <- list(
#   carnivores = readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
#   primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
#   ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])
# 
# margParasTrans_simple <- rbind.fill(margParasTrans$simple)
# margParasTrans_simple$hostGroup <- rep(names(margParasTrans$simple), each = 2)
# margParasTrans_simple$modType <- "simple"
# 
# margParasTrans$complex <- list( 
#   carnivores = readRDS("./Data/JPEK/full/full_brm_carngroup_c_parastrans_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "parRichTransKnown", "groupSizeCar", "estimate__", "lower__", "upper__")],
#   primates = readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
#   ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_parastrans_me.RDS")$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])
# 
# margParasTrans_complex <- rbind.fill(margParasTrans$complex)
# margParasTrans_complex$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
# margParasTrans_complex$modType <- "complex"
# 
# margParasTrans_df <- rbind.fill(margParasTrans_simple, margParasTrans_complex) %>%
#   mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
#          `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
#          `upper__` = round(100*(`upper__`/parRichTransKnown), 1))
# levels(margParasTrans_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
#     
# # Parasite type
# margParasType$simple <- list(
#   carnivores = readRDS("./Data/JPEK/simple/simp_brm_carngroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
#   primates = readRDS("./Data/JPEK/simple/simp_brm_primgroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
#   ungulates = readRDS("./Data/JPEK/simple/simp_brm_unggroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])
# 
# margParasType_simple <- rbind.fill(margParasType$simple)
# margParasType_simple$hostGroup <- rep(names(margParasType$simple), each = 2)
# margParasType_simple$modType <- "simple"
# 
# margParasType$complex <- list( 
#   carnivores = readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_c_ParasType_me.RDS")$`combIUCN:groupSizeCar`$data[c("combIUCN", "parRich_alltypes", "groupSizeCar", "estimate__", "lower__", "upper__")],
#   primates = readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
#   ungulates = readRDS("./Data/JPEK/full/full_brm_unggroup_c_ParasType_me.RDS")$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])
# 
# margParasType_complex <- rbind.fill(margParasType$complex)
# margParasType_complex$hostGroup <- c(rep("carnivores", 4), rep("primates", 2), rep("ungulates", 2))
# margParasType_complex$modType <- "complex"
# 
# margParasType_df <- rbind.fill(margParasType_simple, margParasType_complex) %>%
#   mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
#          `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
#          `upper__` = round(100*(`upper__`/parRich_alltypes), 1))
# levels(margParasType_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
# 
# # Create figure
# FuncPlot1 <- function(dat, y_nam, remove_xlab, remove_faclab) {
#   p <- ggplot(data = dat) +
#     geom_errorbar(data = subset(dat, modType=="simple"), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0, color = "blue", position = position_nudge(x = -0.15)) +
#     geom_point(data = subset(dat, modType=="simple"), aes_string(x = "combIUCN", y = "estimate__"), shape = 21, color = "blue", fill = "white", size = 3, position = position_nudge(x = -0.15)) +
#     geom_errorbar(data = subset(dat, modType=="complex" & groupSizeCar %in% c(NA, "non_group")), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0, color = "black", position = position_nudge(x = 0.15)) +
#     geom_point(data = subset(dat, modType=="complex" & groupSizeCar %in% c(NA, "non_group")), aes_string(x = "combIUCN", y = "estimate__"), shape = 21, color = "black", fill = "white", size = 3, position = position_nudge(x = 0.15)) +
#     scale_y_continuous(breaks= pretty_breaks()) +
#     labs(y = y_nam) +
#     theme_bw(base_size = 10) +
#     facet_grid(. ~ hostGroup)
#   if(remove_faclab) {
#     p <- p + 
#       theme(strip.background = element_blank(),
#             strip.text.x = element_blank()
#             )
#   }
#   if(remove_xlab) {
#     p <- p +
#       theme(axis.title.x=element_blank(),
#             axis.text.x=element_blank(),
#             axis.ticks.x=element_blank())
#   } else {
#     p <- p +
#       labs(x = "Threat Status") 
#   }
# }
# 
# margParasTrans_fig <- FuncPlot1(dat = margParasTrans_df, y_nam = "% Direct-Obligate", remove_faclab = TRUE, remove_xlab = TRUE)
# margParasType_fig <- FuncPlot1(dat = margParasType_df, y_nam = "% Micro-Parasite", remove_faclab = TRUE, remove_xlab = FALSE)
# fig1_final <- plot_grid(margParasTrans_fig, margParasType_fig, ncol = 1, rel_heights = c(0.85, 1))
# 
# pdf("./Results/FINAL-FIGS/fig1_final.pdf")
# fig1_final
# dev.off()
  
### FIG.2 REDO----
# Marginal effects plots for IUCN:group size interaction FOR RICHNESS
# <<<<<< NEED TO REDO WITH CENTERED RESULTS
# carn_groupsize <- readRDS("./Data/JPEK/full/full_brm_carngroup_me.RDS")
carn_groupMarg_plot <- plot(carn_groupsize$`combIUCN:groupSizeCar`) +
  scale_color_grey() +
  scale_fill_grey() +
  ylim(0, 40) +
  labs(x = "Threat Status", y = "Parasite Richness", subtitle = "Carnivores") +
  theme_bw(base_size = 10) +
  scale_x_discrete(limits = c("not_threatened", "threatened"), labels=c("threatened" = "T", "not_threatened" = "NT")) +
  theme(legend.position = c(0.75, 0.85),
        legend.title=element_blank())

# prim_groupsize <- readRDS("./Data/JPEK/full/full_brm_primgroup_me.RDS")
prim_groupMarg_plot <- plot(prim_groupsize$`logGroupSizePriUng:combIUCN`) +
  scale_color_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
  scale_fill_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
  ylim(0, 40) + # <<<<<<<<<<<<<
  labs(x = "log(Group Size)", subtitle = "Primates") +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.75, 0.85),
        legend.title=element_blank(),
        axis.title.y=element_blank())

# ung_groupsize <- readRDS("./Data/JPEK/full/full_brm_unggroup_me.RDS")
ung_groupMarg_plot <- plot(ung_groupsize$`logGroupSizePriUng:combIUCN`) +
  ylim(0, 40) + # <<<<<<<<<<<<<
  scale_color_manual(guide = guide_legend(reverse=TRUE),values = c("tomato2", "seagreen4"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
  scale_fill_manual(guide = guide_legend(reverse=TRUE),values = c("tomato2", "seagreen4"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
  labs(x = "log(Group Size)", subtitle = "Ungulates") +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.75, 0.85),
        legend.title=element_blank(),
        axis.title.y=element_blank()) 
#         
fig2_final <- plot_grid(carn_groupMarg_plot, prim_groupMarg_plot, ung_groupMarg_plot, nrow = 1, rel_widths = c(1, 0.85, 0.85))

pdf("./Results/FINAL-FIGS/FIG_groupsize_rich.pdf", width=7, height=3.5)
fig2_final
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

