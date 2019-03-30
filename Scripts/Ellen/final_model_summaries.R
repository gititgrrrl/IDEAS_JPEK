### FINDINGS:
# > For PRIMATES 'TOTAL RICHNESS' & '% MICROPARASITES', simple model is moderately better than full; for '% CLOSE PARASITES', full is moderately better than simple. For carnivores and ungulates, full model is better than simple for all response types. HOST TRAITS MODULATE THE HOST THREAT STATUS-PARASITE RELATIONSHIP FOR CARNIVORES AND UNGULATES, MORE THAN FOR PRIMATE SPECIES.
# > FOR THREATENED SPECIES IN ALL HOST GROUPS, LARGER GROUPS HAVE LOWER PARASITE RICHNESS. PATTERNS FOR NON-THREATENED SPECIES DIFFER BY HOST GROUP.
# > % CLOSELY TRANSMITTED parasites and % MICROPARASITES are significantly lower for threatened UNGULATES than for non-threatened ungulates. For carnivores and primates, we don't see an effect of host threat status on these response variables. In other words, PARASITE TRAITS DIFFER BTWN THREATENED AND NON-THREATENED SPECIES ONLY FOR UNGULATES. 
# FOR MORE SOCIAL SPECIES (LARGER GROUP SIZE) REGARDLESS OF HOST GROUP, THREATENED SPECIES ARE ASSOCIATED WITH LOWER PARASITE RICHNESS AND LOWER % CLOSELY TRANSMITTED PARASITES.

rm(list=ls())

### LOAD PACKAGES ----
packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "plyr", "brms")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})

### READ IN FILES ----
FuncFilesIn <- function(path_nam, pattern_nam) {
  filenames <- list.files(path = path_nam, pattern = pattern_nam) # <<<<<<<<<< CHANGE AS NEEDED
  
  files_list <- list()
  outname <- NA
  for(i in 1:length(filenames)){
    outname <- unlist(strsplit(x = filenames[i], split = '.RDS'))[1]
    files_list[[i]] <- readRDS(paste0(path_nam, filenames[i]))
    names(files_list)[[i]] <- outname
    assign(x = outname, value = files_list[[i]])
  }
  return(files_list)
  }

# ... full models
fullmod_list <- FuncFilesIn(path_nam = "./FINAL/full/", pattern_nam ="_mod.RDS$")

# ... simple models
simplemod_list <- FuncFilesIn(path_nam = "./FINAL/simple/", pattern_nam ="_mod.RDS$")

# ... full model marginal effects
fullmarginal_list <- FuncFilesIn(path_nam = "./FINAL/full/", pattern_nam ="_marginal.RDS$")

# ... simple model marginal effects
simplemarginal_list <- FuncFilesIn(path_nam = "./FINAL/simple/", pattern_nam ="_marginal.RDS$")

# ... full model group size marginal effects
fullmarginalgroup_list <- FuncFilesIn(path_nam = "./FINAL/full/", pattern_nam ="_marginal_groupsize.RDS$")


### INFLUENTIAL POINTS ----
# > For full & simple models, mostly good and a few okay
lapply(fullmod_list, loo)
lapply(simplemod_list, loo)

### KFOLD MODEL COMPARISON ----
FuncK <- function(mod_list, x_wts) {
  temp_list <- lapply(mod_list, function(x) as.numeric(c(x$kfold$estimates[c("kfoldic"),]))) # kfold ic & SE
  temp_df <- data.frame(do.call("rbind", temp_list))
  temp_summary <- cbind(names(mod_list), temp_df) %>%
    mutate_if(is.factor, as.character)
  rownames(temp_summary) <- NULL
  colnames(temp_summary) = c("model", "kfoldic", "se(kfoldic)")
  Kfoldwts_df <- data.frame(x_wts)
  Kfoldwts_df$model <- rownames(Kfoldwts_df)
  rownames(Kfoldwts_df) <- NULL
  colnames(Kfoldwts_df) = c("weight", "model")
  Kfold_summary <- temp_summary %>% # <<<<< THIS IS THE K-FOLD TABLE
    full_join(Kfoldwts_df, by = "model") %>%
    arrange(kfoldic) %>%
    mutate_if(is.numeric, round, 3)
  return(Kfold_summary)
}

IC_list <- list()
for (r in c("totrich", "parastrans", "parastype")) {
  for (g in c("carngroup", "primgroup", "unggroup")) {
    
    simple_nam <- paste0("simple_", r, "_", g)
    simple_mod <- simplemod_list[[paste0("simple_", r, "_", g, "_mod")]]
    full_nam <- paste0("full_", r, "_", g)
    full_mod <- fullmod_list[[paste0("full_", r, "_", g, "_mod")]]
    mod_list = list(simple_mod, full_mod)
    names(mod_list) = c(simple_nam, full_nam)
    x_wts <- model_weights(simple_mod, full_mod, weights = "kfold")
    names(x_wts) = c(simple_nam, full_nam)
    (IC_list[[paste0(r, "_", g)]] <- FuncK(mod_list = mod_list, x_wts = x_wts))
  }
}

ICtab <- data.frame(do.call("rbind", IC_list))
saveRDS(ICtab, "./FINAL/summaries/IC_tab.RDS")

### MARGINAL EFFECTS PLOTS ----
# Other predictors are set at mean or, FOR CATEGORICAL, AT REFERENCE VALUE

# Plot of threat status marginal effect for simple vs. full ----
# Parasite richness
margRich <- list()
margRich$simple <- list(
  carnivores = simplemarginal_list$simple_totrich_carngroup_marginal$combIUCN[c("combIUCN", "estimate__", "lower__", "upper__")],
  primates = simplemarginal_list$simple_totrich_primgroup_marginal$combIUCN[c("combIUCN", "estimate__", "lower__", "upper__")],
  ungulates = simplemarginal_list$simple_totrich_unggroup_marginal$combIUCN[c("combIUCN", "estimate__", "lower__", "upper__")])
margRich_simple <- rbind.fill(margRich$simple)
margRich_simple$hostGroup <- rep(names(margRich$simple), each = 2)
margRich_simple$modType <- "simple"

margRich$full <- list(
  carnivores = fullmarginal_list$full_totrich_carngroup_marginal$combIUCN[c("combIUCN", "estimate__", "lower__", "upper__")],
  primates = fullmarginal_list$full_totrich_primgroup_marginal$combIUCN[c("combIUCN", "estimate__", "lower__", "upper__")],
  ungulates = fullmarginal_list$full_totrich_unggroup_marginal$combIUCN[c("combIUCN", "estimate__", "lower__", "upper__")])
margRich_full <- rbind.fill(margRich$full)
margRich_full$hostGroup <- rep(names(margRich$full), each = 2)
margRich_full$modType <- "full"

margRich_df <- rbind.fill(margRich_simple, margRich_full)
margRich_df$respType <- "Parasite Richness"
levels(margRich_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
margRich_df$modType <- factor(margRich_df$modType, levels = c("simple", "full"))

# Parasite transmission
margParasTrans <- list()
margParasTrans$simple <- list(
  carnivores = simplemarginal_list$simple_parastrans_carngroup_marginal$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  primates = simplemarginal_list$simple_parastrans_primgroup_marginal$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  ungulates = simplemarginal_list$simple_parastrans_unggroup_marginal$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])

margParasTrans_simple <- rbind.fill(margParasTrans$simple)
margParasTrans_simple$hostGroup <- rep(names(margParasTrans$simple), each = 2)
margParasTrans_simple$modType <- "simple"

margParasTrans$full <- list(
  carnivores = fullmarginal_list$full_parastrans_carngroup_marginal$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  primates = fullmarginal_list$full_parastrans_primgroup_marginal$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")],
  ungulates = fullmarginal_list$full_parastrans_unggroup_marginal$combIUCN$data[c("combIUCN", "parRichTransKnown", "estimate__", "lower__", "upper__")])

margParasTrans_full <- rbind.fill(margParasTrans$full)
margParasTrans_full$hostGroup <- rep(names(margParasTrans$full), each = 2)
margParasTrans_full$modType <- "full"

margParasTrans_df <- rbind.fill(margParasTrans_simple, margParasTrans_full) %>%
  mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
         `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
         `upper__` = round(100*(`upper__`/parRichTransKnown), 1))
margParasTrans_df$parRichTransKnown <- NULL
margParasTrans_df$respType <- "% Closely Transmitted"
levels(margParasTrans_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
margParasTrans_df$modType <- factor(margParasTrans_df$modType, levels = c("simple", "full"))
        
# Parasite type
margParasType <- list()
margParasType$simple <- list(
  carnivores = simplemarginal_list$simple_parastype_carngroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  primates = simplemarginal_list$simple_parastype_primgroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  ungulates = simplemarginal_list$simple_parastype_unggroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])

margParasType_simple <- rbind.fill(margParasType$simple)
margParasType_simple$hostGroup <- rep(names(margParasType$simple), each = 2)
margParasType_simple$modType <- "simple"

margParasType$full <- list(
  carnivores = fullmarginal_list$full_parastype_carngroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  primates = fullmarginal_list$full_parastype_primgroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
  ungulates = fullmarginal_list$full_parastype_unggroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])

margParasType_full <- rbind.fill(margParasType$full)
margParasType_full$hostGroup <- rep(names(margParasType$full), each = 2)
margParasType_full$modType <- "full"

margParasType_df <- rbind.fill(margParasType_simple, margParasType_full) %>%
  mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
         `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
         `upper__` = round(100*(`upper__`/parRich_alltypes), 1))
margParasType_df$parRich_alltypes <- NULL
margParasType_df$respType <- "% Microparasites"
levels(margParasType_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")
margParasType_df$modType <- factor(margParasType_df$modType, levels = c("simple", "full"))

marginal_df <- rbind.fill(margRich_df, margParasTrans_df, margParasType_df)
saveRDS(marginal_df, "./FINAL/summaries/marginal_df.RDS")

refplot_marginal_simpleVfull <- 
  ggplot(data = marginal_df, aes(color = modType, group = modType)) +
  geom_errorbar(data = subset(marginal_df, modType=="simple"), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0, color = "skyblue3", position = position_nudge(x = -0.15)) +
  geom_line(data = subset(marginal_df, modType=="simple"), aes_string(x = "combIUCN", y = "estimate__"), color = "skyblue3", size = 0.1, position = position_nudge(x = -0.15)) +
  geom_point(data = subset(marginal_df, modType=="simple"), aes_string(x = "combIUCN", y = "estimate__"), color = "skyblue3", size = 2, position = position_nudge(x = -0.15)) +
  geom_errorbar(data = subset(marginal_df, modType=="full"), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0, color = "black", position = position_nudge(x = 0.15)) + 
  geom_line(data = subset(marginal_df, modType=="full"), aes_string(x = "combIUCN", y = "estimate__"), color = "black", size = 0.1, position = position_nudge(x = 0.15)) +
  geom_point(data = subset(marginal_df, modType=="full"), aes_string(x = "combIUCN", y = "estimate__"), color = "black", size = 2, position = position_nudge(x = 0.15)) + 
  scale_y_continuous(breaks= pretty_breaks()) +
    labs(x = "Threat Status", y = "Expected Value (95% CI)", subtitle = "(blue is simple; black is full)" +
  theme_bw(base_size = 10) +
  facet_grid(respType ~ hostGroup, scales = "free_y")
saveRDS(refplot_marginal_simpleVfull, "./FINAL/summaries/refplot_marginal_simpleVfull.RDS")

pdf("./FINAL/summaries/refplot_marginal_simpleVfull.pdf")
refplot_marginal_simpleVfull
dev.off()

# Same as above, but overlapping for Joy ----
refplot_marginal_simpleVfull_JOY <-
  ggplot(data = marginal_df, aes(color = modType, group = modType)) +
  geom_errorbar(data = subset(marginal_df, modType=="simple"), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0.1, color = "turquoise3", size = 0.2) +
  geom_line(data = subset(marginal_df, modType=="simple"), aes_string(x = "combIUCN", y = "estimate__"), color = "turquoise3", size = 0.1) +
  geom_point(data = subset(marginal_df, modType=="simple"), aes_string(x = "combIUCN", y = "estimate__"), color = "turquoise3", size = 2) +
  geom_errorbar(data = subset(marginal_df, modType=="full"), aes_string(x = "combIUCN", ymin = "lower__", ymax = "upper__"), width = 0.1, color = "red", size = 0.2) +
  geom_line(data = subset(marginal_df, modType=="full"), aes_string(x = "combIUCN", y = "estimate__"), color = "red", size = 0.1) +
  geom_point(data = subset(marginal_df, modType=="full"), aes_string(x = "combIUCN", y = "estimate__"), color = "red", size = 2) +
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = "Threat Status", y = "Expected Value (95% CI)", subtitle = "blue is simple model; red is full model") +
  theme_bw(base_size = 10) +
  facet_grid(respType ~ hostGroup, scales = "free_y")
png("./FINAL/summaries/refplot_marginal_simpleVfull_JOY.png")
refplot_marginal_simpleVfull_JOY
dev.off()

# Plot of threat status marginal effect for best model ----
# The only one for which SIMPLE is best is 'parastype_primgroup'
# Other predictors are set at mean or, FOR CATEGORICAL, AT REFERENCE VALUE
margParasType$best <- list(
carnivores = fullmarginal_list$full_parastype_carngroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
primates = simplemarginal_list$simple_parastype_primgroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")],
ungulates = fullmarginal_list$full_parastype_unggroup_marginal$combIUCN$data[c("combIUCN", "parRich_alltypes", "estimate__", "lower__", "upper__")])

margParasType_best_df <- rbind.fill(margParasType$best)
margParasType_best_df$hostGroup <- rep(names(margParasType$best), each = 2)
margParasType_best_df$modType <- "best"

margParasType_best_df %<>%
  mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
         `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
         `upper__` = round(100*(`upper__`/parRich_alltypes), 1))
margParasType_best_df$parRich_alltypes <- NULL
margParasType_best_df$respType <- "% Microparasites"
levels(margParasType_best_df$combIUCN) <- list(NT = "not_threatened", T = "threatened")

marginal_best_df <- rbind(subset(margParasTrans_df, modType == "full"), margParasType_best_df)
saveRDS(marginal_best_df, "./FINAL/summaries/marginal_best_df.RDS")

plot_marginal_best <- 
  ggplot(data = marginal_best_df) +
  geom_errorbar(data = subset(marginal_best_df, combIUCN=="NT"), aes_string(x = "hostGroup", ymin = "lower__", ymax = "upper__"), width = 0, position = position_nudge(x = -0.15)) +
  geom_point(data = subset(marginal_best_df, combIUCN=="NT"), aes_string(x = "hostGroup", y = "estimate__"), size = 2, position = position_nudge(x = -0.15)) +
  geom_errorbar(data = subset(marginal_best_df, combIUCN=="T"), aes_string(x = "hostGroup", ymin = "lower__", ymax = "upper__"), linetype = "dotted", width = 0, position = position_nudge(x = 0.15)) + # threatened is dotted line
  geom_point(data = subset(marginal_best_df, combIUCN=="T"), aes_string(x = "hostGroup", y = "estimate__"), shape = 1, size = 2, position = position_nudge(x = 0.15)) + # threatened is open circle
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = "Host Group", y = "Expected Value (95% CI)") +
  theme_bw(base_size = 10) +
  facet_grid(respType ~ ., scales = "free_y")
saveRDS(plot_marginal_best, "./FINAL/summaries/plot_marginal_best.RDS")

pdf("./FINAL/summaries/plot_marginal_best.pdf")
plot_marginal_best
dev.off()

# Plot of group size:threat status marginal effect for full model ----
groupmarg_list <- list()
groupmarg_list$totrich <- list(
  carnivores = fullmarginalgroup_list[["full_totrich_carngroup_marginal_groupsize"]]$`groupSizeCar:combIUCN`,
  primates = fullmarginalgroup_list[["full_totrich_primgroup_marginal_groupsize"]]$`logGroupSizePriUng_c:combIUCN`,
  ungulates = fullmarginalgroup_list[["full_totrich_unggroup_marginal_groupsize"]]$`logGroupSizePriUng_c:combIUCN`)

groupmarg_list$parastrans <- list(
  carnivores = fullmarginalgroup_list[["full_parastrans_carngroup_marginal_groupsize"]]$`groupSizeCar:combIUCN` %>%
    mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
           `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
           `upper__` = round(100*(`upper__`/parRichTransKnown), 1)),
  primates = fullmarginalgroup_list[["full_parastrans_primgroup_marginal_groupsize"]]$`logGroupSizePriUng_c:combIUCN`  %>%
    mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
           `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
           `upper__` = round(100*(`upper__`/parRichTransKnown), 1)),
  ungulates = fullmarginalgroup_list[["full_parastrans_unggroup_marginal_groupsize"]]$`logGroupSizePriUng_c:combIUCN`  %>%
    mutate(`estimate__` = round(100*(`estimate__`/parRichTransKnown), 1),
           `lower__` = round(100*(`lower__`/parRichTransKnown), 1),
           `upper__` = round(100*(`upper__`/parRichTransKnown), 1)))

groupmarg_list$parastype <- list(
  carnivores = fullmarginalgroup_list[["full_parastype_carngroup_marginal_groupsize"]]$`groupSizeCar:combIUCN` %>%
    mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
           `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
           `upper__` = round(100*(`upper__`/parRich_alltypes), 1)),
  primates = fullmarginalgroup_list[["full_parastype_primgroup_marginal_groupsize"]]$`logGroupSizePriUng_c:combIUCN`  %>%
    mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
           `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
           `upper__` = round(100*(`upper__`/parRich_alltypes), 1)),
  ungulates = fullmarginalgroup_list[["full_parastype_unggroup_marginal_groupsize"]]$`logGroupSizePriUng_c:combIUCN`  %>%
    mutate(`estimate__` = round(100*(`estimate__`/parRich_alltypes), 1),
           `lower__` = round(100*(`lower__`/parRich_alltypes), 1),
           `upper__` = round(100*(`upper__`/parRich_alltypes), 1)))

plot_groupmarg_list <- list()
for (r in c("totrich", "parastrans", "parastype")) {
  carn <- ggplot(
    groupmarg_list[[r]]$carnivores, aes_string(x="groupSizeCar", y = "estimate__", color = "effect2__")) + 
    geom_errorbar(data = subset(groupmarg_list[[r]]$carnivores, effect2__=="not_threatened"), aes_string(ymin = "lower__", ymax = "upper__"), width = 0, position = position_nudge(x = -0.15)) +
    geom_point(data = subset(groupmarg_list[[r]]$carnivores, effect2__=="not_threatened"), size = 1, position = position_nudge(x = -0.15)) +
    geom_errorbar(data = subset(groupmarg_list[[r]]$carnivores, effect2__=="threatened"), aes_string(ymin = "lower__", ymax = "upper__"), width = 0, position = position_nudge(x = 0.15)) + 
    geom_point(data = subset(groupmarg_list[[r]]$carnivores, effect2__=="threatened"), size = 1, position = position_nudge(x = 0.15)) +
    scale_y_continuous(breaks= pretty_breaks()) +
    scale_color_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
    scale_fill_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
    labs(y = ifelse(r=="totrich", "Expected Richness", ifelse(r=="parastrans", "Expected % Close Trans.", ifelse(r=="parastype", "Expected % Micro.", NA))), subtitle = if(r=="totrich") {"Carnivores"}) +
    theme_bw(base_size = 7) +
    theme(legend.title = element_blank()) +
    {if(r=="totrich") {
      theme(legend.position = c(0.85, 0.95))
    } else {
      theme(legend.position = "none")}} +
    {if(r %in% c("totrich", "parastrans")) {
      theme(axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank())}}
      
  prim <- ggplot(groupmarg_list[[r]]$primates, aes_string(x="logGroupSizePriUng_c", y = "estimate__", color = "effect2__")) + geom_line() + geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = effect2__), alpha=0.1, color = NA) +
    scale_color_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
    scale_fill_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
    labs(x = if(r=="parastype") {"log(Group Size)"}, subtitle = if(r=="totrich") {"Primates"}) +
    theme_bw(base_size = 7) +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank()) +
    {if(r=="totrich") {
      theme(legend.position = c(0.85, 0.95))
    } else {
      theme(legend.position = "none")}} +
      {if(r %in% c("totrich", "parastrans")) {
        theme(axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank())}}
  
  ung <- ggplot(groupmarg_list[[r]]$ungulates, aes_string(x="logGroupSizePriUng_c", y = "estimate__", color = "effect2__")) + geom_line() + geom_ribbon(aes(ymin=lower__, ymax=upper__, fill = effect2__), alpha=0.1, color = NA) +
    scale_color_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
    scale_fill_manual(values = c("seagreen4", "tomato2"), labels = c("threatened" = " T", "not_threatened" = " NT")) +
    labs(x = if(r=="parastype") {"log(Group Size)"}, subtitle = if(r=="totrich") {"Ungulates"}) +
    theme_bw(base_size = 7) +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank()) +
    {if(r=="totrich") {
      theme(legend.position = c(0.85, 0.95))
    } else {
      theme(legend.position = "none")}} +
      {if(r %in% c("totrich", "parastrans")) {
        theme(axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank())}}
  
  plot_groupmarg_list[[r]] <- plot_grid(carn, prim, ung, rel_widths = c(1.1, 1, 1), nrow = 1)
}
plot_groupmarginal_full <- plot_grid(plotlist =  plot_groupmarg_list, rel_heights = c(1, 0.9, 1.05), ncol=1)
     
saveRDS(plot_groupmarginal_full, "./FINAL/summaries/plot_groupmarginal_full.RDS")

pdf("./FINAL/summaries/plot_groupmarginal_full.pdf")
plot_groupmarginal_full
dev.off()
 
# GENERATE COEFFICIENT PLOTS ----
# NOTE: Only look at coefficients that both simple and full models estimate.
# UPSHOT: For all three response variables, the full vs. simple (using same data) coefficients are pretty different, especially for primates. But the confidence intervals for full model coefficients are quite large, so often still overlap the simple model coefficients.
FuncPlotCoef <- function(mod_list, terms_vec = NULL, exclude_terms = FALSE, plot_title) {
  # Coefficient plots, comparing across models
  #
  # Args:
  #   mod_list:  List of (named) fitted regression models
  #   terms_vec: The model terms to plot. To plot ALL terms, leave as NULL
  #   exclude_terms: TRUE if parameters in 'terms_vec' should be excluded from the plot, rather than included
  #   plot_title: Title for plot
  # Returns:
  #   For multiple models, plot of coefficient estimates with 95% credible intervals. Each term is in a separate facet.
  # 
  df <- mod_list %>% purrr::map(function(.) { 
    broom::tidy(., conf.int = TRUE, par_type = "non-varying") }) %>% 
    dplyr::bind_rows(.id = "model") 
  if(exclude_terms) {
    df <- df %>% 
      filter(!term %in% terms_vec)
  } else {
    if(!is.null(terms_vec)) {
      df <- df %>%
        filter(term %in% terms_vec)
    }
  }
  p <- df %>% 
    ggplot(aes(term, estimate, ymin = lower, ymax = upper, color = model)) +
    geom_pointrange(position = position_dodge(width = 0.3)) + coord_flip() +
    labs(y = "Coefficient estimate (95% CI)", x = "Model parameter", title = plot_title) +
    geom_hline(yintercept = 0,  linetype = "dotted") +
    theme_bw() +
    theme(legend.position="top")
  png(paste0("./FINAL/summaries/refplot_coef_", plot_title, ".png"))
  plot(p)
  dev.off()
}

for(r in c("totrich", "parastrans", "parastype")) {
  for(g in c("carngroup", "unggroup", "primgroup")) {
    FuncPlotCoef(mod_list = list(simple_mod = simplemod_list[[paste0("simple_", r, "_", g, "_mod")]], full_mod = fullmod_list[[paste0("full_", r, "_", g, "_mod")]]), exclude_terms = TRUE, terms_vec ="Intercept", plot_title = paste0(r, "_", g))
  }
}

FuncPlotCoef(mod_list = list(simple_fulldat = simple_mod_fulldat, full = full_mod), terms_vec = c("Intercept", "combIUCNthreatened", "hostGroupprimates", "hostGroupungulates", "logNumHostCitations", "combIUCNthreatened:hostGroupprimates", "combIUCNthreatened:hostGroupungulates", "combIUCNthreatened:logNumHostCitations"))

# Response = % close parasite transmission
FuncPlotCoef(mod_list = list(simple_parastrans_fulldat = simple_mod_parastrans_fulldat, full_parastrans = full_mod_parastrans), terms_vec = c("Intercept", "combIUCNthreatened", "hostGroupprimates", "hostGroupungulates", "logNumHostCitations", "combIUCNthreatened:hostGroupprimates", "combIUCNthreatened:hostGroupungulates", "combIUCNthreatened:logNumHostCitations"))

# Response = % micro parasite
FuncPlotCoef(mod_list = list(simple_parastype_fulldat = simple_mod_parastype_fulldat, full_parastype = full_mod_parastype), terms_vec = c("Intercept", "combIUCNthreatened", "hostGroupprimates", "hostGroupungulates", "logNumHostCitations", "combIUCNthreatened:hostGroupprimates", "combIUCNthreatened:hostGroupungulates", "combIUCNthreatened:logNumHostCitations"))


# ### MODEL PREDICTIONS ----
# FuncPredVObs <- function(dat_df, predict_list, resp_colnam, obsID_colnam, color_colnam, coloraxt_title, facet_colnam) {
#   names(dat_df)[names(dat_df) == resp_colnam] <- "resp"
#   names(dat_df)[names(dat_df) == obsID_colnam] <- "group_var"
#   if(!is.null(color_colnam)) {
#     names(dat_df)[names(dat_df) == color_colnam] <- "color_var"
#   }
#   if(!is.null(facet_colnam)) {
#     names(dat_df)[names(dat_df) == facet_colnam] <- "facet_var"
#   }
#   
#   summary_list <- p1_list <- vector("list", length = length(predict_list))
#   
#   for (i in 1:length(predict_list)) {
#     predict_out <- predict_list[[i]]
#     
#     # Model predicted responses w/95% CI
#     predict_out %<>%
#       as_tibble() %>%
#       mutate(pred_Q2.5 = Q2.5,
#              pred_Q97.5 = Q97.5)
#     
#     summary_list[[i]] <- cbind(dat_df, predict_out) %>%
#       mutate(pred_diff_obs = Estimate - resp,
#              CI_overlap = resp>=pred_Q2.5 & resp <=pred_Q97.5) # does the 95% CI overlap the observed value?
#     
#   } # end of predict_list
#   
#   axes_min = floor(min(sapply(summary_list, function(x) min(subset(x, select = c(pred_Q2.5, resp)), na.rm=TRUE))))
#   axes_max = ceiling(max(sapply(summary_list, function(x) max(subset(x, select = c(pred_Q97.5, resp)), na.rm=TRUE))))
#   
#   for (i in 1:length(predict_list)) {
#     plot_dat <- summary_list[[i]]
#     # Plot 1
#     p1 <- ggplot(data = plot_dat, aes(x = resp, y = Estimate)) + 
#       geom_abline(intercept = 0, slope = 1, 
#                   linetype = "dashed", color = "gray") +
#       labs(x =  "Observed parasite richness", y = "Predicted parasite richness", subtitle = names(predict_list[i])) +
#       scale_color_brewer(palette = "Dark2") +
#       scale_fill_brewer(palette = "Dark2") +
#       theme_bw(base_size = 10) +
#       theme(legend.position="top") 
#     
#     if(!is.null(color_colnam)) {
#       p1 <- p1 +       
#         geom_errorbar(aes(ymin = pred_Q2.5, ymax = pred_Q97.5, color = color_var), width = 0) +
#         geom_point(aes(color = color_var), size = 2, alpha = 0.6) +
#         labs(color = coloraxt_title)
#     } else {
#       p1 <- p1 +       
#         geom_errorbar(aes(ymin = pred_Q2.5, ymax = pred_Q97.5), width = 0, color = "lightblue") +
#         geom_point(size = 2, shape = 1)
#     }
#     
#     if(!is.null(facet_colnam)) {
#       p1 <- p1 +  
#         facet_wrap(.~facet_var, scales = "free")
#     }
#     
#     p1_list[[i]] <- p1
#   }
#   
#   # combine p1's
#   title_p1 <- ggdraw() +
#     draw_label("Predicted vs. observed", fontface = "bold", size = 12)
#   subtitle_p1 <- ggdraw() +
#     draw_label("Vertical lines are 95% confidence intervals of the average prediction.\nThe gray dashed line shows perfect prediction.", size = 10)
#   
#   final_p1 <- plot_grid(title_p1, subtitle_p1, plotlist= p1_list, ncol = 1, rel_heights = c(0.07, 0.13, rep(1, length(p1_list))))
#   
#   return_list <- list(summary_list, final_p1)
# }
# 
# # RESPONSE IS PARASITE RICHNESS ----
# simple_brm_all_predict_fulldat <- readRDS("./FINAL/simple/simple_brm_all_predict_fulldat.RDS")
# full_brm_all_predict <- readRDS("./FINAL/full/full_brm_all_predict.RDS")
# 
# temp_out <- FuncPredVObs(dat_df = fullDat, predict_list = list(`SIMPLE MODEL (using full model data) FOR PARASITE RICHNESS` = simple_brm_all_predict_fulldat, `FULL MODEL FOR PARASITE RICHNESS` = full_brm_all_predict), resp_colnam = "parRich", obsID_colnam = "hostName", color_colnam = "combIUCN", coloraxt_title = "Host IUCN status", facet_colnam = "hostGroup")
# saveRDS(temp_out[[1]], "./Results/model-predictions/predict_observ_summary_simple_fulldat_totrich.RDS")
# pdf(paste0("./Results/model-predictions/predict_observ_plot_simple_fulldat_totrich.pdf")) 
# temp_out[[2]]
# dev.off()
# 
# # RESPONSE IS PARASITE TRANSMISSION ----
# simple_brm_parastrans_predict_fulldat <- readRDS("./FINAL/simple/simple_brm_parastrans_predict_fulldat.RDS")
# full_brm_parastrans_predict <- readRDS("./FINAL/full/full_brm_parastrans_predict.RDS")
# 
# temp_out <- FuncPredVObs(dat_df = fullDat_parastrans, predict_list = list(`SIMPLE MODEL (using full model data) FOR PARASITE TRANSMISSION` = simple_brm_parastrans_predict_fulldat, `FULL MODEL FOR PARASITE TRANSMISSION` = full_brm_parastrans_predict), resp_colnam = "parRichCloseOnly", obsID_colnam = "hostName", color_colnam = "combIUCN", coloraxt_title = "Host IUCN status", facet_colnam = "hostGroup")
# saveRDS(temp_out[[1]], "./Results/model-predictions/predict_observ_summary_simple_fulldat_parastrans.RDS")
# pdf(paste0("./Results/model-predictions/predict_observ_plot_simple_fulldat_parastrans.pdf")) 
# temp_out[[2]]
# dev.off()
# 
# # RESPONSE IS PARASITE TYPE ----
# simple_brm_parastype_predict_fulldat <- readRDS("./FINAL/simple/simple_brm_parastype_predict_fulldat.RDS")
# full_brm_parastype_predict <- readRDS("./FINAL/full/full_brm_parastype_predict.RDS")
# 
# temp_out <- FuncPredVObs(dat_df = fullDat_parastype, predict_list = list(`SIMPLE MODEL (using full model data) FOR PARASITE TRANSMISSION` = simple_brm_parastype_predict_fulldat, `FULL MODEL FOR PARASITE TRANSMISSION` = full_brm_parastype_predict), resp_colnam = "parRich_micro", obsID_colnam = "hostName", color_colnam = "combIUCN", coloraxt_title = "Host IUCN status", facet_colnam = "hostGroup")
# saveRDS(temp_out[[1]], "./Results/model-predictions/predict_observ_summary_simple_fulldat_parastype.RDS")
# pdf(paste0("./Results/model-predictions/predict_observ_plot_simple_fulldat_parastype.pdf")) 
# temp_out[[2]]
# dev.off()

