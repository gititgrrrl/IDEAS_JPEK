### SCRIPT FOR MODEL DIAGNOSTICS

# STILL TO DO:
# > NEED TO DO MODEL SELECTION ON PARAS TYPE AND PARAS TRANS MODELS --WAITING FOR THE IC INFORMATION
# > NEED TO DO RESIDUAL DIAGNOSTICS --NEED RESIDUALS TO BE OBTAINED, USING THE FAST COMPUTER
# > For predict vs. observe plots, the axis limits for each host group should be the same for simple vs. full model (for that host group, but should be allowed to vary among host groups)--figure out how to set within facet levels
# > Figure out patterns in the observations that were not predicted well. Put in the group size information, for example
# > Identify which predictors in the full model probably were not important

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
options (mc.cores=parallel::detectCores ()) # Run chains on multiple cores# Model diagnostics

### READ IN DATA FOR PARASITE RICHNESS ----
allDat <- read_csv("Data/JPEK/script4.csv")
# ...carnivores
dat_rich_carn <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, groupSizeCar, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() # 148 records
dat_rich_carn <- dat_rich_carn[complete.cases(dat_rich_carn),] # only 74 records
# ...ungulates
dat_rich_ung <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() # 102 records
dat_rich_ung <- dat_rich_ung[complete.cases(dat_rich_ung),] # only 60 records
# ...primates
dat_rich_prim <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, groupSizePriUng, parRich, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() # 142 records
dat_rich_prim <- dat_rich_prim[complete.cases(dat_rich_prim),] # only 73 records

### READ IN DATA FOR PARASITE TRANSMISSION ----
allDat <- read_csv("Data/JPEK/script4.csv")
# ...carnivores
dat_parastrans_carn <- allDat %>%
  filter(hostGroup == "carnivores") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizeCar, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

dat_parastrans_carn <- dat_parastrans_carn[complete.cases(dat_parastrans_carn),]

# ...ungulates
dat_parastrans_ung <- allDat %>%
  filter(hostGroup == "ungulates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

dat_parastrans_ung <- dat_parastrans_ung[complete.cases(dat_parastrans_ung),] # only 60 records

# ...primates
dat_parastrans_prim <- allDat %>%
  filter(hostGroup == "primates") %>%
  select(hostName, parRichCloseOnly, parRichTransKnown, groupSizePriUng, logNumHostCitations, combIUCN, logHostSpeciesRange, logHostMass, hostMaxLifespan, absHostMeanLat) %>%
  mutate(logGroupSizePriUng = log(groupSizePriUng)) %>%
  distinct() %>%
  filter(parRichTransKnown > 0) # 139

dat_parastrans_prim <- dat_parastrans_prim[complete.cases(dat_parastrans_prim),]

### READ IN MODELS ----
# Parasite richness
# ...allInt
mod_rich_allInt_carn <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_all.RDS")
mod_rich_allInt_ung <- readRDS("./Data/JPEK/allInt/allInt_brm_unggroup_all.RDS")
mod_rich_allInt_prim <- readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_all.RDS") 
# ...full
mod_rich_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_all.RDS")
mod_rich_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_all.RDS")
mod_rich_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_all.RDS") 
# ...simple
mod_rich_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_all.RDS")
mod_rich_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_all.RDS") 
mod_rich_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_all.RDS")

# Parasite transmission
# ...allInt
mod_parastrans_allInt_carn <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_parastrans.RDS")
mod_parastrans_allInt_ung <- readRDS("./Data/JPEK/allInt/allInt_brm_unggroup_parastrans.RDS")
mod_parastrans_allInt_prim <- readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_parastrans.RDS")
# ...full
mod_parastrans_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_parastrans.RDS")
mod_parastrans_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_parastrans.RDS")
mod_parastrans_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_parastrans.RDS")
# ...simple
mod_parastrans_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_parastrans.RDS")
mod_parastrans_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_parastrans.RDS")
mod_parastrans_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_parastrans.RDS")

# Parasite type
# ...allInt
mod_parastype_allInt_carn <- readRDS("./Data/JPEK/allInt/allInt_brm_carngroup_parastype.RDS")
mod_parastype_allInt_ung <- readRDS("./Data/JPEK/allInt/allInt_brm_unggroup_parastype.RDS")
mod_parastype_allInt_prim <- readRDS("./Data/JPEK/allInt/allInt_brm_primgroup_parastype.RDS")
# ...full
mod_parastype_full_carn <- readRDS("./Data/JPEK/full/full_brm_carngroup_parastype.RDS")
mod_parastype_full_ung <- readRDS("./Data/JPEK/full/full_brm_unggroup_parastype.RDS")
mod_parastype_full_prim <- readRDS("./Data/JPEK/full/full_brm_primgroup_parastype.RDS")
# ...simple
mod_parastype_simple_carn <- readRDS("./Data/JPEK/simple/simp_brm_carngroup_parastype.RDS")
mod_parastype_simple_ung <- readRDS("./Data/JPEK/simple/simp_brm_unggroup_parastype.RDS")
mod_parastype_simple_prim <- readRDS("./Data/JPEK/simple/simp_brm_primgroup_parastype.RDS")

### --- PP_CHECKS ----
# # SIMPLE RICHNESS
# summary(simpleBrm)
# 
# for(g in c("carnivores", "ungulates", "primates")) {
#   pdf(paste0("./Results/model-diagnostics/richness_dens_plot_simple_", g, ".pdf")) 
#   plot(pp_check(simpleBrm, type = "dens_overlay", nsamples = 300, 
#                 newdata = subset(simpleDat, hostGroup == g)) +
#          ggtitle(paste0(g, " density plot")) +
#          theme(plot.title = element_text(size = 12, face = "bold")))
#   dev.off()
# }
# pdf("./Results/model-diagnostics/richness_dens_plot_simple.pdf") 
# pp_check(simpleBrm, type = "violin_grouped", nsamples = 300, group = "combIUCN")
# dev.off()

# SIMPLE (ON FULL DATA) RICHNESS
summary(simpleBrm_fulldat)
for(g in c("carnivores", "ungulates", "primates")) {
  pdf(paste0("./Results/model-diagnostics/richness_dens_plot_simple_fulldat_", g, ".pdf")) 
  plot(pp_check(simpleBrm_fulldat, type = "dens_overlay", nsamples = 300, 
                newdata = subset(fullDat, hostGroup == g)) +
         ggtitle(paste0(g, " density plot")) +
         theme(plot.title = element_text(size = 12, face = "bold")))
  dev.off()
}
pdf("./Results/model-diagnostics/richness_dens_plot_simple_fulldat_iucn.pdf")
pp_check(simpleBrm_fulldat, type = "violin_grouped", nsamples = 300, group = "combIUCN")
dev.off()

# FULL RICHNESS
summary(fullBrm)
for(g in c("carnivores", "ungulates", "primates")) {
  pdf(paste0("./Results/model-diagnostics/richness_dens_plot_full_", g, ".pdf")) 
  plot(pp_check(fullBrm, type = "dens_overlay", nsamples = 300, 
                newdata = subset(fullDat, hostGroup == g)) +
         ggtitle(paste0(g, " density plot")) +
         theme(plot.title = element_text(size = 12, face = "bold")))
  dev.off()
}
pdf("./Results/model-diagnostics/richness_dens_plot_full_iucn.pdf")
pp_check(fullBrm, type = "violin_grouped", nsamples = 300, group = "combIUCN")
dev.off()

### INFLUENTIAL POINTS ----
loo(simpleBrm_fulldat) # ALL OK
loo(fullBrm$loo) # ALL OK

# >>>>>> CAN'T DO THESE UNTIL HAVE LOO-IC DATA
# loo(simpleBrm_parastrans_fulldat)
# loo(fullBrm_parastrans)
# 
# loo(simpleBrm_parastype_fulldat) 
# loo(fullBrm_parastype)

### RESIDUAL PLOTS ----
### >>>>>>>>>>>> PICK UP FROM HERE--ALSO CHECK RESIDUALS AGAINST THE TOTAL TRIALS
FuncPlotScatter <- function(dat_df, x_nam, y_nam, color_by = NULL, facet_by = NULL, free_y = fixed, trans_y = NULL, jitter = FALSE, jitter_width = NULL, jitter_height = NULL, add_loess = TRUE, plot_labs = list(x = NULL, y = NULL, color = NULL, title = NULL)) {
  # Function to generate scatterplot of data
  #
  # Args:
  #   dat_df:  A data frame with the raw data
  #   x_nam: Column name for x-axis variable
  #   y_nam: Column name for y-axis variable
  #   color_by:  Column name for variable to color by
  #   facet_by:  Column name for variable to facet by
  #   free_y:  Are scales shared across all facets?
  #   trans_y: Transformation to apply to y axis
  #   jitter: Jitter data points?
  #   jitter_width: Width to jitter, if jitter = TRUE
  #   jitter_height: Height to jitter, if jitter = TRUE
  #   add_loess:  Add loess smooth? (will be biased by censored data)
  #   plot_labs: list with x, y, title for plots
  #   
  # Returns:
  #   Scatterplot
  # 
  axes_min = min(dat_df[y_nam], na.rm=TRUE)
  axes_max = ceiling(max(dat_df[y_nam], na.rm=TRUE))
  
  if(is.null(color_by)) {
    p1 <- ggplot(data = dat_df, aes_string(x = x_nam, y = y_nam))
  } else {
    p1 <- ggplot(data = dat_df, aes_string(x = x_nam, y = y_nam, color = color_by), alpha = 0.5) + 
      scale_colour_brewer(palette = "Set1")
  }
  
  p1 <- p1 + 
  {if(jitter) {geom_jitter(width = jitter_width, height = jitter_height, size = 2)} else {geom_point(size = 2)}} + # can't use if-else within ggplot
    labs(title = plot_labs$title, x = plot_labs$x, y = plot_labs$y, color = plot_labs$color) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (!is.null(facet_by)) {
    p1 <- p1 +
      facet_grid(as.formula(paste0(facet_by, "~.")))
  }
  
  if(add_loess) {
    p1 <- p1 + geom_smooth(data = dat_df, aes_string(x = x_nam, y = y_nam), method = "loess")
  }
  
  if(!is.null(trans_y)) {
    p1 <- p1 + 
      scale_y_continuous(trans = trans_y, labels = function(x) as.character(round(x,2))) +
      labs(subtitle = paste0("Y-axis is on ", trans_y, "-scale."))
  }
  return(p1)
}

FuncPlotBox <- function(dat_df, x_nam, y_nam, trans_y = NULL, facet_nam = NULL, plot_labs = list(x = NULL, y = NULL, title = NULL)) {
  # Function to generate boxplots of data
  #
  # Args:
  #   dat_df:  A data frame with the raw data
  #   x_nam: Column name for x-axis variable
  #   y_nam: Column name for y-axis variable
  #   trans_y: Transformation to apply to y axis
  #   facet_nam: Column name for facet variable
  #   plot_labs: list with x, y, title for plots
  #   
  # Returns:
  #   Boxplot
  # 
  names(dat_df)[names(dat_df) == y_nam] <- "y"
  names(dat_df)[names(dat_df) == x_nam] <- "x"
  
  ggplot(data = dat_df, aes(x = as.factor(x), y = y)) +
    geom_boxplot() +
    labs(title = plot_labs$title) +
    theme_bw(base_size = 11) +
    theme(axis.title = element_text()) +
    labs(x = plot_labs$x, y = plot_labs$y)
  
  if(!is.null(facet_nam)) {
    p_box <- p_box +
      facet_grid(as.formula(paste(facet_nam, "~ .")))
  }
  
  if(!is.null(trans_y)) {
    p_box <- p_box + 
      scale_y_continuous(trans = trans_y, labels = function(x) as.character(round(x,2))) +
      labs(subtitle = paste0("Y-axis is on ", trans_y, "-scale."))
  }
  return(p_box)
}

# # FULL RICHNESS ----
# dat <- fullDat
# mod <- fullBrm
# mu <- readRDS("~/Google Drive/IDEAS_JPEK/Data/JPEK/full/full_brm_all_mu.RDS")
# resid_full_totrich <-cbind(dat, as_tibble(residuals(mod, type = "pearson")[, "Estimate"])) %>%
#   rename(PResid = value) %>%
#   cbind(as_tibble(fitted(mod)[, "Estimate"])) %>%
#   rename(Fitted = value) %>%
#   cbind(enframe(predict(mod)[, "Estimate"])) %>%
#   rename(Predicted = value)
# 
# # Residuals against XX
# FuncPlotScatter(dat_df = Resid_df, x_nam = "xx", y_nam = "PResid", color_by = "combIUCN", trans_y = NULL, add_loess = TRUE, plot_labs = list(x = "XXXXX", y = "Pearson Residual", title = NULL))
# 
# # Residual (boxplots) by YY and ZZ
# FuncPlotBox(dat_df = Resid_df, x_nam = "f_xx", y_nam = "PResid", trans_y = NULL, facet_nam = "zz", plot_labs = list(x = "xx", y = "Pearson Residual", title = NULL))

### -- MODEL PREDICTIONS ----
FuncPredVObs <- function(dat_df, predict_list, resp_colnam, obsID_colnam, color_colnam, coloraxt_title, facet_colnam) {
  names(dat_df)[names(dat_df) == resp_colnam] <- "resp"
  names(dat_df)[names(dat_df) == obsID_colnam] <- "group_var"
  if(!is.null(color_colnam)) {
    names(dat_df)[names(dat_df) == color_colnam] <- "color_var"
  }
  if(!is.null(facet_colnam)) {
    names(dat_df)[names(dat_df) == facet_colnam] <- "facet_var"
  }
  
  summary_list <- p1_list <- vector("list", length = length(predict_list))
  
  for (i in 1:length(predict_list)) {
    predict_out <- predict_list[[i]]
    
    # Model predicted responses w/95% CI
    predict_out %<>%
      as_tibble() %>%
      mutate(pred_Q2.5 = Q2.5,
             pred_Q97.5 = Q97.5)
    
    summary_list[[i]] <- cbind(dat_df, predict_out) %>%
      mutate(pred_diff_obs = Estimate - resp,
             CI_overlap = resp>=pred_Q2.5 & resp <=pred_Q97.5) # does the 95% CI overlap the observed value?
    
  } # end of predict_list
  
  axes_min = floor(min(sapply(summary_list, function(x) min(subset(x, select = c(pred_Q2.5, resp)), na.rm=TRUE))))
  axes_max = ceiling(max(sapply(summary_list, function(x) max(subset(x, select = c(pred_Q97.5, resp)), na.rm=TRUE))))
  
  for (i in 1:length(predict_list)) {
    plot_dat <- summary_list[[i]]
    # Plot 1
    p1 <- ggplot(data = plot_dat, aes(x = resp, y = Estimate)) + 
      geom_abline(intercept = 0, slope = 1, 
                  linetype = "dashed", color = "gray") +
      labs(x =  "Observed parasite richness", y = "Predicted parasite richness", subtitle = names(predict_list[i])) +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2") +
      theme_bw(base_size = 10) +
      theme(legend.position="top") 
    
    if(!is.null(color_colnam)) {
      p1 <- p1 +       
        geom_errorbar(aes(ymin = pred_Q2.5, ymax = pred_Q97.5, color = color_var), width = 0) +
        geom_point(aes(color = color_var), size = 2, alpha = 0.6) +
        labs(color = coloraxt_title)
    } else {
      p1 <- p1 +       
        geom_errorbar(aes(ymin = pred_Q2.5, ymax = pred_Q97.5), width = 0, color = "lightblue") +
        geom_point(size = 2, shape = 1)
    }
    
    if(!is.null(facet_colnam)) {
      p1 <- p1 +  
        facet_wrap(.~facet_var, scales = "free")
    }
    
    p1_list[[i]] <- p1
  }
  
  # combine p1's
  title_p1 <- ggdraw() +
    draw_label("Predicted vs. observed", fontface = "bold", size = 12)
  subtitle_p1 <- ggdraw() +
    draw_label("Vertical lines are 95% confidence intervals of the average prediction.\nThe gray dashed line shows perfect prediction.", size = 10)
  
  final_p1 <- plot_grid(title_p1, subtitle_p1, plotlist= p1_list, ncol = 1, rel_heights = c(0.07, 0.13, rep(1, length(p1_list))))
  
  return_list <- list(summary_list, final_p1)
}

# RICHNESS, SIMPLE (ON FULL DATA) COMPARED TO FULL ----
# Both simple & full models seem to perform about the same
simple_brm_all_predict_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_all_predict_fulldat.RDS")
full_brm_all_predict <- readRDS("./Data/JPEK/full/full_brm_all_predict.RDS")

temp_out <- FuncPredVObs(dat_df = fullDat, predict_list = list(`SIMPLE MODEL (using full model data) FOR PARASITE RICHNESS` = simple_brm_all_predict_fulldat, `FULL MODEL FOR PARASITE RICHNESS` = full_brm_all_predict), resp_colnam = "parRich", obsID_colnam = "hostName", color_colnam = "combIUCN", coloraxt_title = "Host IUCN status", facet_colnam = "hostGroup")
saveRDS(temp_out[[1]], "./Results/model-predictions/predict_observ_summary_simple_fulldat_totrich.RDS")
pdf(paste0("./Results/model-predictions/predict_observ_plot_simple_fulldat_totrich.pdf")) 
temp_out[[2]]
dev.off()

# PARASITE TRANSMISSION, SIMPLE (ON FULL DATA) COMPARED TO FULL ----
# # Full model is a little better for carnivores & primates, not a clear difference for ungulates
simple_brm_parastrans_predict_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastrans_predict_fulldat.RDS")
full_brm_parastrans_predict <- readRDS("./Data/JPEK/full/full_brm_parastrans_predict.RDS")

temp_out <- FuncPredVObs(dat_df = fullDat_parastrans, predict_list = list(`SIMPLE MODEL (using full model data) FOR PARASITE TRANSMISSION` = simple_brm_parastrans_predict_fulldat, `FULL MODEL FOR PARASITE TRANSMISSION` = full_brm_parastrans_predict), resp_colnam = "parRichCloseOnly", obsID_colnam = "hostName", color_colnam = "combIUCN", coloraxt_title = "Host IUCN status", facet_colnam = "hostGroup")
saveRDS(temp_out[[1]], "./Results/model-predictions/predict_observ_summary_simple_fulldat_parastrans.RDS")
pdf(paste0("./Results/model-predictions/predict_observ_plot_simple_fulldat_parastrans.pdf")) 
temp_out[[2]]
dev.off()

# PARASITE TYPE, SIMPLE (ON FULL DATA) COMPARED TO FULL ----
# Full model is better for carnivores & ungulates, not a clear difference for primates
simple_brm_parastype_predict_fulldat <- readRDS("./Data/JPEK/simple/simple_brm_parastype_predict_fulldat.RDS")
full_brm_parastype_predict <- readRDS("./Data/JPEK/full/full_brm_parastype_predict.RDS")

temp_out <- FuncPredVObs(dat_df = fullDat_parastype, predict_list = list(`SIMPLE MODEL (using full model data) FOR PARASITE TRANSMISSION` = simple_brm_parastype_predict_fulldat, `FULL MODEL FOR PARASITE TRANSMISSION` = full_brm_parastype_predict), resp_colnam = "parRich_micro", obsID_colnam = "hostName", color_colnam = "combIUCN", coloraxt_title = "Host IUCN status", facet_colnam = "hostGroup")
saveRDS(temp_out[[1]], "./Results/model-predictions/predict_observ_summary_simple_fulldat_parastype.RDS")
pdf(paste0("./Results/model-predictions/predict_observ_plot_simple_fulldat_parastype.pdf")) 
temp_out[[2]]
dev.off()


### -- KFOLD MODEL COMPARISON ----
# Model comparison ----
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

# RESPONSE IS PARASITE RICHNESS ----
# RESULTS ----
# ...carnivores, ALLINT IS BARELY BETTER
# model             kfoldic   se(kfoldic) weight
# rich_allInt_carn  590.107   47.866      0.688
# rich_full_carn    591.688   46.474      0.312
# rich_simple_carn  639.944   54.542      0.000
# ...ungulates, FULL WINS BUT SE IS LARGE. COULD NOT CALCULATE K-FOLD FOR 'allInt', but LOO-IC showed that model to be worse than 'full' or 'simple'
# model             kfoldic   se(kfoldic) weight
# rich_full_ung     616.040   77.176      1
# rich_simple_ung   640.017   45.459      0
# ...primates, SIMPLE & FULL PRETTY MUCH HAVE SAME KFOLD-IC (FULL IS BARELY BETTER)
# model             kfoldic   se(kfoldic) weight
# rich_full_prim    456.150   49.169      0.659
# rich_simple_prim  457.470   52.208      0.341
# rich_allInt_prim  494.726   56.304      0.000

# ...carnivores
mod_list = list(rich_simple_carn = mod_rich_simple_carn, rich_full_carn = mod_rich_full_carn, rich_allInt_carn = mod_rich_allInt_carn)
x_wts <- model_weights(mod_rich_simple_carn, mod_rich_full_carn, mod_rich_allInt_carn, weights = "kfold")
names(x_wts) = c("rich_simple_carn", "rich_full_carn", "rich_allInt_carn")
(FuncK(mod_list = mod_list, x_wts = x_wts))
# ...ungulates (CAN'T CALCULATE K-FOLD FOR ALL-INT, FOR SOME REASON)
mod_list = list(rich_simple_ung = mod_rich_simple_ung, rich_full_ung = mod_rich_full_ung, rich_allInt_ung = mod_rich_allInt_ung)
x_wts <- model_weights(mod_rich_simple_ung, mod_rich_full_ung, mod_rich_allInt_ung, weights = "kfold")
names(x_wts) = c("rich_simple_ung", "rich_full_ung", "rich_allInt_ung")
(FuncK(mod_list = mod_list, x_wts = x_wts))
# ...primates
mod_list = list(rich_simple_prim = mod_rich_simple_prim, rich_full_prim = mod_rich_full_prim, rich_allInt_prim = mod_rich_allInt_prim)
x_wts <- model_weights(mod_rich_simple_prim, mod_rich_full_prim, mod_rich_allInt_prim, weights = "kfold")
names(x_wts) = c("rich_simple_prim", "rich_full_prim", "rich_allInt_prim")
(FuncK(mod_list = mod_list, x_wts = x_wts))

# RESPONSE IS PARASITE TRANSMISSION ----
# RESULTS ----
# ...carnivores, FULL IS BARELY BETTER THAN SIMPLE
# model                   kfoldic   se(kfoldic) weight
# parastrans_full_carn    229.037   20.206      0.567
# parastrans_simple_carn  229.573   16.732      0.433
# parastrans_allInt_carn  249.698   23.663      0.000
# ...ungulates, SIMPLE IS BEST
# model                   kfoldic   se(kfoldic) weight
# parastrans_simple_ung   221.451   24.689      1
# parastrans_full_ung     237.718   24.746      0
# parastrans_allInt_ung   273.955   35.313      0
# ...primates, ALLINT IS A BIT BETTER THAN SIMPLE
# model                   kfoldic   se(kfoldic) weight
# parastrans_allInt_prim  167.325   22.707      0.703
# parastrans_simple_prim  169.048   26.259      0.297
# parastrans_full_prim    188.964    29.199     0.000

# ...carnivores
mod_list = list(parastrans_simple_carn = mod_parastrans_simple_carn, parastrans_full_carn = mod_parastrans_full_carn, parastrans_allInt_carn = mod_parastrans_allInt_carn)
x_wts <- model_weights(mod_parastrans_simple_carn, mod_parastrans_full_carn, mod_parastrans_allInt_carn, weights = "kfold")
names(x_wts) = c("parastrans_simple_carn", "parastrans_full_carn", "parastrans_allInt_carn")
(FuncK(mod_list = mod_list, x_wts = x_wts))
# ...ungulates
mod_list = list(parastrans_simple_ung = mod_parastrans_simple_ung, parastrans_full_ung = mod_parastrans_full_ung, parastrans_allInt_ung = mod_parastrans_allInt_ung)
x_wts <- model_weights(mod_parastrans_simple_ung, mod_parastrans_full_ung, mod_parastrans_allInt_ung, weights = "kfold")
names(x_wts) = c("parastrans_simple_ung", "parastrans_full_ung", "parastrans_allInt_ung")
(FuncK(mod_list = mod_list, x_wts = x_wts))
# ...primates
mod_list = list(parastrans_simple_prim = mod_parastrans_simple_prim, parastrans_full_prim = mod_parastrans_full_prim, parastrans_allInt_prim = mod_parastrans_allInt_prim)
x_wts <- model_weights(mod_parastrans_simple_prim, mod_parastrans_full_prim, mod_parastrans_allInt_prim, weights = "kfold")
names(x_wts) = c("parastrans_simple_prim", "parastrans_full_prim", "parastrans_allInt_prim")
(FuncK(mod_list = mod_list, x_wts = x_wts))

# RESPONSE IS PARASITE TYPE ----
# RESULTS ----
# ...carnivores, 'ALLINT' IS BEST BY FAR
#                   model kfoldic se(kfoldic) weight
# parastype_allInt_carn 397.796      37.077      1
# parastype_full_carn 431.285      47.007      0
# parastype_simple_carn 443.640      68.331      0
# ...ungulates, 'FULL' IS BEST BY FAR
# model kfoldic se(kfoldic) weight
# parastype_full_ung 555.892      83.852      1
# parastype_simple_ung 625.129      94.893      0
# parastype_allInt_ung 794.577     150.260      0
# ...primates, 'SIMPLE' IS BEST--ALL ARE SIMILAR
# model kfoldic se(kfoldic) weight
# parastype_simple_prim 341.349      33.540  0.999
# parastype_allInt_prim 355.293      33.426  0.001
# parastype_full_prim 364.141      34.204  0.000
#
# ...carnivores
mod_list = list(parastype_simple_carn = mod_parastype_simple_carn, parastype_full_carn = mod_parastype_full_carn, parastype_allInt_carn = mod_parastype_allInt_carn)
x_wts <- model_weights(mod_parastype_simple_carn, mod_parastype_full_carn, mod_parastype_allInt_carn, weights = "kfold")
names(x_wts) = c("parastype_simple_carn", "parastype_full_carn", "parastype_allInt_carn")
(FuncK(mod_list = mod_list, x_wts = x_wts))
# ...ungulates
mod_list = list(parastype_simple_ung = mod_parastype_simple_ung, parastype_full_ung = mod_parastype_full_ung, parastype_allInt_ung = mod_parastype_allInt_ung)
x_wts <- model_weights(mod_parastype_simple_ung, mod_parastype_full_ung, mod_parastype_allInt_ung, weights = "kfold")
names(x_wts) = c("parastype_simple_ung", "parastype_full_ung", "parastype_allInt_ung")
(FuncK(mod_list = mod_list, x_wts = x_wts))
# ...primates
mod_list = list(parastype_simple_prim = mod_parastype_simple_prim, parastype_full_prim = mod_parastype_full_prim, parastype_allInt_prim = mod_parastype_allInt_prim)
x_wts <- model_weights(mod_parastype_simple_prim, mod_parastype_full_prim, mod_parastype_allInt_prim, weights = "kfold")
names(x_wts) = c("parastype_simple_prim", "parastype_full_prim", "parastype_allInt_prim")
(FuncK(mod_list = mod_list, x_wts = x_wts))
