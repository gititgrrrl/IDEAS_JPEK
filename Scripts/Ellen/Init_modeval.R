### CHECK THE MODEL ----
mod <- mod7_brm

plot(mod)
stanplot(mod) # quick coefficient plot

# Posterior predictive checks ----
pp_check(mod, type = "dens_overlay", nsamples = 300)

for(g in c("carnivores", "ungulates", "primates")) {
  plot(pp_check(mod, type = "dens_overlay", nsamples = 300, newdata = subset(ModDat, HostGroup == g)) +
    ggtitle(paste0(g, " density plot")) +
    theme(plot.title = element_text(size = 12, face = "bold")))
}

for(s in c("threatened", "not_threatened")) {
  plot(pp_check(mod, type = "dens_overlay", nsamples = 300, newdata = subset(ModDat, HostIUCNcomb == s)) +
         ggtitle(paste0(s, " density plot")) +
         theme(plot.title = element_text(size = 12, face = "bold")))
}

pp_check(mod, type = "violin_grouped", nsamples = 300, group = "HostGroup")
pp_check(mod, type = "violin_grouped", nsamples = 300, group = "HostIUCNcomb")
pp_check(mod, type = "loo_pit_qq", nsamples = 300) 

# Marginal effects ----
plot_me <- marginal_effects(mod, method = "fitted", plot = FALSE)

# Residual plots ----
Resid_df <-cbind(ModDat, enframe(residuals(mod, type = "pearson")[, "Estimate"])) %>%
  rename(PResid = value) %>%
  select(-name) %>%
  cbind(enframe(fitted(mod)[, "Estimate"])) %>%
  rename(Fitted = value) %>%
  select(-name) %>%
  cbind(enframe(predict(mod)[, "Estimate"])) %>%
  rename(Predicted = value) %>%
  select(-name) %>%
  left_join(FinalDat[,c("HostName", "HostActivCycle", "HostHomeRange", "HostTrophic", "HostMeanLat")], by = "HostName") %>% # add predictors NOT used in the model
  mutate_if(is.character, as.factor) %>%
  mutate(logHostHomeRange = log(HostHomeRange))

# Residual (boxplots) by IUCN and group
FuncPlotBox(dat_df = Resid_df, x_nam = "HostGroup", y_nam = "PResid", trans_y = NULL, facet_nam = "HostIUCNcomb", cens = FALSE, plot_labs = list(x = "Host Group", y = "Pearson Residual", title = "Residual Plots"))

# Residuals against continuous variables
for (x_var in c("NumHostCitations", "logHostSpeciesRange", "logHostMass", "HostMaxLifespan", "propParCloseTrans")) {
  plot(FuncPlotScatter(dat_df = Resid_df, x_nam = x_var, y_nam = "PResid", color_by = "HostGroup", trans_y = NULL, add_loess = FALSE, cens = FALSE, plot_labs = list(x = x_var, y = "Pearson Residual", title = "Residual Plots")))
}

# Residuals against variables not used in full model
# > ... categorical
FuncPlotBox(dat_df = Resid_df, x_nam = "HostActivCycle", y_nam = "PResid", trans_y = NULL, facet_nam = "HostTrophic", cens = FALSE, plot_labs = list(x = "Host Activity Cycle", y = "Pearson Residual", title = "Residual Plots"))
# > ... continuous
for (x_var in c("logHostHomeRange", "HostMeanLat")) {
  plot(FuncPlotScatter(dat_df = Resid_df, x_nam = x_var, y_nam = "PResid", color_by = "HostGroup", trans_y = NULL, add_loess = FALSE, cens = FALSE, plot_labs = list(x = x_var, y = "Pearson Residual", title = "Residual Plots")))
}

# Coefficient plot ----
FuncPlotCoef(mod_list, terms_vec = NULL, exclude_terms = FALSE)

# Model comparison ----
x <- compare_ic(mod2_brm, mod3_brm, mod4_brm, mod5_brm, mod7_brm, ic= "kfold")
x_wts <- model_weights(mod2_brm, mod3_brm, mod4_brm, mod5_brm, mod7_brm, weights = "kfold")
names(x_wts) = c("mod2", "mod3", "mod4", "mod5", "mod7")
mod_list <- list(mod2 = mod2_brm, mod3 = mod3_brm, mod4 = mod4_brm, mod5 = mod5_brm, mod7 = mod7_brm)

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
saveRDS(Kfold_summary, "kfold_table.RDS")

### SUMMARIES TO PRESENT ----
mod_list <- list(mod2 = mod2_brm, mod3 = mod3_brm, mod4 = mod4_brm, mod5 = mod5_brm)

# KFOLD table
Kfold_summary <- readRDS("kfold_table.RDS")
  
# Model conditional predictions -- DONE 
# plot_me <- plot(marginal_effects(mod), method = "fitted", plot = FALSE)
# plot_me <- readRDS("~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen0/mod7_plot_me.RDS")

pdf("plot_CE_SpeciesRange.pdf")
plot_me$`logHostSpeciesRange` +
  labs(title = "Conditional effects: Species range", x = "ln(Species range)", y = "Parasite richness") +
  # scale_x_continuous(trans='log', breaks= pretty_breaks()) +
  # scale_y_continuous(trans='log', breaks= pretty_breaks()) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw(base_size = 11)
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



# Posterior predictive fit
dat_df = ModDat; mod_list = list(mod7 = mod7_brm); resp_colnam = "ParRich"; obsID_colnam = "HostName"; color_colnam = "HostIUCNcomb"; coloraxt_title = "Host IUCN status"; facet_colnam = "HostGroup"

mu_temp <- readRDS("~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen0/mod7_mu_temp.RDS")
predict_summary <- readRDS("~/Google Drive/GMPD/GlobalParasites/Scripts/Ellen0/mod7_predict_summary.RDS")

  names(dat_df)[names(dat_df) == resp_colnam] <- "resp"
  names(dat_df)[names(dat_df) == obsID_colnam] <- "group_var"
  if(!is.null(color_colnam)) {
    names(dat_df)[names(dat_df) == color_colnam] <- "color_var"
  }
  if(!is.null(facet_colnam)) {
    names(dat_df)[names(dat_df) == facet_colnam] <- "facet_var"
  }
  
  summary_list <- p1_list <- p2_list <- vector("list", length = length(mod_list))
  
  for (i in 1:length(mod_list)) {
    mod <- mod_list[[i]]
    
    # Model fits w/95% CI
    mu_temp %<>%
      as_tibble() %>%
      dplyr::rename(mu_Q2.5 = Q2.5,
                    mu_Q97.5 = Q97.5) %>%
      dplyr::select(-Est.Error)
    mu_temp <- cbind(dat_df, mu_temp)
    
    # Model predicted responses w/95% CI
    predict_summary %<>%
      as_tibble() %>%
      transmute(pred_Q2.5 = Q2.5,
                pred_Q97.5 = Q97.5)
    
    summary_list[[i]] <- cbind(mu_temp, predict_summary) %>%
      dplyr::mutate(mu_err = Estimate - resp,
                    mu_err_Q2.5 = mu_Q2.5 - resp,
                    mu_err_Q97.5 = mu_Q97.5 - resp) %>%
      dplyr::mutate(group_var = factor(group_var, levels = group_var[order(mu_err)]))
  } # end of mod_list
  
  axes_min = floor(min(sapply(summary_list, function(x) min(subset(x, select = c(mu_Q2.5, pred_Q2.5, resp)), na.rm=TRUE))))
  axes_max = ceiling(max(sapply(summary_list, function(x) max(subset(x, select = c(mu_Q97.5, pred_Q97.5, resp)), na.rm=TRUE))))
  
  for (i in 1:length(mod_list)) {
    plot_dat <- summary_list[[i]]
    # Plot 1
    p1 <- ggplot(data = plot_dat, aes(x = resp, y = Estimate)) + 
      # lims(x = c(axes_min, axes_max), y = c(axes_min, axes_max)) +
      geom_abline(intercept = 0, slope = 1, 
                  linetype = "dashed", color = "gray") +
      labs(x =  "Observed parasite richness", y = "Predicted parasite richness", subtitle = names(mod_list[i])) +
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
    
  
  # combine p1's
  title_p1 <- ggdraw() +
    draw_label("Predicted vs. observed", fontface = "bold", size = 12)
  subtitle_p1 <- ggdraw() +
    draw_label("Vertical lines are 95% confidence intervals of the average prediction.\nThe gray dashed line shows perfect prediction.", size = 10)
  
  final_p1 <- plot_grid(title_p1, subtitle_p1, plotlist= p1_list, ncol = 1, rel_heights = c(0.07, 0.13, rep(1, length(p1_list))))
  
  }
  
  pdf("plot_predict_obs.pdf")
  final_p1
  dev.off()
  