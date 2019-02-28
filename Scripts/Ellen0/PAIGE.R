# Please run this code. It's not like it's a complicated model, except that it takes a lot of memory for the computer to be able to calculate model fits & predictions for a zero-truncated Poisson (which is what I ended up using)

packages <- c("magrittr", "cowplot", "GGally", "scales", "tidyverse", "rstan", "brms", "broom", "tidybayes", "purrr")
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    library(x, character.only = TRUE)
  }
})

mod <-  readRDS("mod4_brm.RDS") # Change to wherever the .RDS is saved <<<<<<<<<<
modDat <- 
  

plot_me <- plot(marginal_effects(mod, effects = "c_Precip"), method = "fitted", plot = FALSE)
# I'll need this output

mu_temp <- fitted(mod)
predict_summary <- predict(mod)

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

# I'll need these outputs for a bunch of things <<<<<<<<
saveRDS(plot_me, "mod4_plot_me.RDS")
saveRDS(mu_temp, "mod4_mu_temp.RDS")
saveRDS(predict_summary, "mod4_predict_summary.RDS")
saveRDS(Resid_df, "mod4_Resid_df.RDS")