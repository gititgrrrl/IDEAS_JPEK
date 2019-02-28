# fits and predict values for mod 4 & 7

library(magrittr)
library(tidyverse)
library(brms)
library(tidybayes)

modDat <- readRDS("ModDat.RDS")

mod4 <-  readRDS("mod4_brm.RDS") 
mod7 <-  readRDS("mod7_brm.RDS") 

# plots of marginal effects of predictors

# throws an error for mod4
plot_4 <- plot(marginal_effects(mod4), method = "fitted", plot = FALSE)

# also throws an error for mod7 
plot_7 <- plot(marginal_effects(mod7), method = "fitted", plot = FALSE)


          # fit and predict model 4
mu_4 <- fitted(mod4)
predict_4 <- predict(mod4)


# fit and predict model 7
mu_7 <- fitted(mod7)
predict_7 <- predict(mod7)

# residuals model 4
Resid_4 <-cbind(modDat, enframe(residuals(mod4, type = "pearson")[, "Estimate"])) %>%
  rename(PResid = value) %>%
  select(-name) %>%
  cbind(enframe(fitted(mod)[, "Estimate"])) %>%
  rename(Fitted = value) %>%
  select(-name) %>%
  cbind(enframe(predict(mod)[, "Estimate"])) %>%
  rename(Predicted = value) %>%
  select(-name) #%>%
  # left_join(FinalDat[,c("HostName", "HostActivCycle", "HostHomeRange", "HostTrophic", "HostMeanLat")], by = "HostName") %>% # add predictors NOT used in the model
  # mutate_if(is.character, as.factor) %>%
  # mutate(logHostHomeRange = log(HostHomeRange))


# residuals model 7
Resid_7 <-cbind(modDat, enframe(residuals(mod7, type = "pearson")[, "Estimate"])) %>%
  rename(PResid = value) %>%
  select(-name) %>%
  cbind(enframe(fitted(mod)[, "Estimate"])) %>%
  rename(Fitted = value) %>%
  select(-name) %>%
  cbind(enframe(predict(mod)[, "Estimate"])) %>%
  rename(Predicted = value) %>%
  select(-name) #%>%
  # left_join(FinalDat[,c("HostName", "HostActivCycle", "HostHomeRange", "HostTrophic", "HostMeanLat")], by = "HostName") %>% # add predictors NOT used in the model
  # mutate_if(is.character, as.factor) %>%
  # mutate(logHostHomeRange = log(HostHomeRange))


# Marginal effects plots  
saveRDS(plot_4, "mod4_plot_me.RDS")
saveRDS(plot_7, "mod7_plot_me.RDS")

# model 4
saveRDS(mu_4, "mod4_mu_temp.RDS")
saveRDS(predict_4, "mod4_predict_summary.RDS")

# model 7
saveRDS(mu_7, "mod7_mu_temp.RDS")
saveRDS(predict_7, "mod7_predict_summary.RDS")

# Residuals : error:  Error: object of type 'builtin' is not subsettable
# saveRDS(Resid_df, "mod4_Resid_df.RDS")
# saveRDS(Resid_4, "mod4_Resid_df.RDS")
# saveRDS(Resid_7, "mod7_Resid_df.RDS")

