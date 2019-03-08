# visualize marginal effects
library(cowplot)

simple_me <- readRDS("./Data/JPEK/simple/simple_brm_me.RDS")
plot(simple_me[[4]])

simple_me_close <- readRDS("./Data/JPEK/simple/simple_brm_me_close.RDS")
plot(simple_me_close[[4]])

simple_me_nonclose <- readRDS("./Data/JPEK/simple/simple_brm_me_nonclose.RDS")
plot(simple_me_nonclose[[4]])

simple_me_micro <- readRDS("./Data/JPEK/simple/simple_brm_me_micro.RDS")
plot(simple_me_micro[[4]])

simple_me_macro <- readRDS("./Data/JPEK/simple/simple_brm_me_macro.RDS")
plot(simple_me_macro[[4]])

allp <- simple_me[[4]] + ylab("Richness")
close <- simple_me_close[[4]] + ylab("Richness") + theme(legend.position = "none") 
nonclose <- simple_me_nonclose[[4]] + ylab("Richness") + theme(legend.position = "none") 
micro <- simple_me_micro[[4]] + ylab("Richness") + theme(legend.position = "none") 
macro <- simple_me_macro[[4]] + ylab("Richness") + theme(legend.position = "none") 


top <- allp
bot <- plot_grid(plotlist = list(close, nonclose, micro, macro), 
          labels = c("B", "C", "D", "E"), 
           ncol=2, nrow=2)

plot_grid(allp, bot, nrow=2, labels=c("A"))
