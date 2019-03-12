# visualize marginal effects
library(cowplot)

#### --- SIMPLE models of parasite richness on simple dataset --- #### 

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
pdf("./Results/richness_me_plot.pdf") 
plot_grid(allp, bot, nrow=2, labels=c("A"))
dev.off()

#### --- SIMPLE models of parasite proportions on simple dataset --- #### 

simple_me_parastrans <- readRDS("./Data/JPEK/simple/simple_brm_parastrans_me.RDS")
plot(simple_me_parastrans[[4]])

simple_me_parastype <- readRDS("./Data/JPEK/simple/simple_brm_parastype_me.RDS")
plot(simple_me_parastype[[4]])

simpTrans <- simple_me_parastrans[[4]] + ylab("Odds close") + theme(legend.position = "none") 
simpType <- simple_me_parastype[[4]] + ylab("Odds micro") + theme(legend.position = "none") 


botProp <- plot_grid(plotlist = list(simpTrans, simpType), 
                     #labels = c("B", "C", "D", "E"), 
                     ncol=1, nrow=2)

pdf("./Results/propSimple_me_plot.pdf") 
plot_grid(allp, botProp, nrow=2
          #, labels=c("A")
          )
dev.off()

#### --- FULL models of parasite proportions on FULL dataset --- #### 

full_me <- readRDS("./Data/JPEK/full/full_brm_me.RDS")
plot(full_me[[8]])

full_me_parastrans <- readRDS("./Data/JPEK/full/full_brm_parastrans_me.RDS")
plot(full_me_parastrans[[8]])

full_me_parastype <- readRDS("./Data/JPEK/full/full_brm_parastype_me.RDS")
plot(full_me_parastype[[8]])

allp <- full_me[[8]] + ylab("Richness")
fullTrans <- full_me_parastrans[[8]] + ylab("Odds close") + theme(legend.position = "none") 
fullType <- full_me_parastype[[8]] + ylab("Odds micro") + theme(legend.position = "none") 

botFull <- plot_grid(plotlist = list(simpTrans, simpType), 
                     #labels = c("B", "C", "D", "E"), 
                     ncol=1, nrow=2)

pdf("./Results/propFull_plot.pdf") 
plot_grid(allp, botProp, nrow=2
          #, labels=c("A")
)
dev.off()


#### --- SIMPLE models of parasite proportions on FULL dataset --- #### 

simple_me_f <- readRDS("./Data/JPEK/simple/simple_brm_me_fulldat.RDS")
plot(simple_me_f[[4]])

simple_me_f_parastrans <- readRDS("./Data/JPEK/simple/simple_brm_parastrans_me_fulldat.RDS")
plot(simple_me_f_parastrans[[4]])

simple_me_f_parastype <- readRDS("./Data/JPEK/simple/simple_brm_parastype_me_fulldat.RDS")
plot(simple_me_f_parastype[[4]])

allp <- simple_me_f[[4]] + ylab("Richness")
simpleTransFull <- simple_me_f_parastrans[[4]] + ylab("Odds close") + theme(legend.position = "none") 
simpleTypeFull <- simple_me_f_parastype[[4]] + ylab("Odds micro") + theme(legend.position = "none") 

botSimple_fulldat <- plot_grid(plotlist = list(simpleTransFull, simpleTypeFull), 
                     #labels = c("B", "C", "D", "E"), 
                     ncol=1, nrow=2)

pdf("./Results/propsimple_fulldat_plot.pdf") 
plot_grid(allp, botSimple_fulldat, nrow=2
          #, labels=c("A")
)
dev.off()