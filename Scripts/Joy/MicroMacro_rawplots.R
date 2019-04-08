#Micro vs. Macro raw data for each host clade

carn1 <- data.frame(unggTYPE$parRich_micro, unggTYPE$parRich_macro, unggTYPE$combIUCN)
carn2 <- melt(df1, id.vars='unggTYPE.combIUCN') 
#df2 <- rename('Host Clade' = 'unggTYPE.Host_Order', unggTYPE.parRich_macro = Macroparasite Richness, unggTYPE.parRich_mcro = Microparasite Richness )
head(carn2)

carn_type <-  ggplot(carn2, aes(x=unggTYPE.combIUCN, y=value, fill=variable)) +
  labs(x = 'Threat Status', y = 'Parasite Richness (uncorrected)', fill = 'Parasite Type') +
  geom_bar(stat='identity', position='dodge')

carn_type + 
  scale_fill_brewer(palette="Paired", labels  = c('Microparasite', "Macroparasite")) + 
  #scale_fill_manual(values=c("#999999", "#E69F00"), labels  = c('Non T', "T")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

prim1 <- data.frame(primTYPE$parRich_micro, primTYPE$parRich_macro, primTYPE$combIUCN)
prim2 <- melt(prim1, id.vars='primTYPE.combIUCN') 
#df2 <- rename('Host Clade' = 'primTYPE.Host_Order', primTYPE.parRich_macro = Macroparasite Richness, primTYPE.parRich_mcro = Microparasite Richness )
head(prim2)

prim_type <-  ggplot(prim2, aes(x=primTYPE.combIUCN, y=value, fill=variable)) +
  labs(x = 'Threat Status', y = 'Parasite Richness (uncorrected)', fill = 'Parasite Type') +
  geom_bar(stat='identity', position='dodge')

prim_type + 
  scale_fill_brewer(palette="Paired", labels  = c('Microparasite', "Macroparasite")) + 
  #scale_fill_manual(values=c("#999999", "#E69F00"), labels  = c('Non T', "T")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

ungg1 <- data.frame(unggTYPE$parRich_micro, unggTYPE$parRich_macro, unggTYPE$combIUCN)
ungg2 <- melt(ungg1, id.vars='unggTYPE.combIUCN') 
#df2 <- rename('Host Clade' = 'unggTYPE.Host_Order', unggTYPE.parRich_macro = Macroparasite Richness, unggTYPE.parRich_mcro = Microparasite Richness )
head(ungg2)

ungg_type <-  ggplot(ungg2, aes(x=unggTYPE.combIUCN, y=value, fill=variable)) +
  labs(x = 'Threat Status', y = 'Parasite Richness (uncorrected)', fill = 'Parasite Type') +
  geom_bar(stat='identity', position='dodge')

ungg_type + 
  scale_fill_brewer(palette="Paired", labels  = c('Microparasite', "Macroparasite")) + 
  #scale_fill_manual(values=c("#999999", "#E69F00"), labels  = c('Non T', "T")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


