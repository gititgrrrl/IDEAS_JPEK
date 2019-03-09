# POINT 1:
# Many parasites are classified as both close & non-close. What we are actually interested in determining is how the relative abundance of parasites that REQUIRE close contact for transmission is influenced by threat status, etc. 
# Definitions:
# Close = Parasite transmissible by close contact (HYPOTHESIZE LOW TRANSMISSION AND IT DEPENDS ON GROUP SIZE)
# NonClose = Parasite transmissible by non-close means (HYPOTHESIZE EASY TRANSMISSION)
# Intermediate = Complex life-cycle (HYPOTHESIZE DEPENDS ON TYPE OF THREAT)
# Vector = Vectors such as arthropods, etc. (HYPOTHESIZE EASY TRANSMISSION)
# We find the largest overlap of nonclose+close and of nonclose+intermediate; CloseOnly = 230
# >>>>> I changed it so that CloseOnly = parasites that can ONLY be transmitted by direct contact, and NonClose = parasites that can be transmitted in some way that doesn't require direct <<<<<

# POINT 2:
# Once we add in the parasites coded only to genus, then 31% of the records don't have information on parasite transmission (whether it's direct transmission only or not). That means that when we calculated nonclose (for model) as parasite rich-close, then we counted all those 31% as nonclose. So actually in that model, should have the CloseOnly as 0, 1, NA and remove the NA's (but note in methods that we are only using a subset of the data) 

# POINT 3:
# Carnivore group categories are:
# > unique(carnDat$Groupsize)
# [1] "Solitary"              "Pairs"                
# [3] "Groups"                NA                     
# [5] "Solitary or Pairs"     "Solitary?"            
# [7] "Pairs to family clans" "Solitary or Groups"   
# [9] "Solitary to Pairs"     "Pairs?"               
# [11] "Pairs to Groups"  
# In the script, these were counted as groups (everything else non-group): "Groups", "Pairs to family clans", "Solitary or groups". And I think the next line was... carnDat$carnGrp[carnDat$carnGrp!="Group"] <- "non_group" (<<< CHECK because it should have been lowercase !="group")
# >>>>> I changed it so that these were counted as groups: "Groups", "Solitary or Groups" (NOTE the capitalization), "Pairs to family clans", "Pairs to Groups" <<<<<

# POINT 4:
# The simple model should be: bf(ParRich | trunc(lb = 1) ~ HostIUCNcomb*HostGroup + logNumHostCitations*HostIUCNcomb <<<<<< ADD THE INTERACTION AND LOG-TRANS. CITATIONS

# POINT 5: 
# > The simple subset models should be logistic. Otherwise, we can't really say anything about how threat status affects the relative abundance of subtypes--we only look at plots. So for example, we look at ungulate plots. Non-threatened ungulates have much higher any parasite richness and we see in the subsets that non-threatened still have higher. We can't tell from these plots if the relative % changes b/c anyway the overall richness is so much higher for non-threatened. [THIS IS WHAT ANDREW PARK SUGGESTED]
# > Also, use the new definitions of CloseOnly and NonClose
# > For parasite type, since we are just splitting into micro- and macroparasites, then we should include the prion data (as microparasites)--also, should not have excluded fungal parasites from the start anyway b/c they should contribute to total parasite richness analysis

# POINT 6:
# > Next step is full models. For each of the three analyses, we want to add in all the predictors in various combinations and then use model selection. Except remember that we should remove %close transmit from the predictors.
# > Sonya suggested using absolute latitude, so I changed that in the script
# > Remember that if we are going to do model selection comparing simple with more complex models, we have to use the same dataset. So the ones we are doing right now can be used to see results with the more complete dataset, but for model selection we have to rerun simple models on the subset of data that is used for the full models also.
# > Are we going to try to do something with social group size? If so, comGrpSize seems odd. It's the mean of social & population group size, but seems like we should use the max of those (>>>> I created maxGrpSize. In final data it's 'groupSizePriUng_max' <<<<). Also note that there are some social group sizes that are larger than population group size, which seems wonky

# -----
# MISC:
# Summary of simpleDat (missing IUCN = no_data)
#                  carnivores primates ungulates
# no_data                 7       11        12
# not_threatened         98       61        61
# threatened             43       70        29
# NOTE: The IUCN no_data tended to have low parasite richness (a bit lower than threatened)
# 
# FOR ACTUAL PAPER:
# > Might make more sense to compare results of analyses in which all parasites not identified to species are excluded vs. the opposite extreme
