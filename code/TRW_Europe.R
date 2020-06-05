################################################################################
# TREE RING WIDTH data analyses across Europe (ITRDB)
# Christoforos Pappas
# last update: 03.06.2020
################################################################################


####################
# Calculate slopes #------------------------------------------------------------
####################

# load libraries
library(reshape2)
library(plyr)
library(dplyr)

# clean memory
rm(list=ls())

# set path
path = "D:/TRW_EU/"

# load functions
source("D:/VAR_decay_function_fast.r")
`%!in%` = Negate(`%in%`)

setwd(path)

trw_PCAB = read.table("PCAB_endyr_1980_lser_150_coff_50_zercu_6_min_3.txt", header=T)    
trw_PCAB$ID = paste(trw_PCAB$site.id, trw_PCAB$tree_nr, sep="-")   
trw_PCAB_sf = dcast(trw_PCAB,year+sp.code ~ ID, value.var = "rw_mm")
    
test1 = apply(trw_PCAB_sf[,-c(1,2)], 2, VAR_decay_fast)
test2 = ldply(test1, rbind)  
test2$v = sqrt(test2$v)

slopes1 = test2 %>% 
 group_by(.id) %>% 
 do({
  mod = lm(log10(v) ~ log10(d), data = .)
  data.frame(Slope = coef(mod)[2],
  # StdError = summary(mod)$coefficients[2,2],      
  R2adj = summary(mod)$adj.r.squared)
  })
slopes1$Species = "Picea abies"

slopes_PCAB = slopes1 
colnames(slopes_PCAB) = c("siteIDtreeNR", "Slope", "R2adj", "Species")


trw_LADE = read.table("LADE_endyr_1980_lser_150_coff_50_zercu_6_min_3.txt", header=T)    
trw_LADE$ID = paste(trw_LADE$site.id, trw_LADE$tree_nr, sep="-")   
trw_LADE_sf = dcast(trw_LADE,year+sp.code ~ ID, value.var = "rw_mm")
    
test1 = apply(trw_LADE_sf[,-c(1,2)], 2, VAR_decay_fast)
test2 = ldply(test1, rbind)  
test2$v = sqrt(test2$v)

slopes1 = test2 %>% 
 group_by(.id) %>% 
 do({
  mod = lm(log10(v) ~ log10(d), data = .)
  data.frame(Slope = coef(mod)[2],
  # StdError = summary(mod)$coefficients[2,2],      
  R2adj = summary(mod)$adj.r.squared)
  })
slopes1$Species = "Larix decidua"

slopes_LADE = slopes1 
colnames(slopes_LADE) = c("siteIDtreeNR", "Slope", "R2adj", "Species")


# load metadata of LADE and PCAB in Europe
metadata = read.table("D:/TRW_EU/metadata.txt", header=T)

slopes_df = rbind(data.frame(slopes_PCAB), data.frame(slopes_LADE))
l = strsplit(slopes_df$siteIDtreeNR,"-")
df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
colnames(df) = c("siteID", "treeNR")
slopes_df = cbind(slopes_df, df)
slopes_df$id = as.numeric(as.character(slopes_df$siteID))

# merge slopes and metadata
slopes_df = merge(metadata, slopes_df, by="id",  all.x=T)
slopes_df = subset(slopes_df, species_short %in% c("LADE", "PCAB")) 

TRW_slopes_EU = slopes_df
save(TRW_slopes_EU, file = "D:/postpro/TRW_slopes_EU.RData")
