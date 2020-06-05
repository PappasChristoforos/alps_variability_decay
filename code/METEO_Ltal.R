################################################################################
# METEOROLOGICAL data analyses
# 1. read the data and save the METEO data 
# 2. calulate the Decayvars and save dataframes
# 3. plot the results

# Christoforos Pappas
# last update: 03.06.2020
################################################################################

# load libraries
library(parallel)
library(reshape2)
library(openair)
library(dplyr)

library(ggplot2)
#library(ggsignif)
library(grid)
library(gtable)
library(gridExtra)
library(scales)

# functions
'%!in%' <- function(x,y)!('%in%'(x,y))
# load functions
source("D:/VAR_decay_function_fast.r")


#########
# 1. read METEO data ################################################
#########

setwd("D:/meteo/")

# Soil moisture
sm_h = read.table("Soil_moisture.txt", header=T, sep="\t") # hourly (07/14/08 14:00:00 to 11/18/15 15:00:00 )

# Air temperature
tmp = read.table("Temperature.txt", header=T, sep="\t") # 15 min (10/19/06 16:15:00 to 11/18/15 15:15:00)
tmp$date = as.POSIXct(strptime(tmp$Date, '(%m/%d/%y %H:%M:%S)'), tz = "GMT") 
tmp_h =  timeAverage(tmp, avg.time = "hour", data.thresh = 0, statistic = "mean", type = "default")

# Vapour pressure deficit 
vpd = read.table("Vapour_pressure_deficit.txt", header=T, sep="\t") # hourly (10/19/06 16:15:00 to 11/18/15 15:15:00)
vpd$date = as.POSIXct(strptime(vpd$Date, '(%m/%d/%y %H:%M:%S)'), tz = "GMT") 
vpd_h =  timeAverage(vpd, avg.time = "hour", data.thresh = 0, statistic = "mean", type = "default")


# save it ----------------------------------------------------------------------
save(list=c("sm_h", "tmp_h", "vpd_h"), 
     file="D:/postpro/METEO.RData")
# ------------------------------------------------------------------------------



#########
# 2. calculate the Decayvars ######################################################
#########


# parallel computations for DecayVAR
detectCores()
cl <- makeCluster(8)

DecayVAR_sm_h = parApply(cl, as.data.frame(sm_h[,-1]), 2, VAR_decay_fast)
DecayVAR_tmp_h = parApply(cl, as.data.frame(tmp_h[,-1]), 2, VAR_decay_fast)
DecayVAR_vpd_h = parApply(cl, as.data.frame(vpd_h[,-1]), 2, VAR_decay_fast)

        
# transform from 'list' to 'matrix'
rows_diff <- max(sapply(DecayVAR_sm_h, function(x) max(x[[1]])))
columns <- length(DecayVAR_sm_h)
# make an empty matrix
DecayVAR_sm_h_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_sm_h_matr) = 1:rows_diff
colnames(DecayVAR_sm_h_matr) = colnames(sm_h[,-1])

rows_diff <- max(sapply(DecayVAR_tmp_h, function(x) max(x[[1]])))
columns <- length(DecayVAR_tmp_h)
# make an empty matrix
DecayVAR_tmp_h_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_tmp_h_matr) = 1:rows_diff
colnames(DecayVAR_tmp_h_matr) = colnames(tmp_h[,-1])

rows_diff <- max(sapply(DecayVAR_vpd_h, function(x) max(x[[1]])))
columns <- length(DecayVAR_vpd_h)
# make an empty matrix
DecayVAR_vpd_h_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_vpd_h_matr) = 1:rows_diff
colnames(DecayVAR_vpd_h_matr) = colnames(vpd_h[,-1])


# fill the matrix
for (i in colnames(DecayVAR_sm_h_matr)){# loop through trees
  DecayVAR_sm_h_matr[DecayVAR_sm_h[[i]]$d,i] <- DecayVAR_sm_h[[i]]$v
}# loop through trees

for (i in colnames(DecayVAR_tmp_h_matr)){# loop through trees
  DecayVAR_tmp_h_matr[DecayVAR_tmp_h[[i]]$d,i] <- DecayVAR_tmp_h[[i]]$v
}# loop through trees

for (i in colnames(DecayVAR_vpd_h_matr)){# loop through trees
  DecayVAR_vpd_h_matr[DecayVAR_vpd_h[[i]]$d,i] <- DecayVAR_vpd_h[[i]]$v
}# loop through trees



x11()
matplot(x=log10(1:nrow(DecayVAR_sm_h_matr)), y=log10(DecayVAR_sm_h_matr), type="l")
x11()
matplot(x=log10(1:nrow(DecayVAR_tmp_h_matr)), y=log10(DecayVAR_tmp_h_matr), type="l")
x11()
matplot(x=log10(1:nrow(DecayVAR_vpd_h_matr)), y=log10(DecayVAR_vpd_h_matr), type="l")



# calculate standard deviation
Decaystdev_sm_h_matr = sqrt(DecayVAR_sm_h_matr)
Decaystdev_tmp_h_matr = sqrt(DecayVAR_tmp_h_matr)
Decaystdev_vpd_h_matr = sqrt(DecayVAR_vpd_h_matr)

# save it ----------------------------------------------------------------------
save(list=c("Decaystdev_sm_h_matr",
            "Decaystdev_tmp_h_matr",
            "Decaystdev_vpd_h_matr"),
     file="D:/postpro/METEO_Decaystdev.RData")
# ------------------------------------------------------------------------------



# organise data ----------------------------------------------------------------
sm10_sc = as.data.frame(Decaystdev_sm_h_matr[,c("N08_10", "N13_10", "N16_10", "N19_10", "N22_10")])
sm70_sc = as.data.frame(Decaystdev_sm_h_matr[,c("N08_70", "N13_70", "N16_70", "N19_70", "N22_70")])
tmp_sc = as.data.frame(Decaystdev_tmp_h_matr[,c("N08", "N13", "N16", "N19", "N22")])
vpd_sc = as.data.frame(Decaystdev_vpd_h_matr[,c("N08", "N13", "N16", "N19", "N22")])


sm10_sc$Hours = as.numeric(rownames(sm10_sc))
gg_sm10_sc = melt(sm10_sc, id='Hours')

lm_sm10_sc = subset(gg_sm10_sc, Hours >= 24) %>% 
    group_by(variable) %>% 
    do({
      mod = lm(log10(value) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_sm10_sc$Type  = "sm10"

lm_sm10_sc$Elevation = NA
lm_sm10_sc[substr(lm_sm10_sc$variable, 2,3) == "08", "Elevation"] = 800
lm_sm10_sc[substr(lm_sm10_sc$variable, 2,3) == "13", "Elevation"] = 1300
lm_sm10_sc[substr(lm_sm10_sc$variable, 2,3) == "16", "Elevation"] = 1600
lm_sm10_sc[substr(lm_sm10_sc$variable, 2,3) == "19", "Elevation"] = 1900
lm_sm10_sc[substr(lm_sm10_sc$variable, 2,3) == "22", "Elevation"] = 2200


sm70_sc$Hours = as.numeric(rownames(sm70_sc))
gg_sm70_sc = melt(sm70_sc, id='Hours')

lm_sm70_sc = subset(gg_sm70_sc, Hours >= 24) %>% 
    group_by(variable) %>% 
    do({
      mod = lm(log10(value) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_sm70_sc$Type  = "sm70"

lm_sm70_sc$Elevation = NA
lm_sm70_sc[substr(lm_sm70_sc$variable, 2,3) == "08", "Elevation"] = 800
lm_sm70_sc[substr(lm_sm70_sc$variable, 2,3) == "13", "Elevation"] = 1300
lm_sm70_sc[substr(lm_sm70_sc$variable, 2,3) == "16", "Elevation"] = 1600
lm_sm70_sc[substr(lm_sm70_sc$variable, 2,3) == "19", "Elevation"] = 1900
lm_sm70_sc[substr(lm_sm70_sc$variable, 2,3) == "22", "Elevation"] = 2200


tmp_sc$Hours = as.numeric(rownames(tmp_sc))
gg_tmp_sc = melt(tmp_sc, id='Hours')

lm_tmp_sc = subset(gg_tmp_sc, Hours >= 24 & Hours<=24*30*2 ) %>% 
    group_by(variable) %>% 
    do({
      mod = lm(log10(value) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_tmp_sc$Type  = "tmp"

lm_tmp_sc$Elevation = NA
lm_tmp_sc[substr(lm_tmp_sc$variable, 2,3) == "08", "Elevation"] = 800
lm_tmp_sc[substr(lm_tmp_sc$variable, 2,3) == "13", "Elevation"] = 1300
lm_tmp_sc[substr(lm_tmp_sc$variable, 2,3) == "16", "Elevation"] = 1600
lm_tmp_sc[substr(lm_tmp_sc$variable, 2,3) == "19", "Elevation"] = 1900
lm_tmp_sc[substr(lm_tmp_sc$variable, 2,3) == "22", "Elevation"] = 2200


vpd_sc$Hours = as.numeric(rownames(vpd_sc))
gg_vpd_sc = melt(vpd_sc, id='Hours')

lm_vpd_sc = subset(gg_vpd_sc, Hours >= 24) %>% 
    group_by(variable) %>% 
    do({
      mod = lm(log10(value) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_vpd_sc$Type  = "vpd"

lm_vpd_sc$Elevation = NA
lm_vpd_sc[substr(lm_vpd_sc$variable, 2,3) == "08", "Elevation"] = 800
lm_vpd_sc[substr(lm_vpd_sc$variable, 2,3) == "13", "Elevation"] = 1300
lm_vpd_sc[substr(lm_vpd_sc$variable, 2,3) == "16", "Elevation"] = 1600
lm_vpd_sc[substr(lm_vpd_sc$variable, 2,3) == "19", "Elevation"] = 1900
lm_vpd_sc[substr(lm_vpd_sc$variable, 2,3) == "22", "Elevation"] = 2200



# Put all variables together ---------------------------------------------------
slopes_df = rbind(data.frame(lm_sm10_sc)[,-1], 
                  data.frame(lm_sm70_sc)[,-1],
                  data.frame(lm_tmp_sc)[,-1],
                  data.frame(lm_vpd_sc)[,-1])

slopes_df = subset(slopes_df, Type %in% c("sm10", "tmp", "vpd"))

slopes_df$Type[slopes_df$Type == "vpd"] = "c. D"
slopes_df$Type[slopes_df$Type == "sm10"] = "b. S"
slopes_df$Type[slopes_df$Type == "tmp"] = "a. T"

slopes_df$Type = factor(slopes_df$Type)




bx_meteo = ggplot(data= slopes_df, aes(x=Type, y=Slope)) + 
  geom_boxplot(aes(fill=Type), outlier.shape = NA, show.legend = FALSE) +   

#  geom_signif(comparisons = list(c("Larix decidua", "Picea abies")),
#  test="t.test", test.args=list(alternative = "less"), map_signif_level=T)+
  scale_y_continuous("Slope [-]", limits=c(-0.5, 0))+
  xlab("")+
  theme_bw()+
             theme(aspect.ratio = 1.5,
                   panel.grid.minor.x = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   strip.background = element_blank())

ggsave("D:/Figs/METEO_slopes_bx.eps", plot=bx_meteo, device="eps", height=8, units="cm")

grad_meteo = ggplot(data=slopes_df, aes(group=Type, x=Elevation, y=Slope, colour=Type, fill=Type))+
  
  geom_hline(yintercept = c(-0.5, 0), linetype = "dashed", colour="gray80")+ 
  geom_smooth( method="loess", se=T, show.legend = FALSE, alpha=1)+
  geom_point(shape=21, colour="white", show.legend = FALSE, size=2)+
                                                                         
#  geom_tile(mapping=aes(group=Type, x=Elevation, y=Slope, fill=Type))+

              
  scale_x_continuous("Elevation [m asl.]", 
                     limits=c(0,2500), 
                     breaks=c(0, 800, 1300, 1600, 1900, 2200))+
             
  scale_y_continuous("Slope [-]", 
                     limits=c(-0.5,0))+
                     
  theme_bw()+
  
  theme(aspect.ratio = 1.4,
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_blank())  +
  
  
  coord_flip()
  
ggsave("D:/Figs/METEO_slopes_gradient.eps", plot=grad_meteo, device="eps", width=8, units="cm")



## Summaries
slopes_df %>%
  group_by(Type) %>%
  summarize(mean_slpe = mean(Slope, na.rm=T),
            sd_slpe = sd(Slope, na.rm=T))
            