################################################################################
# DENDROMETER data analyses
# 1. read the data and save the DENDROMETER data (including hourly data, hourly differences, and growth-, hydro-components)in RData format
# 2. calulate the Decayvars and save dataframes

# Christoforos Pappas
# last update: 03.06.2020
################################################################################

# load libraries
library(parallel)
library(reshape2)

# functions
'%!in%' <- function(x,y)!('%in%'(x,y))

#########
# 1. read and save dendrometer data ############################################
#########

#  working directory
setwd("D:/DEN/")

# read the data file; this refers to stem RADIUS in [um]; data are clipped to the growing season
DEN_data <- read.table("clipped.dendrometer.txt", header=T, sep="\t")

dates = as.character(DEN_data$Date)
dates = substr(dates, 2, 18)           

timestep <- data.frame(timestamp = as.POSIXct(strptime(dates, '%m/%d/%y %H:%M:%S'), tz = "GMT"))

rownames(DEN_data) = timestep$timestamp 
DEN_data=DEN_data[,-1]

# remove the data from the wet site and also "N08Ad_L5p" and "S19Ad_L5p" cause they have very few years of data and N08Ad_L1p because of low quality data
#DEN_selection = DEN_data[,colnames(DEN_data)%!in% c("N13Wd_L1p", "N13Wd_S1p", "N13Wd_S2p", "N13Wd_L2p", "N13Wd_L3p", "N13Wd_S3p", "N08Ad_L5p")]
 
which_dendro = c("N08Ad_L4p", "N08Bd_S3p", "N08Bd_S4p", "N13Ad_S1p", "N13Ad_S2p", "N13Bd_L1p", "N13Bd_L2p", "N16Ad_L1p", "N16Ad_L2p", "N16Bd_S1p", "N16Bd_S2p", "N19Ad_L1p", "N19Ad_L2p", "N19Bd_S1p", "N19Bd_S2p", "N22Ad_L1p", "N22Ad_L2p", "S16Ad_L1p", "S16Ad_S2p", "S16Bd_L1p", "S16Bd_S2p", "S19Ad_L1p", "S19Ad_S2p", "S19Bd_L1p", "S19Bd_S2p", "S22Ad_L1p", "S22Ad_L2p")

DEN_selection = DEN_data[,colnames(DEN_data)%in% which_dendro]

# dates
years = substr(rownames(DEN_selection), 1, 4)
months = substr(rownames(DEN_selection), 6, 7)
days = substr(rownames(DEN_selection), 9, 10)
hours = substr(rownames(DEN_selection), 12, 13)

##
# Hourly data ------------------------------------------------------------------
##
DEN_h_df = data.frame(Year= years, 
                    Month= months,
                    Day= days,
                    Hour=hours,
                    DEN_selection)
                    
##                    
# Hourly Differences -----------------------------------------------------------
##
DEN_h_df_diff = apply(DEN_selection, 2, diff, lag = 1, differences = 1)

            
##
# Zweifel et al 2016 NP --------------------------------------------------------
##
DEN_selection_z = DEN_selection

for (jj in 1:ncol(DEN_selection)){

  w = matrix(DEN_selection[,jj])
  colnames(w) = colnames(DEN_selection)[jj]
  rownames(w) = rownames(DEN_selection)
  
  w_cc = w[complete.cases(w[,1]),]
  w_z = as.matrix(w_cc)
  colnames(w_z) = colnames(w)
    
    for (ii in 2:nrow(w_z)){
  
      if(w_z[ii]<=w_z[ii-1]){
      w_z[ii,1] = w_z[ii-1,1]}
    
    }
  
  w_f = data.frame(dates=rownames(w))
  
  w_z = data.frame(w_z)
  w_z$dates = rownames(w_z)
  
  test = merge(w_f, w_z, by="dates", all=T)
  
  DEN_selection_z[,jj] = test[,2]

}

GRO_zwf = DEN_selection_z
             
             
# growth differences
GRO_zwf_diff = apply(GRO_zwf, 2, diff, lag = 1, differences = 1)



##
# King et al 2013 AFM ----------------------------------------------------------
##
DEN_d_df=c()
DEN_dh_df=c()
DR_hD_df = c()

for (i in 1:ncol(DEN_selection)){# loop through trees
  
  Den_h = cbind(DEN_h_df[,1:4],Den_h=DEN_selection[,i])
  
  Den_d = aggregate(Den_h$Den_h, by=list(Y=years, M=months, D=days), mean, na.rm=T)
  colnames(Den_d) = c("Year", "Month", "Day", "Den_d")

  Den_d = Den_d[order(Den_d[,1],Den_d[,2],Den_d[,3], decreasing=F),]
  
  DEN_df_c = merge(Den_d, Den_h)    
  
  DEN_d_df = cbind(DEN_d_df, Den_d$Den_d)
  DEN_dh_df = cbind(DEN_dh_df, DEN_df_c$Den_d)     
  # DR
  DEN_df_c$DR_hD = round((DEN_df_c$Den_h - DEN_df_c$Den_d), 0)
  DR_hD_df = cbind(DR_hD_df, DEN_df_c$DR_hD) 
      
      x11()
      par(mfrow=c(2,1))
      
      plot.ts(DEN_df_c$Den_h)
      lines(DEN_df_c$Den_d, col=2)
      
      plot.ts(DEN_df_c$DR_hD)

}# loop through trees

colnames(DEN_d_df) = colnames(DR_hD_df) = 
colnames(DEN_dh_df) = colnames(DEN_selection)
rownames(DEN_dh_df) = rownames(DEN_h_df) 

DEN_d_df=data.frame(DEN_d_df)
DEN_d_diff_df = apply(DEN_d_df, 2, diff, lag = 1, differences = 1)

DR_hD_df=data.frame(DR_hD_df)



# save it ----------------------------------------------------------------------
save(list=c("DEN_h_df",  # mean hourly dendrometer fluctuation (original measurements) [um h-1]
            "DEN_h_df_diff", # hourly differences [um h-1]
            "DR_hD_df",  # King et al AFM: DEN_h-DEN_d [um h-1]
            "GRO_zwf",   # Zweifel's growth [um h-1]
            "GRO_zwf_diff", # first order differences from Zweifel [um h-1]
            "DEN_d_df", # daily stem fluctuations
            "DEN_d_diff_df", # lag-1 differences in daily stem fluctuations
            "DEN_dh_df"
            ), 
     file="D:/Dropbox/2019DendroVAR/Data/FINAL/postpro/Dendro_Data_GS_clipped.RData")
# ------------------------------------------------------------------------------


#########
# 2. calculate the Decayvars ######################################################
#########

rm(list=ls())

load("D:/postpro/Dendro_Data_GS_clipped.RData")

# load functions
source("D:/VAR_decay_function_fast.r")

# parallel computations for DecayVAR
detectCores()
cl <- makeCluster(8)

DecayVAR_R_h_diff = parApply(cl, DEN_h_df_diff, 2, VAR_decay_fast)    
DecayVAR_GRO_zwf_diff = parApply(cl, GRO_zwf_diff, 2, VAR_decay_fast)    
DecayVAR_DR_h = parApply(cl, DR_hD_df, 2, VAR_decay_fast)    
DecayVAR_R_d_diff = parApply(cl, DEN_d_diff_df, 2, VAR_decay_fast)    
                                                               
# transform from 'list' to 'matrix'
rows_diff <- max(sapply(DecayVAR_R_h_diff, function(x) max(x[[1]])))
columns <- length(DecayVAR_R_h_diff)
# make an empty matrix
DecayVAR_R_h_diff_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_R_h_diff_matr) = 1:rows_diff
colnames(DecayVAR_R_h_diff_matr) = colnames(DEN_h_df_diff)

rows_diff <- max(sapply(DecayVAR_GRO_zwf_diff, function(x) max(x[[1]])))
columns <- length(DecayVAR_GRO_zwf_diff)
# make an empty matrix
DecayVAR_GRO_zwf_diff_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_GRO_zwf_diff_matr) = 1:rows_diff
colnames(DecayVAR_GRO_zwf_diff_matr) = colnames(GRO_zwf_diff)

rows_diff <- max(sapply(DecayVAR_DR_h, function(x) max(x[[1]])))
columns <- length(DecayVAR_DR_h)
# make an empty matrix
DecayVAR_DR_h_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_DR_h_matr) = 1:rows_diff
colnames(DecayVAR_DR_h_matr) = colnames(DR_hD_df)

rows_diff <- max(sapply(DecayVAR_R_d_diff, function(x) max(x[[1]])))
columns <- length(DecayVAR_R_d_diff)
# make an empty matrix
DecayVAR_R_d_diff_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_R_d_diff_matr) = 1:rows_diff
colnames(DecayVAR_R_d_diff_matr) = colnames(DEN_d_diff_df)

# fill the matrix
for (i in colnames(DecayVAR_R_h_diff_matr)){# loop through trees
  DecayVAR_R_h_diff_matr[DecayVAR_R_h_diff[[i]]$d,i] <- DecayVAR_R_h_diff[[i]]$v
  DecayVAR_GRO_zwf_diff_matr[DecayVAR_GRO_zwf_diff[[i]]$d,i] <- DecayVAR_GRO_zwf_diff[[i]]$v 
  DecayVAR_DR_h_matr[DecayVAR_DR_h[[i]]$d,i] <- DecayVAR_DR_h[[i]]$v
  DecayVAR_R_d_diff_matr[DecayVAR_R_d_diff[[i]]$d,i] <- DecayVAR_R_d_diff[[i]]$v
}# loop through trees


x11()
matplot(x=log10(1:nrow(DecayVAR_R_h_diff_matr)), 
        y=log10(DecayVAR_R_h_diff_matr), type="l")
x11()
matplot(x=log10(1:nrow(DecayVAR_GRO_zwf_diff_matr)), 
        y=log10(DecayVAR_GRO_zwf_diff_matr), type="l")
x11()
matplot(x=log10(1:nrow(DecayVAR_DR_h_matr)), 
        y=log10(DecayVAR_DR_h_matr), type="l")
x11()
matplot(x=log10(1:nrow(DecayVAR_R_d_diff_matr)), 
        y=log10(DecayVAR_R_d_diff_matr), type="l")
        
# calculate standard deviation
Decaystdev_R_h_diff_matr = sqrt(DecayVAR_R_h_diff_matr)
Decaystdev_GRO_zwf_diff_matr = sqrt(DecayVAR_GRO_zwf_diff_matr)
Decaystdev_DR_h_matr = sqrt(DecayVAR_DR_h_matr)
Decaystdev_R_d_diff_matr = sqrt(DecayVAR_R_d_diff_matr)

x11()
matplot(x=log10(1:nrow(Decaystdev_R_h_diff_matr)), 
        y=log10(Decaystdev_R_h_diff_matr), type="l")

x11()
plot(x=log10(1:nrow(Decaystdev_R_h_diff_matr)), 
        y=log10(Decaystdev_R_h_diff_matr[,"S19Ad_L5p"]), type="l")
        
x11()
matplot(x=log10(1:nrow(Decaystdev_GRO_zwf_diff_matr)), 
        y=log10(Decaystdev_GRO_zwf_diff_matr), type="l")
x11()
matplot(x=log10(1:nrow(Decaystdev_DR_h_matr)), 
        y=log10(Decaystdev_DR_h_matr), type="l")
x11()
matplot(x=log10(1:nrow(Decaystdev_R_d_diff_matr)), 
        y=log10(Decaystdev_R_d_diff_matr), type="l")
        
                
# organise the Decaystdev dateset
test_R_h_diff = melt(Decaystdev_R_h_diff_matr)

test_R_h_diff$Species = "Picea abies"
test_R_h_diff[substr(test_R_h_diff$Var2, 7,7) == "L", "Species"] = "Larix decidua"

test_R_h_diff$Elevation = NA
test_R_h_diff[substr(test_R_h_diff$Var2, 2,3) == "08", "Elevation"] = 800
test_R_h_diff[substr(test_R_h_diff$Var2, 2,3) == "13", "Elevation"] = 1300
test_R_h_diff[substr(test_R_h_diff$Var2, 2,3) == "16", "Elevation"] = 1600
test_R_h_diff[substr(test_R_h_diff$Var2, 2,3) == "19", "Elevation"] = 1900
test_R_h_diff[substr(test_R_h_diff$Var2, 2,3) == "22", "Elevation"] = 2200
colnames(test_R_h_diff) = c("Hours", "ID", "Stdev", "Species", "Elevation")

test_GRO_zwf_diff = melt(Decaystdev_GRO_zwf_diff_matr)

test_GRO_zwf_diff$Species = "Picea abies"
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$Var2, 7,7) == "L", "Species"] = "Larix decidua"

test_GRO_zwf_diff$Elevation = NA
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$Var2, 2,3) == "08", "Elevation"] = 800
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$Var2, 2,3) == "13", "Elevation"] = 1300
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$Var2, 2,3) == "16", "Elevation"] = 1600
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$Var2, 2,3) == "19", "Elevation"] = 1900
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$Var2, 2,3) == "22", "Elevation"] = 2200
colnames(test_GRO_zwf_diff) = c("Hours", "ID", "Stdev", "Species", "Elevation")

test_DR_h = melt(Decaystdev_DR_h_matr)

test_DR_h$Species = "Picea abies"
test_DR_h[substr(test_DR_h$Var2, 7,7) == "L", "Species"] = "Larix decidua"

test_DR_h$Elevation = NA
test_DR_h[substr(test_DR_h$Var2, 2,3) == "08", "Elevation"] = 800
test_DR_h[substr(test_DR_h$Var2, 2,3) == "13", "Elevation"] = 1300
test_DR_h[substr(test_DR_h$Var2, 2,3) == "16", "Elevation"] = 1600
test_DR_h[substr(test_DR_h$Var2, 2,3) == "19", "Elevation"] = 1900
test_DR_h[substr(test_DR_h$Var2, 2,3) == "22", "Elevation"] = 2200
colnames(test_DR_h) = c("Hours", "ID", "Stdev", "Species", "Elevation")

test_R_d_diff = melt(Decaystdev_R_d_diff_matr)

test_R_d_diff$Species = "Picea abies"
test_R_d_diff[substr(test_R_d_diff$Var2, 7,7) == "L", "Species"] = "Larix decidua"

test_R_d_diff$Elevation = NA
test_R_d_diff[substr(test_R_d_diff$Var2, 2,3) == "08", "Elevation"] = 800
test_R_d_diff[substr(test_R_d_diff$Var2, 2,3) == "13", "Elevation"] = 1300
test_R_d_diff[substr(test_R_d_diff$Var2, 2,3) == "16", "Elevation"] = 1600
test_R_d_diff[substr(test_R_d_diff$Var2, 2,3) == "19", "Elevation"] = 1900
test_R_d_diff[substr(test_R_d_diff$Var2, 2,3) == "22", "Elevation"] = 2200
colnames(test_R_d_diff) = c("Days", "ID", "Stdev", "Species", "Elevation")

R_h_diff_stdev = test_R_h_diff
GRO_zwf_diff_stdev = test_GRO_zwf_diff
DR_h_stdev = test_DR_h
R_d_diff_stdev = test_R_d_diff

# save it ----------------------------------------------------------------------
save(list=c("Decaystdev_R_h_diff_matr",
            "R_h_diff_stdev",
            "Decaystdev_GRO_zwf_diff_matr",
            "GRO_zwf_diff_stdev",
            "Decaystdev_DR_h_matr",
            "DR_h_stdev",
            "Decaystdev_R_d_diff_matr",
            "R_d_diff_stdev"
            ),
     file="D:/postpro/Dendro_Data_GS_clipped_Decaystdev.RData")
# ------------------------------------------------------------------------------


# prepare data for ggplot

rm(list=ls())
load("D:/postpro/Dendro_Data_GS_clipped_Decaystdev.RData")


test_DR_h = melt(Decaystdev_DR_h_matr)
colnames(test_DR_h) = c("Hours", "ID", "Stdev")
test_DR_h$Species = "Picea abies"
test_DR_h[substr(test_DR_h$ID, 7,7) == "L", "Species"] = "Larix decidua"
test_DR_h$Elevation = NA
test_DR_h[substr(test_DR_h$ID, 2,3) == "08", "Elevation"] = 800
test_DR_h[substr(test_DR_h$ID, 2,3) == "13", "Elevation"] = 1300
test_DR_h[substr(test_DR_h$ID, 2,3) == "16", "Elevation"] = 1600
test_DR_h[substr(test_DR_h$ID, 2,3) == "19", "Elevation"] = 1900
test_DR_h[substr(test_DR_h$ID, 2,3) == "22", "Elevation"] = 2200


test_R_h_diff = melt(Decaystdev_R_h_diff_matr)
colnames(test_R_h_diff) = c("Hours", "ID", "Stdev")
test_R_h_diff$Species = "Picea abies"
test_R_h_diff[substr(test_R_h_diff$ID, 7,7) == "L", "Species"] = "Larix decidua"
test_R_h_diff$Elevation = NA
test_R_h_diff[substr(test_R_h_diff$ID, 2,3) == "08", "Elevation"] = 800
test_R_h_diff[substr(test_R_h_diff$ID, 2,3) == "13", "Elevation"] = 1300
test_R_h_diff[substr(test_R_h_diff$ID, 2,3) == "16", "Elevation"] = 1600
test_R_h_diff[substr(test_R_h_diff$ID, 2,3) == "19", "Elevation"] = 1900
test_R_h_diff[substr(test_R_h_diff$ID, 2,3) == "22", "Elevation"] = 2200


test_GRO_zwf_diff = melt(Decaystdev_GRO_zwf_diff_matr)
colnames(test_GRO_zwf_diff) = c("Hours", "ID", "Stdev")
test_GRO_zwf_diff$Species = "Picea abies"
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$ID, 7,7) == "L", "Species"] = "Larix decidua"
test_GRO_zwf_diff$Elevation = NA
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$ID, 2,3) == "08", "Elevation"] = 800
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$ID, 2,3) == "13", "Elevation"] = 1300
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$ID, 2,3) == "16", "Elevation"] = 1600
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$ID, 2,3) == "19", "Elevation"] = 1900
test_GRO_zwf_diff[substr(test_GRO_zwf_diff$ID, 2,3) == "22", "Elevation"] = 2200


# save it ----------------------------------------------------------------------
save(list=c("test_DR_h", "test_R_h_diff", "test_GRO_zwf_diff"),
     file="D:/postpro/gg_Dendro_Data_GS_clipped_Decaystdev.RData")
# ------------------------------------------------------------------------------
