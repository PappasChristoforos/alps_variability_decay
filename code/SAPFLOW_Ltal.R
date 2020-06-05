################################################################################
# SAPFLOW data analyses
# 1. read the data, create a data table with all the sensor data together, 
#    agreegate to hourly time series, calculate hourly differences and save
#    the SAPFLOW data in .RData format
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
# 1. read and save SAPFLOW data ################################################
#########

#  working directory
setwd("D:/SAP/")

file_names = list.files()
SAP_df = c()

for (j in 1:length(file_names)){
  SAP_data1 <- read.table(file_names[j], header=T, sep="\t")  
  SAP_df = cbind(SAP_df, SAP_data1[,"SFD_dr_co_pe"]) 
}

dates = as.character(SAP_data1$Date)
dates = substr(dates, 2, 18)
timestep <- data.frame(timestamp = as.POSIXct(strptime(dates, '%m/%d/%y %H:%M:%S'), tz = "GMT"))
rownames(SAP_df) = as.character(timestep$timestamp) 
colnames(SAP_df) = substr(file_names, 1, 8)

# aggregate 15min data to 1h 
dd = data.frame(Y=substr(rownames(SAP_df), 1, 4),
                M=substr(rownames(SAP_df), 6, 7),
                D=substr(rownames(SAP_df), 9, 10),
                H=substr(rownames(SAP_df), 12, 13))

SAP_df_h = aggregate(SAP_df, mean, na.rm=T, 
                     by=list(Y=dd$Y,
                             M=dd$M,
                             D=dd$D,
                             H=dd$H))
                             
SAP_df_h = SAP_df_h[order(SAP_df_h[,1],SAP_df_h[,2],SAP_df_h[,3], SAP_df_h[,4], decreasing=F),]

SAP_df_h = SAP_df_h[,-c(1:4)]
ss = as.POSIXct(seq(ISOdate(2012, 04, 17, 15, 00, 00), 
                    ISOdate(2015, 12, 09, 14, 00, 00), 
                    "1 h"))                     
rownames(SAP_df_h) = as.character(ss)

# remove the wet site
SAP_df_h = SAP_df_h[,colnames(SAP_df_h) %!in% c("N13WAd_L", "N13WAd_S", "N13WAd_S.1", "N13WBd_L", "N13WBd_L.1", "N13WBd_S")]

# Hourly Differences -----------------------------------------------------------
SAP_df_h_diff = apply(SAP_df_h, 2, diff, lag = 1, differences = 1)


# save it ----------------------------------------------------------------------
save(list=c("SAP_df_h", "SAP_df_h_diff"), 
     file="D:/postpro/Sapflow_H.RData")
# ------------------------------------------------------------------------------


    
#########
# 2. calculate the Decayvars ######################################################
#########

rm(list=ls())

load("D:/postpro/Sapflow_H.RData")

# load functions
source("D:/VAR_decay_function_fast.r")

# parallel computations for DecayVAR
detectCores()
cl <- makeCluster(8)

DecayVAR_SAP_h_diff = parApply(cl, SAP_df_h_diff, 2, VAR_decay_fast)    
        
# transform from 'list' to 'matrix'
rows_diff <- max(sapply(DecayVAR_SAP_h_diff, function(x) max(x[[1]])))
columns <- length(DecayVAR_SAP_h_diff)
# make an empty matrix
DecayVAR_SAP_h_diff_matr = matrix(NA, nrow=rows_diff, ncol=columns)
rownames(DecayVAR_SAP_h_diff_matr) = 1:rows_diff
colnames(DecayVAR_SAP_h_diff_matr) = colnames(SAP_df_h_diff)


# fill the matrix
for (i in colnames(DecayVAR_SAP_h_diff_matr)){# loop through trees
  DecayVAR_SAP_h_diff_matr[DecayVAR_SAP_h_diff[[i]]$d,i] <- DecayVAR_SAP_h_diff[[i]]$v
}# loop through trees

x11()
matplot(x=log10(1:nrow(DecayVAR_SAP_h_diff_matr)), y=log10(DecayVAR_SAP_h_diff_matr), type="l")


# calculate standard deviation
Decaystdev_SAP_h_diff_matr = sqrt(DecayVAR_SAP_h_diff_matr)

x11()
matplot(x=log10(1:nrow(Decaystdev_SAP_h_diff_matr)), y=log10(Decaystdev_SAP_h_diff_matr), type="l")

# organise data ----------------------------------------------------------------

test_SAP_h_diff = melt(Decaystdev_SAP_h_diff_matr)
colnames(test_SAP_h_diff) = c("Hours", "ID", "Stdev")

test_SAP_h_diff$Species = "Picea abies"
test_SAP_h_diff[substr(test_SAP_h_diff$ID, 7,7) == "L", "Species"] = "Larix decidua"

test_SAP_h_diff$Elevation = NA
test_SAP_h_diff[substr(test_SAP_h_diff$ID, 2,3) == "08", "Elevation"] = 800
test_SAP_h_diff[substr(test_SAP_h_diff$ID, 2,3) == "13", "Elevation"] = 1300
test_SAP_h_diff[substr(test_SAP_h_diff$ID, 2,3) == "16", "Elevation"] = 1600
test_SAP_h_diff[substr(test_SAP_h_diff$ID, 2,3) == "19", "Elevation"] = 1900
test_SAP_h_diff[substr(test_SAP_h_diff$ID, 2,3) == "22", "Elevation"] = 2200

SAP_h_diff_stdev = test_SAP_h_diff

# save it ----------------------------------------------------------------------
save(list=c("Decaystdev_SAP_h_diff_matr",
            "SAP_h_diff_stdev"),
     file="D:/postpro/SAP_Data_GS_Decaystdev.RData")
# ------------------------------------------------------------------------------

