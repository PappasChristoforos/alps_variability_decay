################################################################################
# TREE RING WIDTH data analyses at Ltal
# 1. read the data and save the TRW data after removing age effects
# 2. estimate tree AGE (number of rings)
# 3. calculate slopes and save data.frame

# Christoforos Pappas
# last update: 03.06.2020
################################################################################



###
# 1. Read tree ring width time series and 'remove' age effects =================
###

# clean everything
rm(list=ls())

# set working directory where the files are stored
setwd("D:/TRW/")

# get the file names
trw_files <-list.files(pattern="dia")

# get the ids for the sites
IDs <- paste(substr(trw_files, 19,23), substr(trw_files, 15,18), sep="")

# make a list where you put all sites together
trw_sites <- list()
trw_sites_alldata <- list()

for (ii in 1:length(trw_files)){# loop around sites

  # read the text file
  trw = read.table(trw_files[ii], header=T)

  # remove the first column that has the "ring year"  
  # AND devide by 2 to get the TRW from diameter
  trw_df                <-  data.frame(trw[,-1]/2)   # units [cm]
  rownames(trw_df)      <-  trw[,1]

  # make a list with the original data
  trw_sites_alldata[[IDs[ii]]]  <-  trw_df
  
  # age effects ----------------------------------------------------------------
  # 1. select time window of 1870 to 2011, i.e., 142 years
  trw_df = trw_df[rownames(trw_df)>=1870,]
  trw_df = trw_df[rownames(trw_df)<=2011,]
  # 2. select only the time series with no NAs i.e., 142 values
  SelectedTrees = sapply(trw_df, 
                         function(col) sum(complete.cases(col)) == nrow(trw_df))
  trw_df = trw_df[,SelectedTrees]
  # 3. out of these trees, remove the first 30 yr, i.e., period of 1870-1900
  trw_df = trw_df[rownames(trw_df)>=1900,]
  trw_sites[[IDs[ii]]]  <-  trw_df

}# end loop around sites

trw_sites_112yr = trw_sites


###
# 2. Estimate Tree Age =========================================================
###

# functions
list.as.matrix <- function(x, byrow=FALSE, filler=NA){
	maxlen <- max(sapply(x,length))
    xm <- sapply(x, 
	function(xs){fillen <- maxlen-length(xs)
		if (fillen>0) {c(xs,rep(filler,fillen))} else xs
	})
	if (byrow)return(t(xm)) else return(xm)
}

TreeAge = lapply(trw_sites_alldata, function(col) nrow(col)-colSums(is.na(col)))

TreeAgeSummary =c()

for (u in 1:length(names(TreeAge))){

  dummy = list.as.matrix(TreeAge[u])
  TA = data.frame(TreeCoreId = rownames(dummy), TreeAge=dummy[,1], SiteID=colnames(dummy))
  TreeAgeSummary = rbind(TreeAgeSummary, TA)

}

TreeAgeSummary$Species = substr(TreeAgeSummary$SiteID, 1, 5)
TreeAgeSummary$Elevation = substr(TreeAgeSummary$SiteID, 6, 9)
TreeAgeSummary$Elevation[TreeAgeSummary$Elevation=="0800"] = 800
TreeAgeSummary$Elevation = as.numeric(TreeAgeSummary$Elevation)

TreeAge112yr_names = lapply(trw_sites_112yr, function(col) colnames(col))

TreeAge112yr =c()

for (u in 1:length(names(TreeAge112yr_names))){ # loop through trees

  dummy = list.as.matrix(TreeAge112yr_names[u])
  TA = data.frame(TreeCoreId = dummy[,1], SiteID=colnames(dummy))
  TreeAge112yr = rbind(TreeAge112yr, TA)

} # end loop through trees

TreeAge112yrSummary = TreeAgeSummary[TreeAgeSummary$TreeCoreId %in% TreeAge112yr$TreeCoreId, ]
     

# save dataframes

save(list=c("trw_sites_112yr", "trw_sites_alldata", "TreeAgeSummary", "TreeAge112yrSummary"), file="D:/postpro/TRWs_timeseries.RData")

    
###
# 3. Calculate slopes ==========================================================
###     

# clean everything
rm(list= ls())

#  libraries
library(reshape2)
library(plyr)
library(dplyr)


# functions
source("D:/VAR_decay_function_fast.r")


`%!in%` = Negate(`%in%`)

# load tree ring width time series
load("D:/postpro/TRWs_timeseries.RData")

# create a list
slopes = c()

for (ss in 1:length(names(trw_sites_112yr))){# loop through sites

  trw_df = trw_sites_112yr[[names(trw_sites_112yr)[ss]]]

  # the following tree cores were removed due to the low quality data 
  if (names(trw_sites_112yr)[ss] == "Picea1600"){
     trw_df = trw_df[colnames(trw_df) %!in% c("S16OPS4")]
  }
  if (names(trw_sites_112yr)[ss] == "Picea1900"){
     trw_df = trw_df[colnames(trw_df) %!in% c("S19P1", "S19S2")]
  }
  
  test1 = apply(trw_df, 2, VAR_decay_fast)
  test2 = ldply(test1, rbind)  
  test2$v = sqrt(test2$v)
  test2$d = test2$d*24*365
  
  slopes1 = test2 %>% 
              group_by(.id) %>% 
                do({
                  mod = lm(log10(v) ~ log10(d), data = .)
                  data.frame(Slope = coef(mod)[2],
                             StdError = summary(mod)$coefficients[2,2],      
                             R2adj = summary(mod)$adj.r.squared)
                })
  slopes1$Site =names(trw_sites_112yr)[ss] 
  slopes = rbind(slopes,slopes1)  
} # loop through sites


slopes$Species = substr(slopes$Site, 1, 5)
slopes$Elevation = substr(slopes$Site, 6, 9)
             
slopes$Species[slopes$Species == "Larix"] = "Larix decidua"
slopes$Species[slopes$Species == "Picea"] = "Picea abies"

slopes$Elevation[slopes$Elevation == "0800"] = 800

slopes$Elevation = as.numeric(slopes$Elevation)

# slopes from individual tree ring time series
slopes_TRWs_112yr = slopes

TreeAge112yrSummary_new = TreeAge112yrSummary[,c(1,2)]
colnames(TreeAge112yrSummary_new) = c(".id", "TreeAge")     
slopes_df_full_TRWs_112yr = merge(slopes_TRWs_112yr, TreeAge112yrSummary_new, ".id")
    

# save the df
save(list=c("slopes_df_full_TRWs_112yr"),
     file="D:/postpro/slopes_df_full_TRWs_112yr.RData")
     
  