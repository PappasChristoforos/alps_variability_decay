################################################################################
# A function for estimating how the variabiliy (variance estimate; v) of a time series
# changes with averaging time scale (d)

# Christoforos Pappas
# last update: 03.06.2020
################################################################################

VAR_decay_fast <- function(timeseries){

  # load libraries inside for the 'parallel' package
  library(data.table)

  nna <- sum(!is.na(timeseries)) # actual length without accounting for missing values

  # aggregrion scale up to 10% of the sample size
  # delta <- round(0.1*nrow(site_df),0)
  # do not count NAs
  delta <- round(0.1*nna,0)

  if (delta!=0){
    # Create a data.frame to store all the variances
    vardf <- data.frame(d=1:delta, v=NA)

    # subtract the mean of all data points from each individual data point, then
    # divide those points by the standard deviation of all points
    samp <-data.table(scale(timeseries, center = TRUE, scale = TRUE))
    vardf[1,"v"]    <- var(samp, na.rm=T)

    for (i in 2:delta){# loop around aggregation scale
      dummy <- rep(1:floor(nna/i), each=i)
      ttest <- samp[1:length(dummy),list(AVERAGE=mean(V1)),by=list(dummy)]$AVERAGE

      if (sum(!is.na(ttest))>9){ # have at least 10 values for the variance estimation
       dummyvar <- var(ttest, na.rm=T)
       vardf[i,"v"] <- dummyvar
      }
    }# loop around aggregation scale

  }else{
   vardf <- data.frame(d=NA, v=NA)
  }

 return(vardf)
}
