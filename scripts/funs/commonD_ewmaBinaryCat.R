##########################################################
## Purpose: Calculate ewma and do binary categorization for each point
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: Apr 2023
##########################################################

## -----------------------------------------------------------------------------
## eqEwma = calculate ewma manually
## -----------------------------------------------------------------------------
eqEwma <- function(lambda, z, window){
  # note that here we are letting the ewma "warm up" and then discarding the 
  ## initial values. This is because the ewma is very sensitive with few values
  z <- c(rnorm(window, 0, 1), z)
  y <- rep(NA, length(z))
  y[1] <- z[1]
  for(i in 2:length(z)){
    y[i] <- lambda*z[i] + (1-lambda)*y[i-1]
  }
  return(y[(window+1):length(y)])
}

# NOTE: There is a one-liner version of this by using the following. However, we
## speed-tested this compared to eqEwma, and discovered that while `filter` is
## faster for longer `z` (i.e. >= 10000), it is slower when doing short `z` 
## many times (e.g. >= 1000). This is probably because `filter` contains 
## additional overhead when starting the function each time, but no overhead
## once it's already running.

## ewma <- as.vector(stats::filter(z * lambda, 1 - lambda, 
##                    "recursive", init = mean(z[1:30])))
################################################################################

#' @title applyLambda
#' @description Calculate ewma and do binary categorization for each point
applyLambda <- function(pt=NULL, oneTS, base=NULL, lambda, dates, window, 
                        runType, windowType, stablePeriod=NULL,
                        staticInflation=NULL){
  # Step 0: Disturbance dates
  if(runType=="train"){
    dateDistIndex <- which(dates==as.numeric(base[pointid==pt, dateDist]))
    dateModelIndex <- which(dates==as.numeric(base[pointid==pt, datePre]))
  } else {
    dateDistIndex <- which(dates==(stablePeriod[2]+1))
    dateModelIndex <- dateDistIndex-1
  }
  
  # Step 1: Calculate ewma
  if(runType=="train"){
    z <- oneTS[[pt]]$oneTS$stdRes
    t <- oneTS[[pt]]$oneTS$obsDates
  } else {
    z <- oneTS$stdRes
    t <- oneTS$obsDates
  }
  
  # need to add on the monitoring inflation factor for the training data bc
  ## we will be predicting probabilities on this (i.e. post-log model fitting)
  ## so that way it is directly comparable to the landscape application, which
  ## already has the monitoring inflation incorporated
  ### NOTE: we only apply to dist=1 because we are not including dist=2 in the
  ### accuracy metrics
  if(runType=="train"){
    resIndex <- which(t >= dates[dateDistIndex] & 
                        t <= (dates[dateDistIndex] + window))
    z[resIndex] <- z[resIndex] * 1/staticInflation
    # bin[, `:=` (ewmaMCorr = ifelse(dist==1, ewma*(1/staticInflation), ewma), 
    #             pointid = pt)]
  }
  
  ewma <- eqEwma(lambda, z, window)
  
  # if(windowType=="days"){
  #   ptEWMA <- dates
  #   ptEWMA[!(ptEWMA %in% t)] <- NA
  #   ptEWMA[!is.na(ptEWMA)] <- ewma
  #   
  #   ewma <- nafill(ptEWMA, type="locf")
  #   t <- dates
  # }
  
  # Step 2: Do binary categorization 
  ## with weighted z-scores are 0 (no dist), 1 (dist), or 2 (other)
  distPre <- rep(0, sum(t <= dates[dateModelIndex]))
  distBtwn <- rep(2, sum(t > dates[dateModelIndex] & t < dates[dateDistIndex]))
  
  # Categorize based on both days and obs, then label the col as needed
  ## Categorize based on days
  distPost <- rep(1, sum(t >= dates[dateDistIndex] &
                           t < dates[dateDistIndex] + window))
  distBeyond <- rep(2, sum(t >= dates[dateDistIndex] + window))
  distDays <- c(distPre, distBtwn, distPost, distBeyond)
  
  ## Categorize based on observations
  distPost <- rep(1, sum(t >= dates[dateDistIndex]))
  if(length(distPost) > window) distPost[(window + 1):length(distPost)] <- 2
  distObs <- c(distPre, distBtwn, distPost)
  
  ## Add to output table
  bin <- data.table(date=t, ewma=ewma)
  
  if(windowType=="days") bin[, dist := distDays]
  if(windowType=="obs") bin[, dist := distObs]
  
  if(runType=="train") bin[, `:=` (pointid = pt)]
  
  return(bin)
}
