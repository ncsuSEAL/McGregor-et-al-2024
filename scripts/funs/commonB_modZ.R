##########################################################
## Purpose: Functions to filter data, apply spike filter, fit ts models, and 
##          return the Z-scores per sensor
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: May 2022
##########################################################

## -----------------------------------------------------------------------------
##
## reformat = reformat all data into a per-pixel list
## checkNumberModelObs = Get the number of observations used for the ts models
## createModPars = Create parameters to be used in building the ts models
## fitModels = Fit the 4 ts models and choose the best one via AIC
## getZsc = Calculate z-scores (standardized residuals) from the best model
## getModStats = Get model statistics (RMSE, Adj R-sq, etc) from the best model
## detectOutliers = spike filter as adapted from Izzi's code
## applySpike = apply the spike filter in detectOutliers
## calculateZsc = calc all z-scores and return in list
## combineZ = get z-scores from all points/pixels and combine
##
## -----------------------------------------------------------------------------
checkNumberModelObs <- function(y, t, stable){
  preObs <- as.numeric(sum(!is.na(y)[t <= stable[2]]))
  postObs <- as.numeric(sum(!is.na(y)[t > stable[2]]))
  
  if(is.na(preObs) | is.na(postObs)){
    return(list(threshold=NA, nObs=NA))
  } else {
    pre <- ifelse(preObs >= 10, 1.5, 0)
    post <- ifelse(postObs >= 10, 1, 0)
    return(list(threshold=pre+post, nObs=preObs+postObs))
  }
}
createModPars <- function(y, t, stable=c()){
  
  # Make dates NA where y is NA, then remove them
  t[is.na(y)] <- NA
  y <- y[!is.na(y)]
  t <- t[!is.na(t)]
  
  # Prep data for model runs
  y_mod <- y[t >= stable[1] & t <= stable[2]]
  t_mod <- t[t >= stable[1] & t <= stable[2]]
  
  return(list(y_mod, t_mod))
}
fitModels <- function(modPars, weightVals){
  y_mod <- modPars[[1]]
  t_mod <- modPars[[2]]
  
  # Run different sinusoidal models and choose lowest AIC. If >1 dAIC <=2, then
  ## choose most parsimonious
  fitMean <- lm(y_mod ~ t_mod, weights=weightVals)
  fit2 <- lm(y_mod ~ sin((2*pi/365)*t_mod) + cos((2*pi/365)*t_mod) + 
               t_mod, weights=weightVals)
  fit4 <- lm(y_mod ~ sin((2*pi/365)*t_mod) + cos((2*pi/365)*t_mod) +
               sin((4*pi/365)*t_mod) + cos((4*pi/365)*t_mod) + 
               t_mod, weights=weightVals)
  fit6 <- lm(y_mod ~ sin((2*pi/365)*t_mod) + cos((2*pi/365)*t_mod) +
               sin((4*pi/365)*t_mod) + cos((4*pi/365)*t_mod) +
               sin((6*pi/365)*t_mod) + cos((6*pi/365)*t_mod) + 
               t_mod, weights=weightVals)
  aics <- sapply(list(fitMean, fit2, fit4, fit6), AICc)
  
  aicTab <- data.table(num=1:4, aic=aics)
  aicTab <- aicTab[order(aic), ][, diff := abs(aic-aic[1])][diff <= 2, ]
  numModel <- as.numeric(min(aicTab$num))
  
  best_model <- list(fitMean, fit2, fit4, fit6)[[numModel]]
  
  return(best_model)
}
getZsc <- function(best_model, y, t, modelParams, inflation, runType,
                   staticInflation=NULL){
  # need the sd of the residuals ONLY in stable period (contained in best_model)
  # and use that to standardize residuals over whole series
  
  y_mod <- modelParams[[1]]
  t_mod <- modelParams[[2]]
  
  pred_t <- seq(range(t)[1], range(t)[2])
  
  #if there is too little data to make a model, then 
  ## all we can do is run with the mean
  if(length(best_model)==1){
    predval <- rep(mean(y_mod), length(t_mod))
    stanD <- sd(y_mod-predval)
  } else {
    predval <- predict(best_model, newdata=list("t_mod"=pred_t))
    stanD <- sd(best_model$residuals)
  }
  
  # calculate standardized residuals (z-score)
  ## Embedded in here is way to vary the SD by date, such that there's higher SD
  ## in the dry season (January - May)
  if(any(duplicated(t))){
    tD <- c(unique(t), t[duplicated(t)])
    predD <- c(predvalObs, predval[pred_t %in% t[duplicated(t)]])
    predvalObs <- predD[order(tD)]
    # validDates <- as.Date(tD, origin=as.Date("1970-01-01"))
  } else {
    predvalObs <- as.numeric(predval[pred_t %in% t])
    # validDates <- as.Date(pred_t[pred_t %in% t], origin=as.Date("1970-01-01"))
  }
  
  ## match the inflation factor to the doy
  doy <- yday(as.Date(t, origin=as.Date("1970-01-01")))
  stanDCorr <- stanD * inflation[doy]

  ### adjust monitoring period by static factor w/o or on top of seasonal adj
  if(runType=="app"){
    # stanDCorr[!(t %in% t_mod)] <- rep(stanD*2.29, length(t[!(t %in% t_mod)]))
    stanDCorr[!(t %in% t_mod)] <- stanDCorr[!(t %in% t_mod)] * staticInflation 
  }

  zSc <- sapply(1:length(y), function(valNum){
    return(as.numeric((y[valNum] - predvalObs[valNum]) / stanDCorr[valNum]))
    # return(as.numeric((y[valNum] - predvalObs[valNum]) / stanD))
  })
  
  return(zSc)
}
getModStats <- function(best_model, t, modelParams){
  y_mod <- modelParams[[1]]
  t_mod <- modelParams[[2]]
  
  pred_t <- seq(range(t)[1], range(t)[2])
  
  if(length(best_model)==1){
    predval <- rep(mean(y_mod), length(t_mod))
    
    modVars <- list("RMSE" = round(sqrt(mean((y_mod-predval[pred_t %in% 
                                                              t_mod])^2)),5), 
                    "AdjR2" = round(summary(fitMean)$adj.r.squared,5),
                    "nModelObs" = length(y_mod))
    #note that we can use the mean model itself for the Adj R2 bc the R2 is not
    ##dependent on model fit; it is only reporting the correlation of the predictor
    ##and the response. (Also, I manually checked the math and I got the same answer)
  } else {
    modVars <- list("RMSE" = round(sqrt(mean(best_model$residuals^2)),5), 
                    "AdjR2" = round(summary(best_model)$adj.r.squared, 5),
                    "nModelObs" = length(y_mod))
  }
  
  # for(i in 1:3){
  #   var <- ifelse(i==1, "rmse", ifelse(i==2, "adjr2", "nObs"))
  #   out_file <- paste0(modStatsFile, var, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  #   
  #   # write the data to a file as binary integers
  #   if(!file.exists(out_file)){
  #     ff <- file(out_file, 'wb') # create the file and open for writing
  #   }else{
  #     ff <- file(out_file, 'ab') # file exists, append to the end
  #   }
  #   
  #   if(i==1) vals <- modVars$RMSE else if(i==2) vals <- modVars$AdjR2 else vals <- modVars$nModelObs
  #   writeBin(vals, ff)
  #   close(ff)
  #   
  #   ## to inspect the inside of this file, need to open the connection
  #   # ff <- file(out_file, 'rb')
  #   # dat <- readBin(ff, numeric(), n=10) # where "n=10" is the number of instances you want to see
  #   # close(ff)
  # }
  return(modVars)
}
detectOutliers <- function(dates, vals, spikeWin, spikeThresh, 
                           spikeAmp, spikeTime) {
  w_floor <- floor(spikeWin / 2) # number of obs to consider before/after
  #                              date of interest
  n_obs <- length(vals)
  if (n_obs < spikeWin) {
    warning("Window size is larger than the number of dates
         with observations, thus no spikes detected.")
    return(rep(0, length(vals))) # no spike
  }
  
  spikes <- sapply(1:length(vals), function(i){
    if(i < (w_floor + 1) | i > length(vals) - w_floor){
      return(1)
    } else {
      # Determine observation of interest and window before/after
      center <- vals[i] # Center observation of interest
      pre <- vals[(i - w_floor):(i - 1)] # Obs before date of interest
      post <- vals[(i + 1):(i + w_floor)] # Obs after date of
      amp <- mean(c(abs(pre-post)))
      
      # Calculate diffs between the median values before/after central obs
      pre_diff <- abs(center - median(pre))
      post_diff <- abs(median(post) - center)
      
      # spikeThresh is the deviation factor we've assigned.
      ## spike_thresh now gives us the absolute value of the deviation,
      ### here, representing the max of the standard deviation of the pre
      ### and post values. We now say that for a spike to be considered,
      ### one check is if the difference between the pre and post anomaly is
      ### larger than this spike_thresh.
      spike_thresh <- max(spikeThresh * sd(c(pre, post), na.rm = TRUE), 0.1)
      
      # number codes for filtering
      ## 1 = no spike (black)
      ## 2 = passed sd test, failed amplitude test (red)
      ## 3 = failed sd test, passed amplitude test (green)
      ## 4 = true spike (light blue)
      
      if(max(dates[(i + 1):(i + w_floor)]) -
         min(dates[(i - w_floor):(i - 1)]) <= spikeTime){
        spT <- ifelse(pre_diff + post_diff >= spike_thresh, 1.5, 0.5)
        spA <- ifelse(pre_diff >= amp*spikeAmp &
                        post_diff >= amp*spikeAmp, 2.5, 0.5)
        return(sum(spT, spA))
      } else {
        return(5)
      }
    }
  })
  
  # filter out the spikes
  ## 0 = no spike, not enough values to be filtered over for spike 
  ## 1 = no spike (black)
  ## 2 = passed sd test, failed amplitude test (red)
  ## 3 = failed sd test, passed amplitude test (green)
  ## 4 = true spike (light blue)
  
  return(spikes)
}
calculateZsc <- function(y, t, stablePeriod, returnStats=FALSE,
                         returnHarMod=FALSE, inflation, runType,
                         staticInflation=NULL){
  #must be equal to >1 to be valid for rest of process, otherwise assign NAs
  nObs <- checkNumberModelObs(y, t, stable=stablePeriod) 
  
  if(nObs$threshold > 1) {
    # Fit models the first time and get the best one
    modPars <- createModPars(y, t, stable=stablePeriod)
    mod <- fitModels(modPars, weightVals=rep(1, length(modPars[[1]])))
    
    # Calculate weights for weighted least squares and refit the models
    ## https://rpubs.com/mpfoley73/500818
    # weightMod <- 1 / lm(abs(mod$residuals) ~ mod$fitted.values)$fitted.values^2
    # mod <- fitModels(modPars, weightVals=weightMod)
    
    # Calculate standardized residuals from the best model after weighted LS
    modZ <- getZsc(best_model=mod, y, t, modelParams=modPars, inflation, 
                   runType, staticInflation)
    
    if(returnStats){
      modStats <- getModStats(mod, t, modPars)
    } else {
      modStats <- list("RMSE" = NA, "AdjR2" = NA, "nModelObs" = NA)
    }
    
    if(returnHarMod) modOut <- mod else modOut <- NULL
    
  } else {
    modZ <- rep(NA, length(y))
    modStats <- list("RMSE" = NA, "AdjR2" = NA, "nModelObs" = nObs$nObs)
    modOut <- NULL
  }
  
  return(list("modZ"=modZ, "modStats"=modStats, "model"=modOut))
}
combineZ <- function(U, ptData, dates, runType,
                     spikeWin, spikeThresh, spikeAmp, spikeTime, 
                     stablePeriod, returnStats=FALSE, returnHarMod=FALSE, 
                     saveStats="", makeSumPlots=FALSE, inflation, 
                     staticInflation=NULL){
  # U <- "sentinel1"
  # nPixel=1
  
  t <- dates
  
  # if(combined & !grepl("sentinel1", U)){
  #   dtComb <- data.table(t=t, y=ptData[["landsat8sr"]]$ts, s="landsat8sr")
  #   dtComb <- rbind(dtComb, 
  #                   data.table(t=t, y=ptData[["sentinel2l2a"]]$ts, 
  #                              s="sentinel2l2a"))
  #   
  #   dtComb <- dtComb[order(t), ][!is.na(y), ]
  #   keys <- c("t")
  #   dups <- dtComb[,list(y= mean(y)), keys]
  #   
  #   setkeyv(dups, "y")
  #   setkeyv(dtComb, "y")
  #   dupVals <- dtComb[dups
  #                     ][, `:=` (t=NULL, s=ifelse(is.na(s), "combined", s))
  #                       ][order(i.t)]
  #   setnames(dupVals, old=c("i.t"), new=c("t"))
  #   
  # }
  
  if(runType=="train"){
    y <- ptData[[U]]$ts
  } else if(runType=="app"){
    y <- ptData[[U]]$ts
    t <- t[[U]]
  }
  
  if(grepl("sentinel1", U)){
    d <- ptData[[U]]$dir
    d <- d[!is.na(y)]
    # if(length(d[d=="ASCENDING"]) > length(d[d=="DESCENDING"])){
    if(grepl("ASC", U)){
      dir <- "ASCENDING"
    }   else {
      dir <- "DESCENDING"
    }
  }
  
  # Remove NAs from t and y
  t <- t[!is.na(y)]
  y <- y[!is.na(y)]
    
  # For Sentinel1 again, remove undesired path
  if(grepl("sentinel1", U)){
    # have interesting case where the fill value appears to be 0
    t <- t[y != 0 & d == dir]
    y <- y[y != 0 & d == dir]
  }
  
  ## if use the "combined" tag again, then this can't be done
  # Remove duplicates by date by taking the mean of obs when duplicates occur
  dupVals <- data.table(t=t, y=y)
  keys <- "t"
  dupVals <- dupVals[,list(y=mean(y)), keys]
  
  # Retain original t values for "app"
  # if(runType=="app") tOrig <- t
  
  # Apply spike filter (see function for description of filtering)
  spikes <- detectOutliers(dates=dupVals$t, vals=dupVals$y, spikeWin,
                           spikeThresh, spikeAmp, spikeTime)

  t <- dupVals$t[spikes != 4]
  y <- dupVals$y[spikes != 4]
  
  # Fit models and get z-scores
  outZ <- calculateZsc(y, t, stablePeriod, returnStats=returnStats, 
                       returnHarMod=returnHarMod, inflation, runType,
                       staticInflation)
  
  outList <- list(stdResid=outZ$modZ, obsDates=t)
  
  if(grepl("sentinel1", U)) outList <- append(outList, list(dir=dir))
  
  if(returnStats) outList <- append(outList, list(modStats=outZ$modStats))
  if(returnHarMod) outList <- append(outList, list(model=outZ$model))
  if(makeSumPlots){
    outList <- append(outList, list(spikes=spikes, spikeVals=dupVals))
    
    # if(combined){
    #   dataSource <- dupVals$s
    #   outList <- append(outList, list(dataSource=dataSource))
    # }
  }
  
  # if(runType=="app"){
  #   return(list(zSc=absZ, dates=t, modStats=outZ$modStats, model=outZ$model, 
  #               spikes=spikes, spikeVals=dupVals, tOrig=tOrig, dataSource=dataSource))
  # }
  return(outList)
}
