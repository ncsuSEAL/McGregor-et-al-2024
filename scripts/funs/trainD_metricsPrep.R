##########################################################
## Purpose: Functions for iterating over lambdas, sensors, thresholds to calc
##          lags and prep the data for easily calculating accuracy metrics
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.2.2, May 2022 (adapted from earlier script)
## Last modified: April 2023
##########################################################

## -----------------------------------------------------------------------------
## firstDetectionLag = get first day or obs when disturbance detected
## backfillProbs = convert vector of observations to be daily
## calcIncrease = get the difference btwn fastest and slowest lag
## calcLag = calculate the lag
## retrieveMetrics = main function to make full table with accuracy metrics/lags
## allLamAllPt = wrapper for ret-Metrics, for lambdas and training points
## setPar = define parallel cluster for running across points
## runSensors = wrapper to lapply retrieveMetrics over sensors
## runThresholds = wrapper to lapply using all thresholds
## runAll = full wrapper for runSensors + runThresholds
## -----------------------------------------------------------------------------
firstDetectionLag <- function(postDist, numConsecutivePositive, binPt, dateDist,
                              windowType){
  # returns the index of the first positive detection in postDist where
  # the number of consecutive detections is at least numConsecutivePositive
  # NOTE: this allows us to measure performance in a way that is more comparable 
  ## to other work
  
  ## examples to test on
  # postDist <- c(rep(0,6), rep(1, 2), rep(0, 4), rep(1, 25))
  # postDist <- c(rep(0,6), rep(0, 4), rep(1, 25), rep(1, 2))
  # postDist <- c(1, 1, 0, rep(1, 8))
  # postDist <- rep(1, 11)
  
  postDistL <- rle(postDist)$lengths
  postDistV <- rle(postDist)$values
  postDistDates <- binPt[dist==1, date]
  
  ## make sure we have enough consecutive true positives, otherwise return NA
  validDetections <- postDistL[postDistV != 0]
  if(all(validDetections < numConsecutivePositive)) return(NA)
  
  ## get the index of the rle that gives the first time we have enough
  ## consecutive obs
  validIndex <- which(postDistV != 0 & postDistL >= numConsecutivePositive)[1]
  
  ## use the validIndex to calculate the actual lag
  dateFirstDetect <- postDistDates[numConsecutivePositive]
  
  if(validIndex==1 & dateFirstDetect==dateDist){
    ## if the obs showing disturbance occurred on the same day we have it
    ## recorded from Planet, then this is day 0.
    lagFirstDetect <- 0
  } else if(windowType=="obs"){
    ## the validIndex is already the observation lag, by definition
    lagFirstDetect <- validIndex
  } else if(windowType=="days"){
    ## when using days, we need to get the lag using the dates associated with
    ## the validIndex
    if(validIndex!=1){
      dateFirstDetect <- postDistDates[sum(postDistL[1:(validIndex-1)]) + 
                                         numConsecutivePositive]
    }
    lagFirstDetect <- dateFirstDetect - dateDist
  }
  return(lagFirstDetect)
}
backfillProbs <- function(binPt, dates, dateDistIndex, window){
  missingDays <- dates[!(dates %in% binPt$date)]
  addTab <- data.table(date=missingDays, ewma=NA, dist=NA, pointid=pt, 
                       probs=NA)
  binPt <- rbind(binPt, addTab)
  binPt <- binPt[order(date)]
  binPt <- binPt[date==dates[dateDistIndex], dist := 1]
  binPt <- binPt[date==(dates[dateDistIndex] + window), dist := 2]
  binPt[, `:=` (ewma=nafill(ewma,type="locf"),
                dist=nafill(dist, type="locf"),
                probs=nafill(probs, type="locf"))]
  binPt <- binPt[!is.na(dist), ]
  return(binPt)
}
calcIncrease <- function(dateDistIndex, dateModelIndex, windowType, dataPath,
                         pt){
  if(windowType=="days") increase <- dateDistIndex - dateModelIndex
  if(windowType=="obs"){
    # need to load in all the pre-filtered observations in order to get an
    ## "increase" value that's consistent with how the "days" increase is
    ## calculated. This Rdata object is saved from the `train1` script.
    # res <- oneTS[[pt]]$stdRes
    # if(length(res) > 2) sens <- c(sens, "sentinel1")
    
    sens <- c("landsat8sr", "sentinel2l2a")
    
    ds <- lapply(sens, function(s){
      f <- list.files(dataPath, pattern="allObsDates", full.names = TRUE)
      f <- f[grepl(s, f)]
      load(f)
      return(allDates[pointid==pt, ])
    })
    ds <- rbindlist(ds)
    allDates <- sort(unique(ds$dateNum))
    
    dMod <- dates[dateModelIndex]
    dDist <- dates[dateDistIndex]
    
    increase <- max(length(allDates[allDates >= dMod & allDates < dDist]), 1)
  }
  return(increase)
}
calcLag <- function(pt, dateDistIndex, dateModelIndex, windowType, dataPath,
                    preDist, postDist, binPt, consecObs){
  ### First calculate the difference between fastest and slowest detections
  increase <- calcIncrease(dateDistIndex, dateModelIndex, windowType, dataPath,
                           pt)
  
  ### Now calculate actual lags
  if(grepl("a", pt)){
    lagSlow <- NA
    lagFast <- NA
    lagNeg <- NA
    type <- "undisturbed"
  }
  if(!grepl("a", pt)){
    if(length(postDist)<1){
      postDist <- rep(0, 10)
      type <- "no obs in window"
    } else {
      type <- "normal"
    }
    
    dateDist <- dates[dateDistIndex]
    lagFast <- firstDetectionLag(postDist, numConsecutivePositive=consecObs, 
                                 binPt, dateDist, windowType)
    lagSlow <- lagFast + increase
    
    # account for any "negative" occurrences, when dist signal is present
    ## directly prior to the recorded day of dist itself. In other words, this
    ## counts the number of consecutive 1s starting from the end of preDist.
    end <- tail(preDist, 30)
    lagNeg <- with(rle(rev(end)), lengths * values)[1]
    # lagNeg <- ifelse(lagNeg==0, NA, lagNeg)
  }
  return(list(lagSlow=lagSlow, lagFast=lagFast, lagNeg=lagNeg, type=type))
}
retrieveMetrics <- function(pt, dt, consecObs, probThresh, window, windowType, 
                            bayes, lambda, dataPath){
  binPt <- dt[pointid==pt, ]
  dateDistIndex <- which(dates==as.numeric(base[pointid==pt, dateDist]))
  dateModelIndex <- which(dates==as.numeric(base[pointid==pt, datePre]))
  
  ## Step 0: Prep data 
  if(bayes){
    binPt[, probs := NULL]
    setnames(binPt, old="probsBayes", new="probs")
  } else {
    binPt[, probsBayes := NULL]
  }
  
  ## Step 1: Calculate binary: 
  ### Was probability >= threshold (1) or not (0) for a specific day or 
  ### observation? e.g.
  ### preDist == 1 (false positive)| preDist==0 (true negative)
  ### postDist == 1 (true positive)| postDist==0 (false negative)
  preDist <- ifelse(binPt[dist==0, probs] >= probThresh, 1, 0)
  postDist <- ifelse(binPt[dist==1, probs] >= probThresh, 1, 0)
  
  ## Step 2: Calculate lag
  lags <- calcLag(pt, dateDistIndex, dateModelIndex, windowType, dataPath,
                  preDist, postDist, binPt, consecObs)
  
  ## Step 3: Other data
  ## tn and number observations (stable period)
  tn <- length(preDist[!is.na(preDist) & preDist==0])
  nStableObs <- length(preDist[!is.na(preDist)])
  
  ## Step 4: Create per-day / per-obs table
  ### 1s or 0s for each of day or obs of the window period
  if(windowType=="days"){
    dtFull <- backfillProbs(binPt, dates, dateDistIndex, window)
  } else {
    dtFull <- binPt
  }
  preDist <- ifelse(dtFull[dist==0, probs] >= probThresh, 1, 0)
  postDist <- ifelse(dtFull[dist==1, probs] >= probThresh, 1, 0)
  
  statTab <- as.data.table(t(postDist))
  
  if(nrow(statTab)==0){
    statTab <- data.table("stat1" = NA)
  } else {
    colnames(statTab) <- paste0("stat", 1:ncol(statTab))
  }
  
  if(ncol(statTab) < window){
    missingCol <- setdiff(paste0("stat", 1:window), colnames(statTab))
    missingDT <- data.table(t(rep(NA, length(missingCol))))
    colnames(missingDT) <- missingCol
    statTab <- cbind(statTab, missingDT)
  }
  
  out <- data.table(point=pt, lambda=lambda, distType=lags$type, 
                    lagFast=lags$lagFast, lagSlow=lags$lagSlow, 
                    lagNeg=lags$lagNeg, tn=tn, nStableObs=nStableObs)
  out <- cbind(out, statTab)
  return(out)
}
allLamAllPt <- function(lam, saveFile, probThresh, dataPath, window, windowType,
                        bayes, consecObs, par, cl=NULL){
  saveFile <- gsub("N", lam*100, saveFile)
  dt <- fread(saveFile)
  dt[, probType := NULL]
  lambda <- lam
  
  if(!par){
    out <- lapply(pointids, retrieveMetrics, dt, consecObs, 
                            probThresh, window, windowType, bayes, lambda,
                            dataPath)
  } else {
    out <- rbindlist(parLapply(cl, pointids, retrieveMetrics, dt, consecObs, 
                               probThresh, window, windowType, bayes, lambda,
                               dataPath))
  }
  out[, thresh := probThresh]
  return(out)
}
setPar <- function(base, dates, probThresh, window, windowType, bayes, 
                   saveFile, pointids, consecObs){
  
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("base", "dates", "probThresh", "window", "windowType", 
                      "bayes", "saveFile", "pointids", "consecObs"), 
                envir=environment())
  clusterExport(cl, c("firstDetectionLag", "backfillProbs", "calcIncrease",
                      "calcLag"), envir=.GlobalEnv)
  return(cl)
}
runSensors <- function(X, base, dates, probThresh, window, windowType, bayes, 
                       saveFile, pointids, consecObs, dataPath, cl){
  saveFile <- gsub("S", X, saveFile)
  out <- rbindlist(lapply(l, allLamAllPt, saveFile, probThresh, dataPath,
                          window, windowType, bayes, consecObs, 
                          par=TRUE, cl=cl))
  out[, sensor := X]
  return(out)
}
runThresholds <- function(th, base, dates, probThresh, window, windowType, 
                          bayes, saveFile, pointids, consecObs, dataPath, cl){
  probThresh <- th
  print(paste0("Running threshold ", th))
  clusterExport(cl, "probThresh", envir=environment())
  fullOut <- lapply(c("L8S2", "All"), runSensors, base, dates, probThresh, 
                    window, windowType, bayes, saveFile, pointids, consecObs, 
                    dataPath, cl)
  fullOut <- rbindlist(fullOut)
  return(fullOut)
}
runAll <- function(threshSeq=NULL, base, dates, probThresh, window, windowType, 
                   bayes, saveFile, pointids, consecObs, dataPath, fileOutSave){
  cl <- setPar(base, dates, probThresh, window, windowType, bayes, saveFile, 
               pointids, consecObs)
  
  if(is.null(threshSeq)){
    outFull <- lapply(c("L8S2", "All"), runSensors, base, dates, probThresh, 
                      window, windowType, bayes, saveFile, pointids, consecObs,
                      dataPath, cl)
    fileOut <- rbindlist(outFull)
  } else {
    outFull <- lapply(threshSeq, runThresholds, base, dates, probThresh, 
                      window, windowType, bayes, saveFile, pointids, consecObs,
                      dataPath, cl)
    fileOut <- rbindlist(outFull)
    fileOutSave <- gsub("Metrics_", "MetricsThresholds_", fileOutSave)
  }
  
  stopCluster(cl)
  
  if(bayes) fileOutSave <- gsub("normal", "bayes", fileOutSave)
  fwrite(fileOut, fileOutSave)
  
  return(print("Finished!"))
}
