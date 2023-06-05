##########################################################
## Purpose: Functions to calculate results for a particular lambda
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, June 2022
## Last modified: June 2022
##########################################################
# Main function for the sensor combos (train4)
runCombos <- function(nRun, funs, consecObs, returnOnlyMod, bayes, probType,
                      returnOnlyMetrics){
  combos <- c("sentinel1", "landsat8, sentinel2", 
              "landsat8, sentinel2, sentinel1")
  fileNames <- c("S1", "L8S2", "All")
  
  numRun <- as.numeric(nRun)
  sensKeep <- combos[numRun]
  script <- 3
  source("scripts/args.R", local=TRUE)
  source("scripts/argsTrain.R", local=TRUE)
  
  #Bring in the residuals
  fileLoadLoc <- paste0(dataPath, "/trainingPars/train1_", fileNames[numRun])
  load(paste0(fileLoadLoc, ".Rdata"))
  
  if(bayes|probType=="total"){
    dens <- paste0(dataPath, "/densFun/funStable_N_run", nRun, ".Rdata")
    densFun <- c(dens, gsub("Stable", "Dist", dens))
  } else {
    densFun <- NULL
  }
  
  if(bayes){
    load(paste0(dataPath, "/trainingPars/allStaticProbs.Rdata"))
    allProbs <- cbind(base[, .(pointid)], 
                      data.table(prior=allProbs[,1]/10000))
    saveAdd <- "_bayes"
  } else {
    allProbs <- NA
    saveAdd <- paste0("_", probType)
  }
  
  ## test for normality (evidence against using pnorm)
  ## returned item here is vector of those points whose distribution of ewma
  ## values is considered to be normal.
  # h <- lapply(l, testNormal)
  fileSaveLoc <- paste0(gsub("train1", "train4", fileLoadLoc), 
                        "_conObs", consecObs)
  saveFile <- paste0(gsub("train4", "fullZewmaProbs/train4", fileSaveLoc), 
                    "_", windowType, "_N.csv")
  
  ## run the main code
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("oneTS", "base", "dates", "window", "probThresh", 
                      "runType", "windowType", "pointids", "consecObs", 
                      "dataPath", "staticInflation", "bayes", "allProbs", 
                      "probType", "densFun", "saveFile"), envir=environment())
  clusterExport(cl, c("eqEwma", "applyLambda", "generateProbs",
                      "calcProb", "applyBayes"), envir=.GlobalEnv)
  
  ## get observation metrics or logmod only depending on context
  parallel::parSapply(l, processLambda, oneTS, base, dates, window, probThresh, 
                      runType, windowType, pointids, consecObs, dataPath, 
                      staticInflation, bayes, allProbs, probType, densFun, 
                      saveFile, cl=cl)
  stopCluster(cl)

  # outData <- lapply(l, processLambda, oneTS, base, dates, window, 
  #                   plotProbModel=FALSE, saveProbModel, probModelFile, 
  #                   probThresh, runType, windowType, pointids, consecObs, 
  #                   dataPath, returnOnlyMod, ewmaOnly=FALSE, 
  #                   staticInflation=staticInflation, returnProbs=FALSE, bayes, 
  #                   allProbs, probType, densFun)
  
  # Finally, combine and get observation dates for disturbed training points
  obsDates <- rbindlist(lapply(pointids[!grepl("a", pointids)], getObsDates, 
                               base, dates, oneTS, window))
  obsDates[, sensor := fileNames[nRun]]
  fwrite(obsDates, paste0(fileSaveLoc, "_obsDates.csv"))
  
  return(print(paste0("Finished running ", fileNames[numRun], "!")))
}
################################################################################
#' @title processLambda
#' @description This is a wrapper for smaller functions to do the second part of
#' the analysis, where we calculate results for a particular lambda value for 
#' the ewma. Inner steps include
#' applyLambda = calculate ewma and do binary categorization per point
#' fitLogMod = Fit the logistic model using all the points together
#' firstDetectionLag = Calculate the first lag where we detect disturbance
#' retrieveMetrics = for each point, predict the probabilities and get comparison 
#' accuracy metrics, including detection lag (fast and slow), N true negatives, 
#' N obs in stable period, and the binary vector of post-disturbance obs of 
#' length 30 (either observation or day)
#' @param lambda numeric - the chosen lambda, a number from 0.1-1 (from args)
#' @param oneTS a list of residuals and observation dates by training point
#' @param base the original training data csv (loaded from args)
#' @param dates date - the full date sequence of the timeframe (loaded from args)
#' @param plotProbModel logical - should the logistic model be plotted?
#' @param saveProbModel logical - should we save the logistic model?
#' @param probModelFile string -  file name if saveProbModel = TRUE
#' @param probThresh numeric - probability threshold used to denote if an 
#' observation is considered a disturbance (loaded from args)
#' @param backfill logical - should the observations be backfilled to remove NAs?

processLambda <- function(lambda, oneTS, base, dates, window, probThresh, 
                          runType, windowType, pointids, consecObs, 
                          dataPath=NULL, staticInflation, bayes, allProbs, 
                          probType, densFun=NULL, saveFile){
  # probType = "log" or "total"
  # print(lambda)
  # Calculate ewma and do binary categorization per point
  binaryCat <- rbindlist(lapply(pointids, applyLambda, oneTS, base, lambda, 
                                dates, window, runType, windowType,
                                staticInflation=staticInflation))
  
  # Make the logistic model across all disturbed points
  # if(probType=="log"){
  #   logMod <- fitLogMod(binaryCat, plotProbModel, saveProbModel, 
  #                     probModelFile, pointids)
  # } 
  
  ## calculate probabilities and save
  fullProbs <- rbindlist(lapply(pointids, generateProbs, base, lambda, bayes, 
                                probType, binaryCat, logMod=NA, densFun, 
                                allProbs))
  saveFile <- gsub("N", lambda*100, saveFile)
  fwrite(fullProbs, saveFile)
  # 
  # # if(returnOnlyMod){
  # #   return(logMod)
  # # } else {
  #   # Get the probabilities and output accuracy metrics
  #   metrics <- lapply(pointids, retrieveMetrics, base, logMod=NA, dates, 
  #                     binaryCat, probThresh, lambda, runType, window, 
  #                     windowType, consecObs, dataPath, ewmaOnly, 
  #                     returnProbs, bayes, allProbs, probType, densFun)
  #   if(ewmaOnly){
  #     names(metrics) <- pointids
  #   } else if(!returnProbs){
  #     metrics <- rbindlist(metrics)
  #   } 
  #   return(metrics)
# }
    return("Done")
}

# res <- logMod
# RSS <- c(crossprod(res$residuals))
# MSE <- RSS / length(res$residuals)
# sqrt(MSE)
# with(summary(res), 1 - deviance/null.deviance)
# summary(res)


################################################################################
#' @title fitLogMod
#' @description Fit the logistic model using all the points together
#' @param binaryCat a dt of ewma values and binary disturbed / not disturbed 
#' (output of lapply (applyLambda))
fitLogMod <- function(binaryCat, plotProbModel, saveProbModel, probModelFile, 
                      pointids){
  # Downsample the majority class to address class imbalance
  ## The way that we are doing this is by looking at the data from the
  ## monitoring period for both the disturbed and undisturbed point. Then, we
  ## subset to only the dates that have observations in both points and combine
  ## those.
  
  outList <- lapply(pointids[!grepl("a", pointids)], function(focalPt){
    ## get the time points that are common to both the disturbed training point
    ### and its undisturbed cousin
    datDist <- binaryCat[pointid==focalPt & dist==1]
    datNodist <- binaryCat[pointid==paste0(focalPt, "a") & dist==1][, dist := 0]
    
    datesDist <- c(datDist[, date], datNodist[, date])
    datesDist <- datesDist[duplicated(datesDist)]
    
    ## combine the two 
    dt <- rbind(datDist[date %in% datesDist], datNodist[date %in% datesDist])
    
    return(dt)
  })
  
  modDT <- rbindlist(outList)
  
  ## possibly do logmod based on extremes
  nodist <- modDT[dist == 0, ]
  dist <- modDT[dist == 1, ]

  nodistN <- 0.3
  distN <- 0.7
  
  # modDT <- rbind(nodist, dist)
  
  modDT <- rbind(nodist[ewma >= quantile(ewma, probs=nodistN) |
                          ewma == min(ewma)],
                 dist[ewma <= quantile(ewma, probs=distN) |
                        ewma == max(ewma)])
  
  # Calculate logistic regression
  logMod <- glm(dist ~ ewma, data=modDT, family="binomial")
  
  if(plotProbModel){
    plot(modDT[,ewma], modDT[,dist], col=scales::alpha("black", alpha=0.2),
         xlab="EWMA value", ylab="Probability of disturbance")
    curve(predict(logMod, data.frame(ewma=x), type="response"), add=TRUE, 
          lwd=2, col="blue")
    
    # par(bg="black")
    # plot(modDT[,ewma], modDT[,dist], col="brown",
    #      xlab="EWMA value", ylab="Probability of disturbance",
    #      main="Logistic model relating disturbed (y=1) and \nundisturbed (y=0) values",
    #      col.lab = "white", col.main = "white", col.axis="white", fg="white")
    # curve(predict(logMod, data.frame(ewma=x), type="response"), add=TRUE,
    #       lwd=2, col="#FFCC66")
  }
  
  if(saveProbModel){
    saveRDS(logMod, file=probModelFile)
  }
  
  return(logMod)
}

################################################################################
generateProbs <- function(pt, base, lambda, bayes, probType, binaryCat, 
                          logMod=NULL, densFun=NULL, allProbs=NA){
  if(bayes|probType=="total"){
    load(gsub("N", lambda*100, densFun[1]), envir=.GlobalEnv)
    load(gsub("N", lambda*100, densFun[2]), envir=.GlobalEnv)
  }
  
  # print(pt)
  # Step 0: Disturbance dates
  dateDistIndex <- which(dates==as.numeric(base[pointid==pt, dateDist]))
  dateModelIndex <- which(dates==as.numeric(base[pointid==pt, datePre]))
  
  ## calculate the probabilities from the logistic model
  ## NOTE: We first do this only paying attention to the observations because
  ## that's how the lags are calculated (e.g. with numConsecutivePos argument)
  ## We then apply window type before calculating the 1s and 0s for
  ## false/true positives/negatives
  binPt <- binaryCat[pointid==pt , ]
  
  ## calculate probs
  ### NOTE THIS HAPPENS HERE AND BELOW (bc we backfill for days)
  if(probType=="log"){
    probs <- as.numeric(as.vector(predict(logMod, data.frame(ewma=binPt$ewma), 
                                          type="response")))
  } else if(probType=="total"){
    probs <- sapply(binPt$ewma, calcProb, prior=NULL)
  }
  
  if(bayes){
    probsBayes <- applyBayes(allProbs, pt, binPt)
  } else {
    probsBayes <- NA
  }
  
  binPt[, `:=` (probs = probs, probType=probType, probsBayes=probsBayes)]
  return(binPt)
}

################################################################################
#' @title firstDetectionLag
#' @description Interior to retrieveMetrics. Calculate the first lag where we 
#' detect disturbance
#' @param postDist numeric - a binary vector of post-disturbance observations, 
#' where 1 = the probability is above the desired threshold and 0 = vice-versa
#' @param numConsecutivePositive number of consecutive positive observations we
#' want to have present before we consider a disturbance detected
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
################################################################################
#' @title retrieveMetrics
#' @description for each point, predict the probabilities and get comparison 
#' accuracy metrics, 
#' @param pt focal training point. This is run via lapply
#' @param binaryCat a dt of ewma values and binary disturbed / not disturbed 
#' (output of lapply (applyLambda))
#' @return a DT with accuracy metrics of detection lag (fast and slow), 
#' N true negatives, N obs in stable period, and the binary vector of 
#' post-disturbance obs of length window (either observation or day)
retrieveMetrics <- function(X, base, dates, probThresh, lambda, 
                            runType, window, windowType, consecObs, dataPath, 
                            bayes, oneTS){
  # prep necessary variables
  binPt <- X
  binPt[, probType := NULL]
  pt <- unique(binPt$pointid)
  dateDistIndex <- which(dates==as.numeric(base[pointid==pt, dateDist]))
  dateModelIndex <- which(dates==as.numeric(base[pointid==pt, datePre]))
  
  if(bayes){
    binPt[, probs := NULL]
    setnames(binPt, old="probsBayes", new="probs")
  } else {
    binPt[, probsBayes := NULL]
  }
  
  # define pre and post-disturbance probs as either above or below threshold
  ## Here, a value of "1" means the prob is >= threshold, and "0" is opposite
  ## preDist == 1 (false positive)| preDist==0 (true negative)
  ## postDist == 1 (true positive)| postDist==0 (false negative)
  preDist <- ifelse(binPt[dist==0, probs] >= probThresh, 1, 0)
  postDist <- ifelse(binPt[dist==1, probs] >= probThresh, 1, 0)
  
  if(windowType=="days"){
    increase <- dateDistIndex - dateModelIndex
  } else if(windowType=="obs"){
    # need to load in all the pre-filtered observations in order to get an
    ## "increase" value that's consistent with how the "days" increase is
    ## calculated. This Rdata object is saved from the `train1` script.
    res <- oneTS[[pt]]$stdRes
    sens <- c("landsat8sr", "sentinel2l2a")
    if(length(res) > 2){
      sens <- c(sens, "sentinel1")
    } 
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
    
  } else {
    lagSlow <- NA
    lagFast <- NA
    lagNeg <- NA
    type <- "undisturbed"
  }
  
  # If we're looking at days, then we now want to backfill the probabilities
  ## to have a continuous time series
  if(windowType=="days"){
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
    
    preDist <- ifelse(binPt[dist==0, probs] >= probThresh, 1, 0)
    postDist <- ifelse(binPt[dist==1, probs] >= probThresh, 1, 0)
  }
  
  ## tn and number observations (stable period)
  tn <- length(preDist[!is.na(preDist) & preDist==0])
  nStableObs <- length(preDist[!is.na(preDist)])
  
  ## Create a table with 1s or 0s for each of day or obs of the window period
  statTab <- as.data.table(t(postDist))
  colnames(statTab) <- paste0("stat", 1:ncol(statTab))
  if(ncol(statTab) < window){
    missingCol <- setdiff(paste0("stat", 1:window), colnames(statTab))
    missingDT <- data.table(t(rep(NA, length(missingCol))))
    colnames(missingDT) <- missingCol
    statTab <- cbind(statTab, missingDT)
  }
  
  out <- data.table(point=pt, lambda=lambda, distType=type, 
                    lagFast=lagFast, lagSlow=lagSlow, lagNeg=lagNeg,tn=tn,
                    nStableObs=nStableObs)
  out <- cbind(out, statTab)
  return(out)
}

################################################################################
#' @title getObsDates
#' @description for each point, return the dates of observations over the window
#' @param pt focal training point. This is run via lapply
#' @param base original table of training data
#' @param dates full sequence of dates for the timeframe we're looking at
#' @param oneTS list of aggregated pt timeseries
#' @param window disturbance window to look over
#' @return a DT with accuracy metrics of detection lag (fast and slow), 
#' N true negatives, N obs in stable period, and the binary vector of 
#' post-disturbance obs of length 30 (either observation or day)
getObsDates <- function(pt, base, dates, oneTS, window){
  dateDistIndex <- which(dates==as.numeric(base[pointid==pt, dateDist]))
  dateModelIndex <- which(dates==as.numeric(base[pointid==pt, datePre]))
  
  ptDT <- oneTS[[pt]]$oneTS
  datesObs <- ptDT[obsDates >= dates[dateDistIndex], obsDates][1:window]
  datesFill <- ptDT[obsDates >= dates[dateDistIndex], obsDates][1]
  datesFill <- seq(datesFill, datesFill + (window-1), by=1)
  return(rbind(data.table(pointid=pt, datesObs=datesObs, fill=0),
               data.table(pointid=pt, datesObs=datesFill, fill=1)))
}
################################################################################
## archive
# if(returnOnlyMod){
#   out <- lapply(1:length(outData), function(m){
#     res <- outData[[m]]
#     RSS <- sum((res$y-predict(res, type="response"))^2)
#     MSE <- RSS / length(res$y)
#     
#     sumMod <- summary(res)
#     
#     aicc <- round(AICc(outData[[m]]), 2)
#     rmse <- round(sqrt(MSE), 4)
#     r2 <- round(with(sumMod, 1 - deviance/null.deviance), 4)
#     
#     return(data.table(aic=aicc, rmse=rmse, r2=r2, lambda=l[m]))
#   })
#   outFull <- rbindlist(out)
#   outFull[, `:=` (sens = ifelse(nRun==2, 2, 3), window=window)]
#   return(outFull)
# }