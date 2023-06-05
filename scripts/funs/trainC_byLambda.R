##########################################################
## Purpose: Functions to calculate zEWMA and probabilities for each lambda
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, June 2022
## Last modified: June 2022
##########################################################

## -----------------------------------------------------------------------------
## getObsDates = for each point, return dates of observations over the window
## generateProbs = calc probs, then apply bayes if wanted
## processLambda = for each lambda, get zEWMA, do binaryCat, then generate prob
## runCombos = wrapper to process all lambdas, and do get obs dates
## -----------------------------------------------------------------------------
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
generateProbs <- function(pt, base, lambda, bayes, probType, binaryCat, 
                          logMod=NULL, densFun=NULL, allProbs=NA){
  if(bayes|probType=="total" & !is.null(densFun)){
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
  
  binPt[, `:=` (probs = as.numeric(probs), probType=probType, 
                probsBayes=as.numeric(probsBayes))]
  return(binPt)
}
processLambda <- function(lambda, oneTS, base, dates, window,
                          runType, windowType, pointids, 
                          dataPath=NULL, staticInflation, bayes, allProbs, 
                          probType, ewmaOnly=FALSE, densFun=NULL, saveFile=""){
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
  
  if(ewmaOnly){
    return(binaryCat)
  }
  
  ## calculate probabilities and save
  fullProbs <- rbindlist(lapply(pointids, generateProbs, base, lambda, bayes, 
                                probType, binaryCat, logMod=NA, densFun, 
                                allProbs))
  
  if(length(saveFile) > 0){
    saveFile <- gsub("N", lambda*100, saveFile)
    fwrite(fullProbs, saveFile)
  }
  ### if at some point this is uncommented, will need to bring back the 
  ### consecObs argument.
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
runCombos <- function(nRun, funs, returnOnlyMod, bayes, probType,
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
    dens <- paste0(dataPath, "/trainingPars/train3_densFun/funStable_N_run", 
                   nRun, ".Rdata")
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
  fileSaveLoc <- gsub("train1", "train4", fileLoadLoc)
  saveFile <- paste0(gsub("train4", "train4_fullZewmaProbs/train4", 
                          fileSaveLoc), "_N.csv")
  
  ## run the main code
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("oneTS", "base", "dates", "window", 
                      "runType", "windowType", "pointids", 
                      "dataPath", "staticInflation", "bayes", "allProbs", 
                      "probType", "densFun", "saveFile"), envir=environment())
  clusterExport(cl, c("eqEwma", "applyLambda", "generateProbs",
                      "calcProb", "applyBayes"), envir=.GlobalEnv)
  
  ## get observation metrics or logmod only depending on context
  parallel::parSapply(l, processLambda, oneTS, base, dates, window, 
                      runType, windowType, pointids, dataPath, 
                      staticInflation, bayes, allProbs, probType, densFun, 
                      ewmaOnly=FALSE, saveFile, cl=cl)
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
  fileSaveLoc <- gsub(paste0("conObs", "[[:digit:]]"), "", fileSaveLoc)
  fwrite(obsDates, paste0(fileSaveLoc, "_obsDates.csv"))
  
  return(print(paste0("Finished running ", fileNames[numRun], "!")))
}

#######################################################################
## archive
### fit logistic model to then predict over for probs (saving in case need to
#### come back to it at some point)
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
