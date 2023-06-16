##########################################################
## Purpose: Functions for second step of landscape application
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: June 2022
##########################################################
## -----------------------------------------------------------------------------
##
## reformat = convert large list of sensor matrices into list of each pixel
## with its constituent ts from each sensor
## saveProbs = write out the prob ts for each pixel per day to a binary file
## processPix = processing steps for each pixel, including fitting ts models,
##              calculating reiduals, aggregating to single ts, remove dups,
##              calculate ewma, do binary cat, and calculate probs
##
## -----------------------------------------------------------------------------

reformat <- function(X, dat){
  out <- lapply(names(dat), function(Y, dMat=dat){
    ts <- dMat[[Y]]$ts[X, ]
    
    if(Y=="sentinel1"){
      # have to do conditional here because Chatthin is covered by 1 S1 tile, 
      ## so directions are all the same.
      ## Will need to change this when applying to areas that overlap
      dir <- dMat[[Y]]$dir
      return(list(ts=ts, dir=dir))
    } else {
      return(list(ts=ts))
    }
  })
  names(out) <- names(dat)
  
  return(out)
}
fillDates <- function(value, marker, dates, full){
  outFull <- dates
  outFull[!(dates %in% marker)] <- NA
  outFull[!is.na(outFull)] <- value
  
  if(full) outFull <- nafill(outFull, type="locf")
  return(outFull)
}
createProbsFolders <- function(saveDir, probsDir, dates){
  outFold <- probsDir
  
  folderPath <- paste0(saveDir, outFold)
  if(!dir.exists(folderPath)) dir.create(folderPath)

  fold <- list.files(folderPath)
  days <- paste0("day", seq(1:length(dates)))
  foldToAdd <- days[!(days %in% fold)]
  if(length(foldToAdd) > 0){
    sapply(foldToAdd, function(dayF) dir.create(file.path(folderPath, dayF)))
  }
  return(print("All folders accounted for"))
}
saveProbs <- function(dayN, allDat, outFile, landscapeVal, ewmaFull=FALSE){
  ## We MUST do integers because of how much data we have. So we multiply out 
  ## first, though note that 10000 is the MAX in terms of memory
  ## Here we are saving each day's probabilities (the column dayN) into a vector
  ## that is saved as binary. This way when we build the daily maps, we pull in
  ## the vectors one at a time to create the full landscape matrix.
  
  dat <- allDat[, dayN]
  dat <- as.integer(dat*10000) #need to do this even for binary
  
  # append data to single file
  if(landscapeVal){
    out_file <- gsub("/zEWMAfull", 
                     paste0("/day", dayN, "/zEWMAfull"), outFile)
  } else if(ewmaFull){
    out_file <- gsub("/probsBin", paste0("/day", dayN, "/probsBin"), outFile)
  } else {
    out_file <- gsub("/allProbs", paste0("/day", dayN, "/allProbs"), outFile)
  }
  
  # write the prob data to a file as binary integers
  ffP <- file(out_file, 'wb') # create the file and open for writing
  writeBin(dat, ffP)
  close(ffP)
  
  # write the dates data to a file as binary integers
  # ffD <- file(out_fileDates, 'wb') # create the file and open for writing
  # writeBin(as.vector(range(expand_dates)), ffD)
  # close(ffD)
  return(print(paste0("Done with Day ", dayN, "!")))
}
updateProbs <- function(e, allProbs, ewmaFull, dMonitor, binaryDist=FALSE){
  priorProbPix <- allProbs[e]
  
  if(binaryDist){
    ewmaFullPix <- ewmaFull[, e]
  } else {
    ewmaFullPix <- ewmaFull[[e]]
  }
  
  ## initially we were only applying bayesian updating to the monitoring
  ### period, but for consistency's sake with how we're doing the rest of the
  ### analysis now, we are changing it (17 Apr 23)
  if(!is.na(dMonitor)){
    ewmaFullPix <- ewmaFullPix[dMonitor:length(ewmaFullPix)] #subset to monitoring
  }
  
  pixProbs <- sapply(1:length(ewmaFullPix), function(day){
    outNum <- calcProb(z=ewmaFullPix[day], prior=priorProbPix)
    outNum <- ifelse(outNum < 0, 0, outNum)
    return(outNum)
  })
  return(pixProbs)
}
applyBayesLand <- function(nRun, cl, ewmaFull, saveDir, probsDir, matRuns, nCol, 
                      dMonitor, txtOutFile, binaryDist=FALSE){
  ## bring in prior probs, and subset to only the cells analyzed in
  ## this nRun
  matRuns <- as.data.table(matRuns)
  if(nRun==1){
    cellRange <- c(1, nCol*matRuns[nRun, nR])
  } else {
    cellRange <- c(nCol*matRuns[(nRun-1), endR] + 1, 
                   nCol*matRuns[(nRun), endR])
  }
  
  load(paste0(saveDir, gsub("probs", "allStaticProbs.Rdata", probsDir)))
  allProbs <- allProbs[cellRange[1]:cellRange[2]]
  allProbs <- allProbs/10000
  
  ## do the math for Bayes' Theorem
  write(paste0("Starting to apply Bayes theorem math for run ", nRun, " at ",  
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)

  clusterExport(cl, c("ewmaFull", "allProbs", "dMonitor"), envir=environment())
  clusterExport(cl, c("calcProb"), envir=.GlobalEnv)
  probsUpdated <- parLapply(cl, 1:length(allProbs), updateProbs, allProbs,
                            ewmaFull, dMonitor, binaryDist)

  write(paste0("Finished applying Bayes theorem math for run ", nRun, " at ",  
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)
  return(probsUpdated)
}
processPix <- function(pix, sensorNames, dateSens, runType, window, spikeWin, 
                       spikeThresh, spikeAmp, spikeTime, stablePeriod, 
                       lambda, logMod, returnStats, returnHarMod, saveStats, 
                       makeSumPlots, dates, windowType, inflation,
                       staticInflation, spike2, getSensDates=FALSE, 
                       bayesUpdate, landscapeVal, zEWMAdates, txtOutFile){
  # write(paste0("Processing a pixel at ", 
  #              format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
  #       append=TRUE)
  
  # Step 4: Fit timeseries models and calculate st. residuals
  stdRes <- lapply(sensorNames, combineZ, ptData=pix, 
                   dates=dateSens, runType, spikeWin, spikeThresh, 
                   spikeAmp, spikeTime, stablePeriod, returnStats, returnHarMod,
                   saveStats, makeSumPlots, inflation, staticInflation)
  names(stdRes) <- sensorNames
  
  # Step 5: Aggregate / remove duplicates and NA
  oneTS <- removeDups(stdRes, makeSumPlots)
  
  ## 5a. Do second pass of spike filter on aggregated residuals
  ### This (oneTS) is NOT used in the inflation factor calculation
  if(spike2){
    spikesAgg <- detectOutliers(dates=oneTS$obsDates, vals=oneTS$stdRes, 
                                spikeWin, spikeThresh, spikeAmp, spikeTime)
    
    oneTS <- oneTS[spikesAgg < 4, ]
  }
  
  # Step 6: EWMA and binary categorization
  binaryCat <- applyLambda(pt=NULL, oneTS, base=NULL, lambda, dates, window, 
                           runType, windowType, stablePeriod)
  
  if(getSensDates){
    obsDates <- binaryCat$date
    datesMonitor <- obsDates[obsDates >= as.Date("2019-01-01") &
                               obsDates <= as.Date("2019-12-31")]
    return(datesMonitor)
  } else {
    # Step 7: Calculate non-Bayesian probs
    # probs <- as.numeric(as.vector(predict(logMod, 
    #                                       data.frame(ewma=binaryCat$ewma), 
    #                                       type="response")))
    probs <- sapply(binaryCat$ewma, calcProb, prior=NULL)
    
    ## expand probs to fill in for all dates (i.e. if we didn't have an obs for
    ## a certain date, we fill it by bringing in the last non-NA obs)
    ### NB bayes and normal are set to FALSE here on purpose
    probsFull <- fillDates(value=probs, marker=binaryCat$date, dates, full=TRUE)
    
    ## need to convert all preceding NAs to 0 in order to not have these
    ## be confused with masked, non-forested pixels
    probsFull[is.na(probsFull)] <- 0
    
    if(returnStats) stop("Need to add functionality for returning 
                       individual model stats")
    
    if(zEWMAdates){
      pixelZewma <- binaryCat$ewma
      pixelDate <- binaryCat$date
    } else {
      pixelZewma <- NA
      pixelDate <- NA
    }
    
    if(bayesUpdate|landscapeVal){
      ewmaFull <- fillDates(value=binaryCat$ewma, marker=binaryCat$date, dates, 
                            full=TRUE)
      # ewmaFull <- ewmaFull[dates >= (stablePeriod[2] + 1)]
    } else {
      ewmaFull <- NA
    }
    return(list(probsFull=probsFull, ewmaFull=ewmaFull, pixelZewma=pixelZewma, 
                pixelDate=pixelDate))
  }
}

getResids <- function(){
  # Step 4: Fit timeseries models and calculate st. residuals
  stdRes <- lapply(sensorNames, combineZ, ptData=pix, 
                   dates=dateSens, runType=runType,
                   spikeWin=spikeWin, spikeThresh=spikeThresh, 
                   spikeAmp=spikeAmp, spikeTime=spikeTime, 
                   stablePeriod=stablePeriod, 
                   returnStats=returnStats, returnHarMod=returnHarMod,
                   saveStats=saveStats, combined=combined,
                   makeSumPlots=makeSumPlots)
  names(stdRes) <- sensorNames
}

# allD <- dates
# allD[!(allD %in% oneTS$obsDates)] <- NA
# allD[!is.na(allD)] <- probs
# allD <- nafill(allD, type="locf")
# 
# plot(dates, probsFill, type="l")
# lines(dates, allD, col="red")