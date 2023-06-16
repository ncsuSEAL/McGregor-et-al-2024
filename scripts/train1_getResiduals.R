##########################################################
## Purpose: Wrapper function to read, process, and filter sensor data before
##          fitting ts models, calculating std residuals, and creating ewma.
##          Output is a binary categorization of the ewma observations
## Run medium:
##  - Mac. This barely takes a minute to run. Could do PC or hpc if really want
##  - HPC submitted job, via csh of the same name in SEAL/Ian/dissertation
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: Feb 2023
##########################################################
library(data.table)
library(parallel)

source("scripts/funs/trainA_combineDat.R") # step 1
source("scripts/funs/commonA_processDat.R") # step 1
source("scripts/funs/commonB_modZ.R") #step 2b
source("scripts/funs/commonC_dupNA.R") #step 2c
funs <- ls()

# Main function
getPointRes <- function(pt, dataPath, sensorFold, sensorFiles, bandList, 
                        vegIndex, base, dates, spikeWin, spikeThresh, 
                        spikeAmp, spikeTime, returnStats, returnHarMod, 
                        saveStats, makeSumPlots, inflation, spike2){

  dateDistIndex <- which(dates==as.numeric(base[pointid==pt, dateDist]))
  dateModelIndex <- which(dates==as.numeric(base[pointid==pt, datePre]))
  sensorNames <- names(bandList)
  
  ############## Begin Steps 1-2 #############################
  # Step 1: Read in the data, filter by QAQC, and calculate index
  ## NB this is the slowest step out of the 3 in this function. I think this is
  ### ok, though, because this code was made to mirror the functionality of the
  ### landscape application (and it's not exactly the same).
  ptData <- lapply(sensorNames, getIndex, dataPath, sensorFold, point=pt, 
                        sensorFiles, bandList, dates, vegIndex)
  names(ptData) <- sensorNames
  
  if("sentinel1" %in% sensorNames){
    names(ptData) <- c("landsat8sr", "sentinel1ASC", "sentinel2l2a")
    ptData <- append(ptData, list("sentinel1DES"=ptData$sentinel1ASC))
    ptData <- ptData[order(names(ptData))]
  }
  
  # Step 2: Fit time series models per sensor, get standardized residuals
  stablePeriod <- c(dates[1], dates[dateModelIndex])
  
  stdRes <- lapply(names(ptData), combineZ, ptData, dates, runType,
                   spikeWin, spikeThresh, spikeAmp, spikeTime, 
                   stablePeriod, returnStats, returnHarMod,
                   saveStats, makeSumPlots, inflation)
  names(stdRes) <- names(ptData)
  
  ##2c. Combine Z-scores across sensors, removing duplicates
  oneTS <- removeDups(stdRes, makeSumPlots)
  
  ## 2d. Do second pass of spike filter on aggregated residuals
  ### This (oneTS) is NOT used in the inflation factor code, so doesn't
  ## matter if it's TRUE when doing the initial calculation
  if(spike2){
    spikesAgg <- detectOutliers(dates=oneTS$obsDates, vals=oneTS$stdRes, 
                                spikeWin, spikeThresh, spikeAmp, spikeTime)
    
    oneTS <- oneTS[spikesAgg < 4, ]
  }
  
  if(makeSumPlots){
    return(list(ptData=ptData, stdRes=stdRes, oneTS=oneTS))
  } else {
    return(list(oneTS=oneTS))
  }
}

# Run main function for 3 sensor combos
runSensCombo <- function(nRun, funs, inflationType){
  combos <- c("sentinel1", "landsat8, sentinel2", 
              "landsat8, sentinel2, sentinel1")
  fileNames <- c("S1", "L8S2", "All")
  
  numRun <- as.numeric(nRun)
  sensKeep <- combos[numRun]
  script <- 1
  bayes <- FALSE #necessary for args.R
  source("scripts/args.R", local=TRUE)
  source("scripts/argsTrain.R", local=TRUE)
  
  ## save all observation dates pre-filtering 
  if(nRun==3){
    sapply(names(bandList), function(s){
      dt <- fread(file.path(dataPath, sensorFold, 
                            sensorFiles[grepl(s, sensorFiles)]))
      allDates <- dt[, dateNum := as.numeric(date)][, .(dateNum, pointid)]
      save(allDates, file=paste0(dataPath, "/allObsDates_", s, ".Rdata"))
      return("Done")
    })
  }
  
  ## bring in inflation
  inflationFile <- paste0(dataPath, "/trainingPars/train2_seasonalAdjustment", 
                          fileNames[nRun], ".Rdata")
  if(inflationType=="static" | !file.exists(inflationFile)){
    inflation <- rep(1, 365)
  } else {
    load(inflationFile)
    inflation <- vec
  }
  
  ## run the main code
  fileSaveLoc <- paste0(dataPath, "/trainingPars/train1_", fileNames[numRun])
  
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(MuMIn))
  clusterExport(cl, ls(), envir=environment())
  clusterExport(cl, c(funs), envir=.GlobalEnv)
  oneTS <- parallel::parLapply(pointids, getPointRes, dataPath, sensorFold, 
                  sensorFiles, bandList, vegIndex, base, dates, spikeWin, 
                  spikeThresh, spikeAmp, spikeTime, returnStats,
                  returnHarMod, saveStats, makeSumPlots, inflation, spike2,
                  cl=cl)
  stopCluster(cl)
  
  names(oneTS) <- pointids
  save(oneTS, file=paste0(fileSaveLoc, ".Rdata"))
  return(print(paste0("Finished prepping ", fileNames[numRun], "!")))
}

## Use "static" inflation if haven't run train2_sdinflation.R
## After, use "dynamic" inflation. It's ok to have spike2=TRUE both times bc
## its result is not used in the inflation seasonality calculation
sapply(2:3, runSensCombo, funs, inflationType="dynamic")
