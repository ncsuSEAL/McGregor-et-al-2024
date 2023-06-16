##########################################################
## Purpose: Apply ts models and prob model at landscape scale to get 
##          probabilities of disturbance
## Input: all necessary functions
## Variable created: matrices of daily prob of disturbance per 
##          defined chunks, saved by day
## Run medium:
##  - HPC job, via csh file of same name located in SEAL/Ian/
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Dec 2021
## Last updated: June 2022
##########################################################
library(data.table)
library(parallel)

# setwd("Z:/IanMcGregor/")

source("scripts/funs/appC_analyzeLandscape.R")
source("scripts/funs/commonB_modZ.R")
source("scripts/funs/commonC_dupNA.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")
source("scripts/funs/commonF_probBayes.R")
funs <- ls()

args <- commandArgs(trailingOnly=TRUE)

applyLandProbs <- function(nRun, appRegion, bayesUpdate, writeProbs, 
                           returnDates, landscapeVal, zEWMAdates, nCores){
  applicationStep <- 2
  nRun <- as.numeric(nRun)
  if(bayesUpdate) bayes <- TRUE else bayes <- FALSE
  
  # Step 1. Define and load key variables / define parameters
  ## nRun is used here to load data
  source("scripts/args.R", local=TRUE)
  source("scripts/argsApp.R", local=TRUE)
  
  # Step 2. Write out progress to text file
  write(paste0("Started ", appLoc, outText, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)
  
  # Step 3. Reformat the data for parallelization
  eachPixel <- lapply(1:nrow(allMatrices$landsat8sr$ts), reformat, 
                      dat=allMatrices)
  dateSens <- lapply(allMatrices, function(X) return(X$dates))
  datesRange <- as.vector(sapply(dateSens, function(X) return(range(X))))
  dates <- seq(min(datesRange), max(datesRange), by=1)
  
  rm(allMatrices)
  print(paste0("Data has been reformatted. Moving on to timeseries ", 
                "models and probability calculation"))
  
  # Step 3B. Create probability folders if don't exist
  createProbsFolders(saveDir, probsDir, dates)
  
  if(bayesUpdate){
    createProbsFolders(saveDir, probsDir=paste0(probsDir, "Bayes"), dates)
                      # dates=(endStable+1):dates[length(dates)])
  }
  
  if(landscapeVal){
    createProbsFolders(saveDir, probsDir=gsub("probs", "ewma", probsDir), dates)
  }
  
  if(zEWMAdates){
    zewmaDir <- paste0(saveDir, gsub("probs", "ewmaPixDates", probsDir))
    if(!dir.exists(zewmaDir)) dir.create(zewmaDir)
  }
  
  # Step 3C. Load density functions for L8S2 only
  dens <- paste0(saveDir, "trainingPars/densFun/funStable_", lambda*100, 
                 "_run2.Rdata")
  densFun <- c(dens, gsub("Stable", "Dist", dens))
  load(densFun[1]); load(densFun[2])
  
  # Step 4-7: ts models, std res, aggregate, ewma, binaryCat, get probs
  stablePeriod <- c(dates[1], endStable)
  sensorNames <- names(eachPixel[[1]])

  write(paste0("Started parallel process for run ", nRun, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)
  
  cl <- parallel::makeCluster(nCores)
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(MuMIn))
  clusterExport(cl, c("sensorNames", "dateSens", "runType", "window", 
                      "spikeWin", "spikeThresh", "spikeAmp", "spikeTime",
                      "stablePeriod", "lambda", "logMod", "returnStats",
                      "returnHarMod", "saveStats", "makeSumPlots", "dates",
                      "outFileProbs", "appLoc", "windowType", "inflation",
                      "staticInflation", "spike2", "returnDates", "bayesUpdate", 
                      "landscapeVal", "txtOutFile", "densEDist", "densEStable",
                      "zEWMAdates"), 
                envir=environment())
  clusterExport(cl, c(funs))

  out <- parallel::parLapply(eachPixel, processPix, sensorNames,
                             dateSens, runType, window,
                             spikeWin, spikeThresh, spikeAmp, spikeTime,
                             stablePeriod, lambda, logMod, 
                             returnStats, returnHarMod, saveStats, 
                             makeSumPlots, dates, windowType, 
                             inflation, staticInflation, spike2, 
                             getSensDates=returnDates, bayesUpdate, 
                             landscapeVal, zEWMAdates, txtOutFile, cl=cl)
  print(paste0("Finished initial parallel processing. Moving on to next step ",
                "to save output as binary and/or apply bayes."))
  stopCluster(cl)

  write(paste0("Finished parallel process for run ", nRun, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)
  
  # bleh <- list()
  # for(m in 1:1){
  #   bleh[[m]] <- processPix(eachPixel[[m]], sensorNames,
  #   dateSens, runType, window,
  #   spikeWin, spikeThresh, spikeAmp, spikeTime,
  #   stablePeriod, lambda, logMod, 
  #   returnStats, returnHarMod, saveStats, 
  #   makeSumPlots, dates, windowType, 
  #   inflation, staticInflation, spike2, 
  #   getSensDates=returnDates, bayesUpdate,txtOutFile)
  # }

  # Step 8. Split out the model stats from the probs
  ## NOTE if actually want this then will need to add functionality to
  ### processPix to return the stats as "modStats"
  if(returnStats){
    stats <- lapply(allProbs, function(X) return(X$modStats))
    
    # This creates a list of 3 matrices, 1 for each model statistic, with the
    ## order of rows = Landsat-8, Sentinel-2, Sentinel-1
    ## note the "*9" is because for each pixel, there are 3 sensors each with 3 model statistics
    statsMats <- lapply(1:3, function(Y){
      mat <- matrix(unlist(stats)[seq(Y, length(stats)*9, by=3)], nrow=3)
      return(mat)
    })
    names(statsMats) <- c("rmse", "adjr2", "nObs")
    
    saveRDS(statsMats, paste0(modelFitsOut, "modStats.RDS"))
    # writeStatsRast(dat=stats, endR=endRow, startR=startRow, nR=numberRow, 
    #               nC=nCol, saveFile=paste0(modelFitsOut, "modStats.tif"))
  }
  
  # for validation, we want to return the exact cell numbers along with the
  ## valid observation dates (after QAQC and 2 spike filters). This is bc
  ## each cell has different valid dates post-filtering
  if(returnDates){
    saveFold <- paste0(saveDir, appLoc, "/validation")
    if(!dir.exists(saveFold)) dir.create(saveFold)
    
    ## corroborate cell number (for this run) with valid dates
    validDates <- rbindlist(lapply(1:length(out), function(X){
      return(data.table(date=out[[X]], cell=X))
    }))

    ## use the unique valid dates and reformat into list
    dUnique <- sort(unique(validDates$date))
    cellsPerDate <- lapply(dUnique, function(X){
      return(validDates[date==X, cell])
    })

    names(cellsPerDate) <- paste0("date=", dUnique)

    ## write out as Rdata
    save(cellsPerDate, file=paste0(saveFold, "/dates_start", startRow, "end", 
                                 endRow, ".Rdata"))
  } else {
    # Step 9. Save the probabilities per day in a binary file.
    ## First, bring each pixel's prob values together into a single matrix, 
    ## where each column is a different day
    # allProbs <- t(sapply(out, function(X) return(X$probDist)))
    cl <- parallel::makeCluster(nCores)
    clusterEvalQ(cl, library(data.table))
    
    ## prep output for saving
    if(writeProbs){
      allDat <- t(matrix(unlist(lapply(out, `[[`, 1)), nrow=length(dates)))
      outFile <- outFileProbs
    }
    
    ## write out zewma 
    if(zEWMAdates){
      for(w in c("pixelZewma", "pixelDate")){
        allBin <- lapply(out, '[[', w)
        
        ## make equal length vectors to make matrix
        maxL <- max(lengths(allBin))
        allBinMat <- sapply(allBin, function(X){
          vec <- round(X*10000)
          return(c(vec, rep(NA, maxL-length(vec))))
        })
        outF <- gsub("probs", "ewmaPixDates", outFileProbs)
        save(allBinMat, file=paste0(gsub("allProbs", w, outF), ".Rdata"))
      }
    }

    ## retain all zEWMA values
    if(landscapeVal){
      ewmaFull <- lapply(out, `[[`, 2)
      # dEnd <- which(dates==stablePeriod[2])
      # stableZ <- lapply(ewmaFull, function(z){
      #   return(z[1:dEnd])
      # })

      ## save zEWMA
      allDat <- t(matrix(unlist(ewmaFull), nrow=length(dates)))
      outFile <- gsub("probs/allProbs", "ewma/zEWMAfull", outFileProbs)
    }

    # Step 10. Do Bayesian updating
    if(bayesUpdate){
      ewmaFull <- lapply(out, `[[`, 2)
      dMonitor <- NA #which(dates==endStable) + 1
      if(!is.na(dMonitor)){
        nr <- length(dates[dMonitor:length(dates)])
      } else {
        nr <- length(dates)
      }
      clusterExport(cl, c("densEDist", "densEStable", "dMonitor"), 
                    envir=environment())
      allDat <- applyBayesLand(nRun, cl, ewmaFull, saveDir, probsDir, matRuns, 
                              nCol, dMonitor, txtOutFile)
      allDat <- t(matrix(unlist(allDat), nrow=nr))
      outFile <- gsub("probs/", "probsBayes/", outFileProbs)
      stopCluster(cl)

      ## make new cluster for exporting so we don't keep all the data from
      ## ewmaFull AND allDat on each core
      cl <- parallel::makeCluster(nCores)
      clusterEvalQ(cl, library(data.table))
    }

    ## Then, scale by 10000 (integer) and write to file
    write(paste0("Saving prob files for run ", nRun, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)

    clusterExport(cl, c("saveProbs", "allDat", "outFile", "landscapeVal"), 
                  envir=environment())
    parSapply(1:ncol(allDat), saveProbs, allDat, outFile, landscapeVal, cl=cl)
    
    stopCluster(cl)

    write(paste0("Finished saving prob files for run ", nRun, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)
  }
  
  print("All outputs saved. Finished with this iteration.")
  write(paste0("Finished ", appLoc, outText, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
        append=TRUE)
}

# returnDates = Ch 1: only return dates of observations.
# zEWMAdates = Ch 1: get zewma values and zewma dates for every pixel (can be true at same time as writeProbs)
# writeProbs = Ch 1 and Ch 2: save original or Bayesian probability files
# bayesUpdate = Ch 2: do bayes updating with prior probs and save updated prob files
# landscapeVal = Ch 3: return only zEWMA values (no probs)
applyLandProbs(nRun=args[1], appRegion=5, bayesUpdate=TRUE, writeProbs=TRUE, 
               returnDates=FALSE, landscapeVal=FALSE, zEWMAdates=FALSE, nCores=20)
