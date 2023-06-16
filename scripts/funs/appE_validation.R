##########################################################
## Purpose: Functions for running validation tests for the main algorithm
##          based on the training data
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Dec 2022
## Last modified: Jan 2023
##########################################################
## -----------------------------------------------------------------------------
## getSensorDates = identify observation dates from both sensors
## formatCoords = format coordinates to show both UTM and latlon x-y
## getCellCoords = get focal pixel coordinates from landscape application output
## getValidationPix = full function to create table of validation pixels and
##                    dates to look over
## -----------------------------------------------------------------------------
getSensorDates <- function(nRun, appRegion){
  applicationStep <- 2
  nRun <- as.numeric(nRun)
  
  # Load variables
  ## nRun is used here to load data
  source("dissertation/scripts/args.R")
  source("dissertation/scripts/argsApp.R", local=TRUE)
  
  # get dates
  eachPixel <- lapply(1:nrow(allMatrices$landsat8sr$ts), reformat, 
                      dat=allMatrices)
  dateSens <- lapply(allMatrices, function(X) return(X$dates))
  datesRange <- as.vector(sapply(dateSens, function(X) return(range(X))))
  dates <- seq(min(datesRange), max(datesRange), by=1)
  
  stablePeriod <- c(dates[1], endStable)
  sensorNames <- names(eachPixel[[1]])
  
  out <- processPix(eachPixel[[166000]], sensorNames, dateSens, runType, window, 
                    spikeWin, spikeThresh, spikeAmp, spikeTime, stablePeriod,
                    backfill, lambda, logMod, returnStats, returnHarMod,
                    saveStats, combined, makeSumPlots, dates, windowType, 
                    inflation, staticInflation, spike2, getSensDates=TRUE)
  
  return(out$date)
}
formatCoords <- function(Y, oneDay, plotCRS){
  coord <- xyFromCell(oneDay, Y)
  ptVec <- vect(coord, crs=plotCRS)
  ptVecLL <- project(ptVec, "EPSG:4326")
  coordLL <- crds(ptVecLL)
  
  out <- data.table(cbind(coord, coordLL))
  colnames(out) <- c("utmX", "utmY", "degX", "degY")
  out[, cell := Y]
  return(out)
}
getCellCoords <- function(X, dates, runsDT, nCol, nCuts, probThresh, pixIgnore,
                          numValidate, appRegion, appLocations, nMatRow, 
                          plotCRS, dirPath, validCells){
  set.seed(34586)
  thisDate <- as.numeric(gsub("X", "", X))
  
  ## get day within full ts and read in data
  dayFile <- paste0("day", which(dates==thisDate))
  
  cellsAboveThresh <- calcDailyMetric(i=dayFile, dirPath, runsDT, nCol, nCuts, 
                                      probThresh, pixIgnore, ratio=FALSE, 
                                      index=TRUE, probs=TRUE)
  dailyDat <- as.data.table(do.call(cbind, cellsAboveThresh))
  dailyDat <- dailyDat[index %in% validCells[[X]], ]
  sampN <- sample(1:nrow(dailyDat), numValidate)
  focalCells <- dailyDat[sampN, index]
  focalProbs <- dailyDat[sampN, probs]/10000
  
  # create map so we can get centroids
  region <- appLocations[[appRegion]]
  extVec <- region$ext
  oneDay <- terra::rast(nrows=nMatRow, ncols=nCol, 
                        xmin=extVec[[1]], xmax=extVec[[2]], 
                        ymin=extVec[[3]], ymax=extVec[[4]],
                        crs=plotCRS)
  
  ## extract coords and also reproject to degree
  cellTab <- lapply(focalCells, formatCoords, oneDay, plotCRS)
  cellTab <- rbindlist(cellTab)
  cellTab[, `:=` (date = as.Date(thisDate, origin=as.Date("1970-01-01")),
                  dist = 1, #everything is supposed to be disturbed
                  probsRS = focalProbs)]
  setcolorder(cellTab, c("date", "cell", "utmX", "utmY", "degX", "degY", 
                         "probsRS"))
  return(cellTab)
}
getValidationPix <- function(appRegion, numValidate){
  applicationStep <- 3
  hpc <- TRUE
  bayes <- FALSE
  source("dissertation/scripts/argsApp.R", local=TRUE)
  source("dissertation/scripts/args.R", local=TRUE)
  source("dissertation/scripts/funs/appD_probMaps.R")
  
  # as a reminder, when you assign values to a raster, it starts at top left, 
  ## and goes left-right, top-bottom
  
  # Define timesteps and cells
  f <- list.files(paste0(saveDir, appLoc, "/validation"), full.names=TRUE)
  fShort <- gsub("^.*validation/", "", f)
  fNum <- gsub("dates_start", "", gsub("end.*", "", fShort))
  f <- f[order(as.numeric(fNum))]
  
  runs <- lapply(f, function(X){load(X); return(cellsPerDate)})
  test <- sapply(runs, function(X){
    return(as.numeric(gsub("date=", "", names(X))))
  })
  
  ## get the dates that are common to all parts of the region
  dates <- sort(Reduce(intersect, test))
  
  ## randomly select 5 different dates
  set.seed(203498)
  grp <- seq(1, length(dates), length.out=6)
  focalDates <- sapply(1:5, function(X){
    if(X < 5) out <- sample(dates[grp[X]:(grp[X+1]-1)], 1)
    if(X==5) out <- sample(dates[grp[X]:grp[X+1]], 1)
    return(out)
  })
  
  # Here we are masking out any pixels that don't have >= perCover % forest as
  ## of the beginning of 2019 (based on Copernicus 100m land cover product)
  maskRegion <- createForestMask(saveDir, appRegion, nMatRow, nCol, 
                                 appLocations, perCover=30)
  pixIgnore <- which(is.na(values(maskRegion)))
  
  ## get the cells that correspond to these dates and mask out pixIgnore
  validCells <- sapply(focalDates, function(d){
    cellsByDate <- sapply(1:length(runs), function(X){
      subVec <- runs[[X]][[paste0("date=", d)]]
      
      ## get the absolute pixel number from the overal raster. Note that we can
      ## definitively say 60 here because even for run 36, we're focused on 
      ## the previous rows.
      subVec <- subVec + ((60*nCol) * (X-1))
      return(subVec)
    })
    
    cellsByDate <- unlist(cellsByDate)
    cellsByDate <- cellsByDate[!(cellsByDate %in% pixIgnore)]
    return(cellsByDate)
  })
  
  names(validCells) <- paste0("X", focalDates)
  
  ## now, need to load in the probs and sample pixels based on those
  runsDT <- as.data.table(matRuns)
  nCuts <- runsDT[, nR][1]
  dirPath <- paste0(saveDir, probsDir)
  dates <- getDates(saveDir, fold=paste0(appLoc, "/landscape/"))
  
  # Loop over the dates and randomly select 30 pixels that are >=50% prob
  ## Get centroids and their probabilities, and put them into table
  out <- lapply(names(validCells), getCellCoords, dates, runsDT, nCol, nCuts, 
                probThresh, pixIgnore, numValidate, appRegion, appLocations, 
                nMatRow, plotCRS, dirPath, validCells)
  out <- rbindlist(out)
  out[, `:=` (region = appRegion, val="", dateImg="", notes="")]
  
  return(out)
}
