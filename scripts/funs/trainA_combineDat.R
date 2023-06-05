##########################################################
## Purpose: Combine the satellite data into matrices similar to 
##          landscape application
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: April 2022
##########################################################
## -----------------------------------------------------------------------------
## Function = purpose
##
## expandValues = get full ts of data, substituting NAs for days without obs
## combineBandMat = put all band values into a single matrix
## getIndex = create matrix of index values (e.g. NDVI) per sensor with each 
##            column representing a date
## combineLandsat = combine L8 and L9 data for Schenck
## -----------------------------------------------------------------------------
combineLandsat <- function(dataPath){
  f <- list.files(file.path(dataPath, "sensorData"), pattern="landsat", full.names=TRUE)
  dt <- rbindlist(lapply(f, function(X) fread(X)))
  dt <- dt[order(pointid, date), ]
  fwrite(dt, file.path(dataPath, "sensorData", "landsat8sr.csv"))
}
expandValues <- function(dt, col, dates){
  
  # Get a full time series and match the dates to the actual observation dates
  ## If duplicate dates, retain only the mean of each band. We're doing this because
  ### the differences between dup obs per sensor per day are solely due to a difference
  ### in reprojection over different tiles, so the values should be consistent.
  ### Also because otherwise we'd have to retain date vectors from every single pixel.
  #### This is primarily for Sentinel-2.
  
  dt <- dt[, .(date, get(col))]
  setnames(dt, old=c("V2"), new="focusCol")
  if(col!="direction"){
    dt <- dt[,list(focusCol= mean(focusCol)), by=c("date")]
  } else {
    
    # Duplicated dates for S1 have the same direction, so it's fine to just remove
    ## one of the instances.
    dt <- dt[!duplicated(date), ]
  }
  
  expandVals <- rep(NA, length(dates))
  expandVals[dates %in% as.numeric(dt$date)] <- 
    dt[as.numeric(date) >= min(dates), focusCol]
  return(expandVals)
}
combineBandMat <- function(Y, dt, dates){
  bandCol <- colnames(dt)[grepl(Y, colnames(dt))]
  out <- expandValues(dt=dt, col=bandCol, dates=dates)
  return(out)
}
getIndex <- function(satName, dataPath, sensorFold, sensorFiles, 
                     point, bandList, dates, vegIndex){
  # Read in the sensor data
  dt <- fread(file.path(dataPath, sensorFold, 
                        sensorFiles[grepl(satName, sensorFiles)]))
  dt <- dt[, date := as.Date(date)][pointid==point, ][order(date), ]
  
  # Get the band names necessary for calculating the index and filtering on QAQC
  bands <- names(bandList[[which(grepl(satName, names(bandList)))]])
  
  # Remove rows that don't have data
  ## here just using the first band as representative
  bName <- colnames(dt)[grepl(bands[[1]], colnames(dt))]
  rowsRemove <- which(is.na(dt[, get(bName)]))
  if(length(rowsRemove != 0)) dt <- dt[-rowsRemove, ]
  
  # Create matrix of values for training point for each band across a full ts
  dat <- lapply(bands, combineBandMat, dt=dt, dates=dates)
  names(dat) <- bands
  
  # Filter data by QAQC and calculate index
  if(!grepl("sentinel1", satName)){
    matrixIndex <- filterThenIndex(dat, sat=satName, bandNames=bands, vegIndex)
    return(list(ts=matrixIndex))
  } else {
    matrixIndex <- dat[[bands]]
    dir <- expandValues(dt=dt, col="direction", dates=dates)
    return(list(ts=matrixIndex, dir=dir))
  }
  
  if(max(matrixIndex$ts, na.rm=TRUE) > 1 | min(matrixIndex, na.rm=TRUE) < 0){
    stop(paste0("NDVI values are outside valid range for ", point))
  } else {
    return(matrixIndex)
  }
}
