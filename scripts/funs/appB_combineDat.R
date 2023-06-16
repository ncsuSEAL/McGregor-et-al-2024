##########################################################
## Purpose: Functions for first step of landscape application, namely 
##          getting VRT data and applying QA/QC / ndvi calculation
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last updated: Nov 2021
##########################################################
## -----------------------------------------------------------------------------
## Get weights and apply to Z-scores
##
## reorderVRTs = reorder the VRTs to match sensor names
## createDataVectors = make initial vectors for running getRasterLines
## getRasterLines = main function, extract the data from VRTs and put in matrix
## getMetaDates = get metadata for an image 
## getMetaS1 = get metadata for Sentinel-1
## processMatrix = extract VRT data, filter and calculate index, and get 
##                image metadata

## -----------------------------------------------------------------------------
## Notes:
## Function `getRasterLines` 
##    This function is run for each band per sensor we want to extract. So for 
##    ndvi, we are running this 5 times for L8 (2 bands + 3 QA bands), 3 times 
##    for S2 (SCL + 2 bands), and 1 time for S1 (VH or VV).
##    The "dsets" argument refers to each VRT file (1 per date) in a full 
##    timeseries (per sensor)
## -----------------------------------------------------------------------------

reorderVRTs <- function(vrtList=NULL, sensor="", element=NA){
  dates <- sapply(vrtList[[sensor]], function(X) strsplit(X, "_")[[1]][element])
  dates <- gsub("T*", "", dates)
  dates <- as.Date(dates, format="%Y%m%d")
  correctedVRTs <- vrtList[[sensor]][order(dates)]
  return(correctedVRTs)
}
createDataVectors <- function(base, saveDir, sensors, appLoc){
  vrtFiles <- list.files(paste0(saveDir, appLoc, "/sensorVRT"), 
                         full.names = TRUE)
  
  vrtsInit <- list("landsat8sr" = vrtFiles[grepl("LC08.*_10mproj", vrtFiles)], 
                   "sentinel2l2a" = vrtFiles[grepl("S2A.*_10mproj", vrtFiles)], 
                   "sentinel1" = "")
  elements <- list("landsat8sr" = 5, "sentinel2l2a" = 2, "sentinel1" = 6)
  
  vrts <- lapply(1:3, function(Q){
    tt <- reorderVRTs(vrtList=vrtsInit, sensor=names(vrtsInit)[Q], 
                      element=elements[[Q]])
    return(tt)
  })
  names(vrts) <- names(vrtsInit)
  
  # these come from the .txt files used to create the VRTs
  ## note the S1 bands come from GEE, not the actual VRTs (can't rename the 
  ## bands without reading in the raster)).
  ## Also, with the QAQC bands in a different vrt, we start over at 1.
  bandList <- list(list("B4"=1, "B5"=2, "QA_PIXEL"=1, "QA_RADSAT"=2, 
                        "QA_AEROSOL"=3), 
                   list("NDVI"=1, "SCL"=1), 
                   list("VH"=2)) #("VV" = 1, "VH" = 2) 
  
  names(bandList) <- names(vrts)
  
  if(sensors=="L8S2"){
    return(list(bandList=bandList[c("landsat8sr", "sentinel2l2a")], 
                vrts=vrts[c("landsat8sr", "sentinel2l2a")]))
  } else if(sensors=="All"){
    return(list(bandList=bandList, vrts=vrts))
  }
}
getRasterLines <- function(dsets, band_name, start_row, n, nMatRow, nCol, 
                           datesObs, max_open_datasets=2.75e3){
  # extracts n rows from all dsets starting at start_row
  # returns a matrix of the values. Assumes all dsets have
  # the same pixel dimensions (nrow & ncol).
  
  # construct sensor-specific band names
  if(grepl("LC08", dsets[1])){
    dsets_bands <- bandList[grepl("landsat8", names(bandList))][[1]]
  } else if(grepl("S2A", dsets[1])){
    dsets_bands <- bandList[grepl("l2a", names(bandList))][[1]]
  } else {
    dsets_bands <- bandList[grepl("S1", names(bandList))][[1]]
  }
  
  ## subset to only the numeric or QA vrts
  if(grepl("QA|SCL", band_name)){
    dsets <- dsets[grepl("_qa", dsets)]
    dsets_bands <- dsets_bands[grepl("QA|SCL", names(dsets_bands))]
  } else {
    dsets <- dsets[grepl("_numeric", dsets)]
    dsets_bands <- dsets_bands[!grepl("QA|SCL", names(dsets_bands))]
  }
  
  dsets <- dsets[order(datesObs)]
  
  # get rows/cols from the first non-NA SDS
  # tmp_r <- raster(sds_names[which(!is.na(sds_names))[1]])
  nrows <- nMatRow
  ncols <- nCol
  # rm(tmp_r)
  
  # initialize output matrix
  out_vals <- matrix(NA, nrow=ncols * as.numeric(n), ncol=length(dsets))
  
  # determine how many blocks of datasets to open simultaneously
  ## this equates to 1, but leaving it in here for legacy
  num_blocks <- ceiling(length(dsets) / max_open_datasets)
  
  # loop through all datasets in a block and extract data
  for(dset_block in 1:num_blocks){
    # determine which dataset to start and end on
    dset_start <- ((dset_block - 1) * max_open_datasets) + 1
    dset_end <- min((dset_block * max_open_datasets), length(dsets))
    
    # open all datasets
    gds <- lapply(dset_start:dset_end, function(Z){
      return(GDAL.open(dsets[Z]))
    })
    
    # extract lines from valid datasets, repeat NAs for invalid ones
    for (j in 1:length(gds)){
      print(paste0("Currently getting data for date ", j, " of ", length(gds)))
      if(!is.na(gds[j])){
        val <- getRasterData(gds[[j]], band=dsets_bands[[band_name]], 
                             offset = c((start_row - 1), 0), 
                             region.dim = c(n, ncols), as.is = FALSE)
        out_vals[, (dset_start + j - 1)] <- c(val)    
      }else{
        out_vals[, (dset_start + j - 1)] <- NA
      }
    }
    # close all valid files
    for(i in which(!is.na(gds))) GDAL.close(gds[[i]])
    rm(gds)
  }
  return(out_vals)
}
getMetaDates <- function(f, sat=""){
  #note the first strsplit needs to be 15 if running on PC normally but 13 for hpc
  if(grepl("landsat", sat)){
    d <- as.Date(strsplit(strsplit(f, "/")[[1]][6], "_")[[1]][[3]],
                 format="%Y%m%d")
  } else if(grepl("sentinel2", sat)){
    d <- as.Date(gsub("T.*", "", strsplit(strsplit(f, "/")[[1]][6], "_")[[1]][[3]]),
                 format="%Y%m%d")
  } else if(grepl("sentinel1", sat)){
    d <- as.Date(gsub("T.*", "", strsplit(strsplit(f, "/")[[1]][6], "_")[[1]][[5]]),
                 format="%Y%m%d")
  }
  return(d)
}
getMetaS1 <- function(satName="", s1Files=c()){
  # get S1 direction by matching the file names
  s1meta <- fread("dissertation/data/myanmar/s1metadata.csv")
  s1meta[, id := gsub("_0.*", "", `system:index`)]
  
  s1data <- data.table(file=s1Files)
  s1data[, id := as.vector(sapply(s1Files, function(s){
    gsub(".vrt", "", tail(strsplit(s, "/")[[1]], n=1))
  }))]
  
  setkeyv(s1data, "id")
  setkeyv(s1meta, "id")
  s1data <- s1data[s1meta]
  return(s1data)
}
processMatrix <- function(vrtChoice, vrts, bandList, startRow, 
                          numberRow, nMatRow, nCol, fileOut, funs, appLoc,
                          vegIndex){
  satName <- names(vrts[vrtChoice])
  vrtVec <- vrts[[satName]]
  
  # Get dates for the observations
  #numeric not important, just getting only 1 date per obs
  vrtVec <- vrtVec[grepl("numeric", vrtVec)] 
  datesObs <- as.vector(sapply(vrtVec, getMetaDates, sat=satName))
  
  if(grepl("landsat", satName)){
    bands <- names(bandList[[which(grepl("landsat", names(bandList)))]])
  }
  if(grepl("l2a", satName)){
    bands <- names(bandList[[which(grepl("l2a", names(bandList)))]])
  }
  if(grepl("sentinel1", satName)){
    bands <- names(bandList[[which(grepl("sentinel1", names(bandList)))]])
  }
  
  # have to do this for every band necessary per sensor
  cl <- makeCluster(detectCores(), setup_strategy = "sequential")
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(rgdal))
  clusterExport(cl, c("vrtChoice", "bands", "getRasterLines", "satName", "vrts", 
                      "bandList", "startRow", "numberRow", "nMatRow", "nCol",
                      "datesObs"),
                envir=environment())
  
  # Put all sensor data into matrix, where columns are the dates and rows are 
  ## each pixel.
  ## note this takes a bit due to creating the matrix, opening each VRT, 
  ## then actually getting the data.
  dat <- parLapply(cl, bands, function(Y){
    print(paste0("Retrieving data for ", satName, " for Band ", Y))
    
    # cols are dates (each VRT), rows are pixel numbers
    ## with parLapply and n=50, S2 = ~15 mins, L8 = 2 min, S1 = ~12 mins
    bigMat <- getRasterLines(dsets = vrts[[vrtChoice]], band_name = Y, 
                             start_row = startRow, n = numberRow,
                             nCol=nCol, nMatRow=nMatRow, datesObs=datesObs)
    # number 82 for L8 (2019-01-18) and number 178 for S2 (2019-01-16)
    # cells <- c(114313, 114314, 114315, 117085, 117086, 117087, 119857, 119858, 119859)
    # sapply(cells, function(X) return(bigMat[X, ]))

    return(bigMat)
  })
  stopCluster(cl)
  names(dat) <- bands
  
  # apply QA/QC to the bands individually, then calculate ndvi. This returns a
  ## matrix of ndvi
  if(!grepl("sentinel1", satName)){
    matrixIndex <- filterThenIndex(dat, sat=satName, bandNames=bands, vegIndex) 
  } else {
    matrixIndex <- dat[[bands]]
    rm(dat)
  }
  
  write(paste0("Finished extracting data for ", appLoc, " from ", 
               names(vrts)[vrtChoice]), file=fileOut, append=TRUE)

  # Finally, put the dates in order now that the matrices have been 
  ## done in order as well
  datesObs <- datesObs[order(datesObs)]
  
  if(grepl("sentinel1", satName)){
    s1data <- getMetaS1(sat=satName, s1Files=vrtsVec)
    return(list(ts= matrixIndex, dir = s1data[,orbitProperties_pass], 
                dates=datesObs))
  } else {
    return(list(ts = matrixIndex, dates=datesObs)) 
  }
}
