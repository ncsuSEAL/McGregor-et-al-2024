##########################################################
## Purpose: Functions to extract Sentinel 2 data from downloaded and converted
##          L2A images using sen2cor software.
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 3.6.3, May 2021
## Last modified: R-4.1.1, July 2022
##########################################################
## -----------------------------------------------------------------------------
## Function = purpose
##
## extractBands = Function 0: Overall function to get Sentinel-2 band values
### To do this function, you MUST have already downloaded the necessary 
### Sentinel-2 L1C images and subsequently processed them to be L2A (file path 
### denoted by the L2APath variable). 
### This function creates a csv of Sentinel-2 L2A band values for each of the 
### training points, for each date in their timeseries. The output table mimics 
### the structure of the other tables that were directly obtained from this 
### research's Google Earth Engine scripts for Landsat-8 and Sentinel-1.

## prepData = Fun 1: Get list of images to look through for extraction.
### This function can be used to filter out images already extracted if doing 
### this in batches. Leaving the last 2 arguments as their defaults will include 
### all the imagery (referring to Functions 2A and 2B). Here, prevDT in this 
### case would be the output from an earlier run of the whole process, e.g.
### prev <- fread(file.path(base, "combTiles.csv"))

## noDups = Fun 1A: Filter out images already used to extract data
## specificDate = Fun 1B: Filter out images based on specific dates
## createS2DT = Fun 2: Create dt with extracted band values per training point
## getRastData = Fun 2A: extract data for points in a specific image
## -----------------------------------------------------------------------------
extractBands <- function(prev=NULL, tileStr, L2APath, 
                         txtFile, trainPath, tilePath, 
                         removeNArows=FALSE, outFile){
  if(tileStr=="all") tileStr <- list.files(file.path(L2APath))
  tileStr <- tileStr[!grepl("app", tileStr)]
  
  # Create a table showing all tiles for the training data
  dt <- fread(file.path(trainPath))
  pointTiles <- fread(file.path(tilePath))
  pointTiles <- pointTiles[tile %in% tileStr, ]
  
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(terra))
  clusterExport(cl, list("prev", "tileStr", "L2APath", "pointTiles", 
                         "txtFile", "prepData", "noDups", "specificDate", 
                         "createS2DT", "getRastData"),
                envir=environment())
  
  # Get list of images to run through for extraction
  imgs <- unlist(parLapply(cl, tileStr, prepData, prevDT=prev,
                           dateString1="", L2APath=L2APath))
  eachTile <- parLapply(cl, imgs, createS2DT, pointSub=pointTiles, 
                        txtOut=txtFile, pathL2A=L2APath)
  # eachTile <- lapply(imgs, createS2DT, pointSub=pointTiles, 
  #                       txtOut=txtFile, pathL2A=L2APath)
  stopCluster(cl)
  combTiles <- rbindlist(eachTile)
  
  # Do some data clean-up before saving
  ## some dates have no data, find and replace with NAs 
  ### (but keep the row for now in case need to compare with TOA later)
  colBands <- colnames(combTiles)[!grepl("satellite|date|lat|lon|pointid",
                                         colnames(combTiles))]
  rowN <- c(as.numeric(which(base::rowSums(combTiles[,lapply(.SD, is.na),
                                                     .SDcols=colBands]) == 1)))
  
  if(removeNArows & length(rowN) > 0){
    combTiles <- combTiles[-rowN, ]
  } 
  # else {
  #   combTiles <- combTiles[rowN, c(colBands) := NA]
  # }
  
  ## rbind to any previous data
  fullTab <- rbind(prev, combTiles)
  
  ## THIS WRITES TO THE HPC FOLDER
  fwrite(fullTab, outFile)
  write(paste0("Finished all extractions at: ", Sys.time()), file=txtFile, 
        append = TRUE)
}
prepData <- function(Tile, prevDT=NULL, dateString1="", 
                     L2APath=""){
  #filter out imgs already extracted
  imgsFoc <- noDups(tile=Tile, prevDT, L2APath)
  
  #only look at imgs from certain date (e.g. only 2019 images)
  imgsFilt <- specificDate(imgsFoc, dateString=dateString1)
  return(imgsFilt)
}

noDups <- function(tile, prevDT, L2APath){
  ## get folder and files of L2A images
  imgs <- list.files(file.path(L2APath, tile), full.names = TRUE, pattern=".tif")
  if(!is.null(prevDT)){
    stop("Need to fix this part of `noDups` function")
    satID <- paste0("T", tile, "_", substring(basename(imgs),12,26))
    prevID <- gsub("S2_L2A/", "", unique(prevDT[,satellite]))
    
    imgs <- imgs[!satID %in% prevID]
  }
  
  if(length(imgs)!=0){
    return(imgs)
  }
}

specificDate <- function(imgsFoc, dateString){
  filt <- imgsFoc[grepl(dateString, imgsFoc)]
  return(filt)
}

createS2DT <- function(Q, pointSub=NULL, txtOut="", pathL2A=""){
  # browser()
  # print(paste0("Working on img ", Q))
  tileNm <- strsplit(Q, "/")[[1]][[3]]
  
  # Create a spatial vector (shapefile) for the training points in this tile
  ## Note normally this would be done outside the function because it's the 
  ## same for each tile, but this isn't working for parallel
  tileDT <- pointSub[tile == tileNm, .(lon, lat, pointid)]
  tileVec <- vect(as.data.frame(tileDT), crs="+init=epsg:4326")
  
  # Write out progress statement
  write(paste0("Started extracting bands for tile ", tileNm, ", Image ", 
                 Q, " at: ", Sys.time()), 
          file=txtOut, append = TRUE)

  # extract the data for the points in this specific tile
  allBands <- getRastData(img=Q, tileVec=tileVec) 
  allBandsDT <- cbind(data.table(pointid=tileDT$pointid), allBands)

    ## combine and format column names (add satellite and pointid)
    sat <- gsub(".tif", "", strsplit(Q, "/")[[1]][4])
    satDate <- as.Date(as.POSIXct(substring(sat, 12,19), format="%Y%m%d"))
    
    allBandsDT <- allBandsDT[, `:=` (satellite = sat, date = satDate,
                                     lat=tileDT$lat, lon=tileDT$lon)]
    
    # Reorder columns
    setcolorder(allBandsDT, c("satellite", "date", "lat", "lon", "pointid", 
                              "NDVI", "SCL"))
    
    write(paste0("Finished extracting bands for tile ", tileNm, ", Image ",Q, " at: ", 
                  Sys.time()), file=txtOut, append = TRUE)
  
    return(allBandsDT)
  }

getRastData <- function(img, tileVec){
  # browser()
  ## read in the tif
  r <- rast(img)
  ## only keep the points that lie within the tile
  # if(is.na(crs(r))) crs(r) <- "+init=epsg:32646"
  if(crs(tileVec) != crs(r)) tileVec <- project(tileVec, r)
  extVals <- as.data.table(extract(r, tileVec))
  extVals[, ID := NULL]
  colnames(extVals) <- c("NDVI", "SCL")
  
  return(extVals)
}

# createS2DT <- function(Q, pointSub=NULL, txtOut="", pathL2A=""){
#   # browser()
  
#   # check to see if there are missing images
#   nFolders <- length(list.files(paste0(Q, "/IMG_DATA")))
#   nImgs10 <- length(list.files(paste0(Q, "/IMG_DATA/R10m")))
#   nImgs20 <- length(list.files(paste0(Q, "/IMG_DATA/R20m")))
#   imgSub <- strsplit(Q, "L2A/")[[1]][2]
#   tileNm <- substr(imgSub, 0, 5)
#   imgNm <- strsplit(strsplit(imgSub, paste0(tileNm, "/"))[[1]][2], ".SAFE/")[[1]][1]
  
#   num <- list.files(paste0(pathL2A, "/", tileNm))
#   totalFiles <- length(num)
#   thisImg <- which(num==paste0(imgNm, ".SAFE"))
  
#   # Create a spatial vector (shapefile) for the training points in this tile
#   ## Note normally this would be done outside the function because it's the 
#   ## same for each tile, but this isn't working for parallel
#   tileVec <- pointSub[tile == tileNm, .(lon, lat, pointid)]
#   tileVec <- vect(as.data.frame(tileVec), crs="+init=epsg:4326")
  
#   # Write out progress statement
#   if(nFolders<2 | nImgs10==0 | nImgs20==0){
#     write(paste0("Error for tile ", tileNm, ", Image ", thisImg, " of ", 
#                  totalFiles, " due to missing band TIFs: ", imgNm), 
#           file=txtOut, append = TRUE)
#   } else {
#     write(paste0("Started extracting bands for tile ", tileNm, ", Image ", 
#                  thisImg, " of ", totalFiles, " at: ", Sys.time(), ", ", imgNm), 
#           file=txtOut, append = TRUE)
    
#     # Start extracting the bands
#     # define the file names for the Sentinel-2 images
#     bands <- c(10, 20)
#     bandFolder <- c("R10m", "R20m")
#     bandsList <- list(bands10=c("B02,B03,B04,B08,AOT,WVP"), 
#                       bands20=c("B05,B06,B07,B8A,B11,B12,SCL"))
    
#     # extract the data for the points in this specific tile
#     allBands <- lapply(1:2, getRastData, eachImg=Q, bandFolder=bandFolder, bandsList=bandsList, tileVec=tileVec)
    
#     ## combine and format column names (add satellite and pointid)
#     allBandsDT <- do.call(cbind, allBands)
#     sat <- paste0("S2_L2A/", gsub(".{8}$", "", colnames(allBandsDT)[2]))
#     satDate <- as.Date(as.POSIXct(substring(sat, 15), format="%Y%m%dT%H%M%S"))
#     namesCol <- colnames(allBandsDT)[2:ncol(allBandsDT)]
#     bandNames <- namesCol[!grepl("MSK", namesCol)]
#     bandNames <- gsub("0", "", substring(bandNames, nchar(bandNames)-6, nchar(bandNames)-4))
    
#     # there is no MSK_CLDPRB or MSK_SNWPRB when I only do 10 and 20m outputs
#     # maskNames <- namesCol[grepl("MSK", namesCol)]
#     # maskNames <- substring(maskNames, nchar(maskNames)-13, nchar(maskNames)-4)
#     colnames(allBandsDT) <- c("pointid", bandNames)
#     allBandsDT <- allBandsDT[, `:=` (satellite = sat, date = satDate,
#                                      lat=tileVec$lat, lon=tileVec$lon)]
    
#     # Reorder columns
#     setcolorder(allBandsDT, c("satellite", "date", "lat", "lon", "pointid", 
#                               "B2", "B3", "B4", "B5", "B6", "B7", "B8",
#                               "B8A", "B11", "B12", "AOT", "WVP", "SCL"))
    
#     write(paste0("Finished extracting bands for tile ", tileNm, ", Image ",thisImg, " of ", totalFiles, " at: ", Sys.time(), ", ", imgNm), file=txtOut, append = TRUE)
#     return(allBandsDT)
#   }
# }
# getRastData <- function(resImg, eachImg="", bandFolder="", bandsList=NA, 
#                         tileVec=NA){
#   # browser()
#   ## get the necessary files for this resolution (to match GEE)
#   imgPath <- file.path(eachImg, "IMG_DATA", bandFolder[resImg])
#   files <- list.files(imgPath, full.names = TRUE, pattern=".jp2")
#   files <- files[grepl(gsub(",", "|", bandsList[resImg]), files)] #remove unwanted files
  
#   # this is for if doing 60m bands
#   # else {
#   # imgPath <- file.path(tileName, "QI_DATA")
#   #     files <- list.files(imgPath, full.names = TRUE, pattern=".jp2")
#   #     files <- files[grepl("MSK", files)]
#   #     files <- files[!grepl("6", substring(files, nchar(files)-6))]
#   # }
  
#   ## read in the files as a raster stack
#   r <- rast(files)
#   ## only keep the points that lie within the tile
#   # if(is.na(crs(r))) crs(r) <- "+init=epsg:32646"
#   if(crs(tileVec) != crs(r)) tileVec <- project(tileVec, r)
#   extVals <- as.data.table(extract(r, tileVec))
#   extVals[, ID := NULL]
  
#   if(resImg==1){
#     ids <- data.table(pointid=tileVec$pointid)
#     extVals <- cbind(ids, extVals)
#   }
#   return(extVals)
# }

# noDups <- function(tile="46QHL", prevDT=prev, filePathL2A=""){
#   ## get folder and files of L2A images
#   imgs <- list.files(file.path(filePathL2A, tile), full.names = TRUE)
#   if(!is.null(prevDT)){
#     satID <- paste0("T", tile, "_", substring(basename(imgs),12,26))
#     prevID <- gsub("S2_L2A/", "", unique(prevDT[,satellite]))
    
#     imgs <- imgs[!satID %in% prevID]
#   }
  
#   if(length(imgs)!=0){
#     imgsSpec <- list.files(file.path(imgs, "GRANULE"), full.names=TRUE)
#     return(imgsSpec)
#   }
# }

