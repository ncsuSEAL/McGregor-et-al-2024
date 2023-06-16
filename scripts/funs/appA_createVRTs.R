##########################################################
## Purpose: Functions to create VRTs of sensor data to be used in landscape 
##          application
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: July 2022
##########################################################
## -----------------------------------------------------------------------------
## prepL8text = prep the text file to make L8 vrts
## prepS2text = prep the text file to make S2 vrts
## prepS1text = prep the text file to make S1 vrts
## pepareVRTtext = wrapper for `prepXXtext` functions
## cropTIFs = crop the landscape tifs to just the desired extent
## indVRT = create vrt for each individual image file
## createVRT = wrapper for `indVRT` based on numeric or QAQC bands
## prepLandscapeApp = wrapper for all previous functions
## -----------------------------------------------------------------------------
prepL8text <- function(fileNum, runPath){
  ## ONLY for when the landsat data was downloaded as full tiles from USGS
  ## Bulk Download Application
  
  path <- paste0(runPath, "l8c2Data/", f[fileNum])
  files <- paste0(path, "/", list.files(path, pattern="*.TIF"))
  files <- files[grepl(paste0(paste0(indBands, collapse="|"), 
                              "|SR_Q|QA_"), files)]
  
  # sort them so the bands are first, followed by QA
  files <- c(files[grepl("B", files)], files[grepl("Q", files)])
  
  txtFileName <- paste0("l8c2Data/", f[fileNum], ".txt")
  
  file.create(txtFileName)
  fileConn<-file(txtFileName)
  writeLines(files, fileConn)
  close(fileConn)
  
  return(paste0("Finished group ", Y))
}
prepS2text <- function(fileNum, f, runPath, s2Path, indBands){
  img <- paste0(runPath, s2Path, f[fileNum])
  txtFileName <- gsub("tif", "txt", img)
  
  file.create(txtFileName)
  fileConn<-file(txtFileName)
  writeLines(img, fileConn)
  close(fileConn)
  
  # insidePath <- paste0(path, "/GRANULE/", 
  #                      list.files(paste0(path, "/GRANULE")), "/IMG_DATA/")
  # 
  # # initially we were doing 10m but we aren't anymore to be consistent
  # ## across all bands
  # # files10m <- list.files(paste0(insidePath, "R10m"))
  # # files10m <- files10m[grepl("B", files10m)]
  # # files10m <- paste0("R10m/", files10m)
  # 
  # ## As of Nov 2021, we are only using the 20m bands to resample
  # files20m <- list.files(paste0(insidePath, "R20m"))
  # 
  # # filter out the 10m bands that are duplicated
  # # files20m <- files20m[!grepl(paste0(gsub("10m", "20m", files10m), 
  # #                                    collapse="|"), files20m)]
  # files20m <- files20m[grepl(paste0(paste0(indBands, collapse="|"), "|SCL"),
  #                            files20m)]
  # files20m <- paste0("R20m/", files20m)
  # # files <- c(files10m, files20m)
  # files <- paste0(insidePath, files20m)
  # 
  # txtFileName <- paste0(s2Path, f[fileNum], ".txt")
  # txtFileName <- gsub(".SAFE", "_SAFE", txtFileName)
  
  
  
  return(paste0("Finished group ", fileNum))
}
prepS1text <- function(fileNum, runPath){
  thisFile <- fileName[fileNum]
  
  inputFiles <- f[grepl(thisFile, f)]
  inputFiles <- paste0(runPath, inputFiles)
  
  txtFileName <- paste0("s1ARD/", thisFile, ".txt")
  
  file.create(txtFileName)
  fileConn<-file(txtFileName)
  writeLines(inputFiles, fileConn)
  close(fileConn)
  
  return(paste0("Finished group ", fileNum))
}
prepareVRTtext <- function(sensor, appLoc, s2tile, index, runPath){
  
  ## don't need to make txt files for L8 app because already clipped when
  ### downloaded from GEE
  # if(grepl("landsat8sr", sensor)){
  #   if(index=="ndvi") indBands <- c("B4", "B5")
  #   
  #   #134044 is the row/path combo that fully includes Chatthin
  #   f <- list.files(paste0("l8c2Data/", appLoc), pattern="LC08.*")
  #   
  #   if(appLoc=="chatthin"){
  #     f <- f[!grepl(".txt", f) & grepl("134044", f)]
  #   }
  #   sapply(1:length(f), prepL8text, runPath, appLoc)
  # }
  
  if(grepl("sentinel2sr", sensor)){
    if(index=="ndvi") indBands <- c("B04", "B8A")
    
    ## 46QGM is the tile that contains Chatthin
    s2Path <- paste0("s2Data/L2A/", s2tile, "/")
    f <- list.files(s2Path)
    f <- f[!grepl(".txt", f)]
    
    sapply(1:length(f), prepS2text, f, runPath, s2Path, indBands)
  }
  
  if(grepl("sentinel1", sensor)){
    f <- list.files("s1ARD", pattern="*.tif")
    fileName <- unique(as.vector(sapply(f, function(X){
      return(strsplit(X, "_.{6}_")[[1]][1])
    })))
    
    f <- paste0("s1ARD/", f)
    
    # create txt files with each group of tifs to be made into their own VRTs
    sapply(1:length(fileName), prepS1text, runPath)
  }
  
  print(paste0("All text files created for ", sensor, " for ", appLoc, "."))
}
cropTIFs <- function(X, sensor, filePath, extentCrop, epsg, outFolder, dryrun){
  inFolder <- filePath
  
  ## Replace file extension in prep
  if(sensor=="L8"){
    rootOut <- gsub(inFolder, paste0(inFolder, "/", outFolder), X)
  } else if(sensor=="S2"){
    rootOut <- gsub(paste0(strsplit(inFolder, "/")[[1]][3], "/"), 
                    paste0(outFolder, "/"), X)
  }
  
  ## use gdalwarp to both crop and reproject the tif if needed
  gdalwarp(srcfile=X, dstfile=rootOut, t_srs=paste0("EPSG:", epsg),
           te=extentCrop, overwrite=TRUE, dryrun=dryrun)
}
indVRT <- function(file, filePath, bandNum, appLoc, isS1, type, extentTarget, 
                   dryrun){
  # browser()
  inFolder <- filePath
  outFolder <- paste0("dissertation/data/myanmar/", appLoc, "/sensorVRT")
  
  ## Replace folder name and file extension in prep for vrt conversion
  vrtFileName <- gsub(inFolder, outFolder, file)
  vrtFileName <- gsub(".tif", paste0("_", type, ".vrt"), vrtFileName)
  
  # Build the initial VRTs, resample to 10m, choose resampling method
  gdalbuildvrt(gdalfile=file,  output.vrt=vrtFileName, b=bandNum,
               overwrite=TRUE, dryrun=dryrun)
  
  # Now reproject, resample to 10m, specify resampling method
  vrtProjName <- gsub(".vrt", "_10mproj.vrt", vrtFileName)
  r <- ifelse(type=="numeric", "cubic", "nearest")
  
  gdalwarp(srcfile=vrtFileName, dstfile=vrtProjName, t_srs="EPSG:32646", 
           te=extentTarget, r=r, tr=c(10,10), tap=TRUE,
           overwrite=TRUE, dryrun=dryrun)
  
  return(vrtFileName)
}
createVRT <- function(Y, sensor, fileList, filePath, extentTarget, appLoc,
                      S2Bands){
  rFile <- rast(fileList[Y])
  isS1 <- ifelse(sensor=="S1", TRUE, FALSE)
  
  ## need to rename the bands because S2 loses in the comglomerate process
  if(sensor=="S2") names(rFile) <- s2Bands
  
  ## numeric bands
  bandNum <- which(!grepl("QA|SCL", names(rFile)))
  vrtNum <- indVRT(file=fileList[Y], filePath, bandNum, appLoc, 
                   isS1, type="numeric", extentTarget, dryrun=FALSE)
  
  ## QAQC bands
  bandNum <- which(grepl("QA|SCL", names(rFile)))
  if(length(bandNum) > 0){
    vrtQA <- indVRT(file=fileList[Y], filePath, bandNum, appLoc, 
                    isS1, type="qa", extentTarget, dryrun=FALSE)
  }
  
  # gdalbuildvrt(gdalfile=c(vrtNum, vrtQA), output.vrt="sensorVRT/test.vrt",
  #              separate=TRUE)
  # gdalinfo("sensorVRT/test.vrt")
}
prepLandscapeApp <- function(exts, appLocations, cropTIF, makeVRT, parL, 
                             epsg, s2Bands){
  print(paste0("Start processing region ", names(appLocations)[exts]))
  
  # First bring in the extents of each landscape area
  
  ## NEED TO CHANGE EXTENT ORDER TO BE xmin ymin xmax ymax
  ### this is because the order in gdal is different
  extent <- appLocations[[exts]]$ext
  extentTarget <- c(extent[1], extent[3], extent[2], extent[4])
  s2tile <- appLocations[[exts]]$s2tile
  
  ## need to change the projection because some of the s2 tiles are in diff proj
  pts <- data.frame(lon=c(extent[1], extent[1], extent[2], extent[2]), 
                    lat=c(extent[3], extent[4], extent[4], extent[3]))
  shp <- vect(pts, crs=paste0("EPSG:", epsg))
  shp <- project(shp, appLocations[[exts]]$crs)
  
  extentOrig <- as.numeric(c(ext(shp)[1], ext(shp)[3], ext(shp)[2], ext(shp)[4]))
  extentCrop <- c(extentOrig[1]-1000, extentOrig[2]-1000, 
                  extentOrig[3]+1000, extentOrig[4]+1000)
  
  # Second, either crop tifs and/or make VRTs
  ## could make this parallel if want speed
  appLoc <- names(appLocations)[exts]
  substr(appLoc, 1, 1) <- toupper(substr(appLoc, 1, 1))
  appLoc <- paste0("app", appLoc)
  
  if(cropTIF){
    print("Starting to crop S2 tifs")
    sensor <- c("S2") #abbr to match folder names
    s2Path <- paste0("s2Data/L2A/", appLocations[[exts]]$s2tile)
    
    ## First, crop the tifs
    ### only need to do S2 for new app locations bc l8 already downloaded cropped
    filePath <- s2Path
    fileList <- list.files(filePath, pattern=".tif", full.names=TRUE)
    
    if(parL){
      cl <- parallel::makeCluster(detectCores()-1)
      clusterEvalQ(cl, library(data.table))
      clusterEvalQ(cl, library(gdalUtilities))
      clusterExport(cl, c("sensor", "filePath", "extentCrop", "appLoc",
                          "cropTIFs"), 
                    envir=environment())
      clusterExport(cl, c(origVars))
      parSapply(cl, X=fileList, cropTIFs, sensor, filePath, extentCrop, epsg,
                outFolder=appLoc, dryrun=FALSE)
      stopCluster(cl)
    } else {
      sapply(fileList, cropTIFs, sensor, filePath, extentCrop, epsg, 
             outFolder=appLoc, dryrun=FALSE)
    }
    
    print("Finished cropping S2 tifs")
    
    # sapply(sensType, function(sensor){
    #   filePath <- ifelse(sensor=="L8", "l8c2Data",
    #                      ifelse(sensor=="S2", s2Path, "s1ARD"))
    #   fileList <- list.files(filePath, pattern=".txt", full.names=TRUE)
    #   sapply(fileList, cropTIFs, sensor, filePath, extentCrop, dryrun=FALSE)
    # })
  }
  if(makeVRT){
    sensType <- c("L8", "S2") #abbr to match folder names
    sapply(sensType, function(sensor){
      print(paste("Starting to make aligned and resampled VRTs for", sensor))
      filePath <- ifelse(sensor=="L8", "l8c2Data",
                         ifelse(sensor=="S2", "s2Data/L2A", "s1ARD"))
      filePath <- paste0(filePath, "/", appLoc)
      fileList <- list.files(filePath, pattern=".tif", full.names=TRUE)
      
      if(parL){
        cl <- parallel::makeCluster(detectCores()-1)
        clusterEvalQ(cl, library(data.table))
        clusterEvalQ(cl, library(terra))
        clusterEvalQ(cl, library(gdalUtilities))
        clusterExport(cl, c("sensor", "fileList", "filePath", "extentTarget", 
                            "appLoc", "s2Bands"), envir=environment())
        clusterExport(cl, c(origVars))
        parSapply(X=1:length(fileList), createVRT, sensor, fileList, filePath,
                  extentTarget, appLoc, s2Bands, cl=cl)
        stopCluster(cl)
      } else {
        sapply(1:length(fileList), createVRT, sensor, fileList, filePath,
               extentTarget, appLoc)
      }
      print(paste("Finished making aligned and resampled VRTs for", sensor))
    })
  }
  print(paste0("Finish processing region ", names(appLocations)[exts]))
}
################################################################################
# Archive
cropTIFsArchive <- function(X, sensor, filePath, extentCrop, outFolder, dryrun){
  inFolder <- filePath
  text <- as.vector(unlist(fread(X, header=FALSE)))
  
  ## Replace file extension in prep
  
  if(sensor=="L8"){
    rootOut <- gsub(inFolder, paste0(inFolder, "/", outFolder), X)
  } else if(sensor=="S2"){
    rootOut <- gsub(paste0(strsplit(inFolder, "/")[[1]][3], "/"), 
                    paste0(outFolder, "/"), X)
  }
  
  vrtName <- gsub(".txt", ".vrt", rootOut)
  tifName <- gsub(".txt", ".tif", rootOut)
  
  ## Build a VRT cropped to CropExtent (no other transformation)
  sepBands <- ifelse(sensor=="S1", FALSE, TRUE)
  gdalbuildvrt(gdalfile=text, output.vrt=vrtName, separate=sepBands,
               te=extentCrop, overwrite=TRUE, dryrun=dryrun)
  
  ## Write out the cropped VRT to a tif
  ### we're using a common datatype of UInt16 bc regardless of QA or band val,
  ### numbers are > 0 (aka unsigned) and stored as integers
  ### L8 description: https://www.usgs.gov/media/files/landsat-8-9-collection-2-level-2-science-product-guide
  ### S2 description: https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/resolutions/radiometric
  gdal_translate(src_dataset=vrtName, dst_dataset=tifName, ot="UInt16", 
                 a_nodata="none", dryrun=dryrun)
}