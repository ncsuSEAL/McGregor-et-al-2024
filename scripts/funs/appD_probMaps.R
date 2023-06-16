##########################################################
## Purpose: Functions for third step of landscape application, specifically the
##          creation of daily probability maps and calculating the fraction of
##          the landscape for each day that is disturbed.
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: July 2022
##########################################################
## -----------------------------------------------------------------------------
## getDates = return the specific dates of the landscape timeseries
## calcRatio = calculate ratio of disturbed pixels for a specific day
## dailyDistRatio = run calcRatio over all days of the timeseries
## getDailyData = load the binary probability data for specific day
## createProbMap = create and save probability map for specific day
## makeDailyImgs = wrapper to iterate `createProbMap` for every day of ts
## -----------------------------------------------------------------------------
getDates <- function(saveDir, fold, obs=FALSE){
  f <- list.files(paste0(saveDir, fold))
  allMatrices <- readRDS(paste0(saveDir, fold, 
                                f[grepl("mats1-", f)]))
  dateSens <- lapply(allMatrices, function(X) return(X$dates))
  
  if(obs){
    dates <- sort(unique(unlist(dateSens)))
  } else {
    datesRange <- as.vector(sapply(dateSens, function(X) return(range(X))))
    dates <- seq(min(datesRange), max(datesRange), by=1)
  }
  return(dates)
}
createForestMaskCopernicus <- function(path, perCover, appLocations, appRegion){
  forest19 <- rast(paste0(path, "forestCover/forest19.tif"))
  
  # reclassify to only keep >= percForest
  forest19 <- forest19$`tree-coverfraction`
  forest19[forest19 < perCover] <- NA
  forest19[!is.na(forest19)] <- 1
  
  # create a polygon shapefile to reproject to, then crop original data
  ex <- appLocations[[appRegion]]$ext
  
  polyg <- matrix(rbind(c(ex[1], ex[3]), c(ex[1], ex[4]), 
                        c(ex[2], ex[4]), c(ex[2], ex[3])), ncol=2)
  
  v <- vect(polyg, crs=appLocations[[appRegion]]$crs, type="polygons")
  forest19 <- project(forest19, appLocations[[appRegion]]$crs)
  out <- crop(forest19, v)
  
  # resample to 10m res so we don't mask out pixels that are fine
  fac <- floor(res(out)[1] / 10)
  out <- disagg(out, fac, method="near") 
  return(out)
}
createForestMask <- function(path, appRegion, nMatRow, nCol, appLocations, 
                             perCover, returnGlad=FALSE){
  ## data comes from https://storage.googleapis.com/earthenginepartners-hansen/GFC-2021-v1.9/download.htmlDiss
  
  # bring in the GLAD tif
  loss <- rast(paste0(path, "Hansen_GFC-2021-v1.9_lossyear_30N_090E.tif"))
  cover <- rast(paste0(path, "Hansen_GFC-2021-v1.9_treecover2000_30N_090E.tif"))
  
  # create a polygon shapefile to reproject to, then crop original data
  ex <- appLocations[[appRegion]]$ext
  
  polyg <- matrix(rbind(c(ex[1], ex[3]), c(ex[1], ex[4]), 
                        c(ex[2], ex[4]), c(ex[2], ex[3])), ncol=2)
  
  v <- vect(polyg, crs=appLocations[[appRegion]]$crs, type="polygons")
  v <- project(v, loss)
  
  lossCrop <- crop(loss, ext(v))
  coverCrop <- crop(cover, ext(v))
  
  # only keep pixels with >= perCover forest cover
  coverCrop[coverCrop < perCover] <- NA
  
  # only keep pixels that were either never disturbed or disturbed after 2018
  lossCrop <- lossCrop == 0 | lossCrop > 18
  
  # mask the disturbed label with only those pixels that were forest
  masked <- mask(lossCrop, coverCrop, maskvalues=NA)
  masked[masked == 0] <- NA
  
  target <- terra::rast(nrows=nMatRow, ncols=nCol, 
                        xmin=ex[[1]], xmax=ex[[2]], 
                        ymin=ex[[3]], ymax=ex[[4]],
                        crs=appLocations[[appRegion]]$crs)
  targetLL <- project(target, lossCrop)
  
  masked <- resample(masked, targetLL, method="near")
  
  # reproject to UTM and resample to 10m
  
  if(returnGlad){
    lossCrop[lossCrop==1] <- NA
    out <- project(lossCrop, appLocations[[appRegion]]$crs)
    out <- resample(out, target, method="near")
  } else {
    out <- project(masked, target)
  }
  return(out)
}
calcDailyMetric <- function(i, dirPath, runsDT, nCol, nCuts, probThresh, 
                            pixIgnore, ratio, index, probs=FALSE){
  specDay <- list.files(file.path(dirPath, i))
  # print(i)
  
  ## Read in the data
  allData <- sapply(1:nrow(runsDT), getDailyData, runsDT, specDay, 
                    dirPath, dayF=i, nData=nCol*nCuts)
  
  ## Combine and rescale
  vecData <- unlist(allData)
  vecData[pixIgnore] <- NA
  
  if(ratio){
    vecData <- vecData[!is.na(vecData)]
    
    totalPix <- length(vecData)
    dailyRatio <- length(vecData[vecData >= (probThresh*10000)]) / totalPix
    return(dailyRatio)
  }
  
  if(index){
    distCells <- which(vecData >= (probThresh*10000))
    
    if(probs){
      distCells <- list(index=distCells, 
                        probs=vecData[which(vecData >= (probThresh*10000))])
    }
    
    return(distCells)
  }
}
# 
# firstDistDate <- c(rep(NA, max(matRuns$endR)*nCol))
# distIndex <- which(vecData >= probThresh*10000)
# sapply(distIndex, function(idx){
#   actualDate <- dates[as.numeric(gsub("day", "", i))]
#   if(is.na(firstDistDate[idx]) | actualDate < firstDistDate[idx]){
#     firstDistDate[idx] <- actualDate
#   }
# })

dailyDistMetric <- function(dates, dirPath, matRuns, nCol, nCuts, probThresh,
                           pixIgnore, ratio, index){
  plotDates <- seq(dates[1], dates[length(dates)], by=1)
  
  dayFiles <- list.files(dirPath)
  dayFiles <- dayFiles[order(as.numeric(gsub("day", "", dayFiles)))]
  runsDT <- as.data.table(matRuns)
  nCuts <- runsDT[, nR][1]
  
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("calcDailyMetric","getDailyData", "dirPath", "runsDT", 
                      "nCol", "nCuts", "probThresh", "pixIgnore", "ratio",
                      "index"), 
                envir=environment())
  distMetric <- parSapply(dayFiles, calcDailyMetric, dirPath, runsDT, nCol, 
                          nCuts, probThresh, pixIgnore, ratio, index, cl=cl)
  stopCluster(cl)
  
  if(index) distMetric <- sort(unique(unlist(distMetric)))
  return(distMetric)
}
dailyDistMonitor <- function(dates, dirPath, matRuns, nCol, nCuts, val){
  plotDates <- seq(dates[1], dates[length(dates)], by=1)
  dayFiles <- list.files(dirPath)
  dayFiles <- dayFiles[order(as.numeric(gsub("day", "", dayFiles)))]
  
  dayFiles <- dayFiles[plotDates >= as.numeric(as.Date("2019-01-01"))]
  plotDates <- plotDates[plotDates >= as.numeric(as.Date("2019-01-01"))]
  
  runsDT <- as.data.table(matRuns)
  nCuts <- runsDT[, nR][1]
  
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("calcDailyMetric","getDailyData", "dirPath", "runsDT", 
                      "nCol", "nCuts", "val"), 
                envir=environment())
  probsMonitor <- parSapply(cl, dayFiles, function(i){
    specDay <- list.files(file.path(dirPath, i))
    
    ## Read in the data
    allData <- sapply(1:nrow(runsDT), getDailyData, runsDT, specDay, 
                      dirPath, dayF=i, nData=nCol*nCuts)
    
    ## Combine and rescale
    vecData <- unlist(allData)
    
    out <- vecData[val$cell]
    return(out)
  })
  stopCluster(cl)
}
getDailyData <- function(chunk, runsDT, specDay, dirPath, dayF, nData){
  startRow <- as.numeric(unlist(runsDT[chunk, 1]))
  
  rowData <- specDay[grepl(paste0("start", startRow, "end"), specDay)]
  
  ## Read in data for particular chunk
  ffP <- file(file.path(dirPath, dayF, rowData), 'rb')
  dat <- readBin(ffP, 'integer', n=nData) # where "n=XX" is the number of instances you want to see
  close(ffP)
  
  return(dat)
}
createProbMap <- function(i, dayFiles, dirPath, runsDT, extVec, crsStr, 
                          nRow, nCol, resN, shapePath, nCuts, pal, plotDates,
                          fullOutDir, appLoc, pixIgnore, binaryDist){
  specDay <- list.files(file.path(dirPath, dayFiles[i]))
  
  ## Read in the data
  allData <- sapply(1:nrow(runsDT), getDailyData, runsDT, specDay, 
                    dirPath, dayF=dayFiles[i], nData=nCol*nCuts)
  
  ## Combine and rescale
  vecData <- unlist(allData)
  
  vecData <- round(vecData/10000, 2)
  
  oneDay <- terra::rast(nrows=nRow, ncols=nCol, 
                        xmin=extVec[[1]], xmax=extVec[[2]], 
                        ymin=extVec[[3]], ymax=extVec[[4]],
                        crs=crsStr)
  
  if(binaryDist | pal == "full"){
    vecData[pixIgnore] <- NA
    values(oneDay) <- vecData
    
    if(binaryDist){
      col <- c("#333333", "#0099CC")
      breaks <- NULL
      type <- "classes"
    } else {
      col <- viridis::plasma(11)
      breaks <- rev(seq(0,1,by=0.1))
      type <- "interval"
    }
  } else if(pal=="half"){
    col <- c("#333333", viridis::plasma(5))
    breaks <- rev(c(0, seq(0.5,1,by=0.1)))
    vecData <- ifelse(vecData < 0.5, 0, vecData)
    vecData[pixIgnore] <- NA
    values(oneDay) <- vecData
    type <- "interval"
    
    ## for edfm graphical abstract
    # makePal <- colorRampPalette(c("coral4", "orange", "yellow"))
    # col <- makePal(5)
    # base <- fread("dissertation/data/myanmar/trainingDataPoints.csv")
    # pt <- vect(base[pointid=="g1p0_7"], geom=c("coordX", "coordY"), crs="EPSG: 4326")
    # pt <- project(pt, "EPSG: 32646")
    # pt <- buffer(pt, 350)
    # test <- crop(oneDay, ext(pt))
    
    ## for edfm and agu poster
    # makePal <- colorRampPalette(c("coral4", "orange", "yellow"))
    # col <- c("#333333", makePal(5))
    # breaks <- rev(c(0, seq(0.5,1,by=0.1)))
    # vecData <- ifelse(vecData < 0.5, 0, vecData)
    # vecData[pixIgnore] <- NA
    # values(oneDay) <- vecData
  } else {
    stop("Must specify palette size for legend, either 'half' or 'full'")
  }

  if(grepl("Chatthin", appLoc)){
    pax <- list(yat=seq(2605000, 2620000, by=5000), 
                ylabs=seq(2605000, 2620000, by=5000))
  } else {
    pax <- list()
  }
  
  if(binaryDist){
    colNA <- "grey"
    levels <- c("No", "Yes")
    plg <- list(title="Deforested", cex=0.8)
  } else {
    colNA <- "black"
    levels <- c(0, 1)
    plg <- list(title="Probability", cex=0.8)
  }
  
  png(paste0(fullOutDir, "/", dayFiles[i], ".png"), units="cm", res=350,
      width=14, height=9.5)
  plot(oneDay, col=col, type=type, levels=levels, all_levels=TRUE,
       main=as.character(as.Date(plotDates[i], origin=as.Date("1970-01-01"))), 
       cex.main=0.8, colNA=colNA, pax=pax, breaks=breaks, plg=plg)
  sbar(d=2500, xy="bottomleft", label="2.5 km", col="white", cex=0.8, 
       halo=FALSE)
  
  ## for edfm and agu poster
  # par(bg="black") # only use for when looking at plot in rstudio
  # plot(oneDay, col=col, type="interval",
  #      main=paste0("Landscape probability of disturbance for ",
  #                  as.character(as.Date(plotDates[i],
  #                                       origin=as.Date("1970-01-01")))),
  #      colNA=NA,
  #      breaks=breaks, plg=list(text.col="white"),
  #      pax=list(col.axis="white", col.ticks="white", fg="white"),
  #      col.main="white", fg="white", bty="n")
  # box(col="white")
  
  if(grepl("Chatthin", appLoc)){
    chat <- vect(shapePath)
    bound <- terra::project(chat, crsStr)
    lines(as.lines(bound), col="white")
  }
  dev.off()
  
  return(print(paste0("Finished plotting Day ", i, "!")))
}
makeDailyImgs <- function(shapePath, plotCRS, extent, nRow, nCol, dirPath, 
                          timeFrame, matRuns, saveName, pal, dates, fullOutDir,
                          region, appLoc, pixIgnore, binaryDist){
  nCuts <- matRuns$nR[1]
  
  ## define other vars
  dateBreaks <- timeFrame
  runsDT <- as.data.table(matRuns)
  extVec <- extent
  resN <- 10
  crsStr <- plotCRS
  
  plotDates <- seq(dates[1], dates[length(dates)], by=dateBreaks)
  
  dayFiles <- list.files(dirPath)
  dayFiles <- dayFiles[order(as.numeric(gsub("day", "", dayFiles)))]
  dayFiles <- dayFiles[seq(1, length(dayFiles), by=dateBreaks)]
  
  # Create png maps for each day in the (sequenced) timeseries
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterEvalQ(cl, library(terra))
  clusterExport(cl, c("createProbMap", "getDailyData", "dayFiles", "dirPath",
                      "runsDT", "extVec", "crsStr", "nRow", "nCol",
                      "resN", "shapePath", "nCuts", "pal", "plotDates", 
                      "fullOutDir", "appLoc", "pixIgnore", "binaryDist"),
                envir=environment())
  parSapply(1:length(dayFiles), createProbMap, dayFiles, dirPath, runsDT,
         extVec, crsStr, nRow, nCol, resN, shapePath, nCuts, pal, plotDates,
         fullOutDir, appLoc, pixIgnore, binaryDist, cl=cl)
  
  stopCluster(cl)
  
  # createProbMap(i=100, dayFiles, dirPath, runsDT, extVec, crsStr, nRow, nCol, 
  #               resN, bound, nCuts, pal, plotDates, fullOutDir)
  
  # sapply(1:length(dayFiles), createProbMap, dayFiles, dirPath, runsDT,
  #        extVec, crsStr, nRow, nCol, resN, bound, nCuts, pal, plotDates,
  #        imgFolder)
  
  return(print("Done!"))
}

################################################################################
# Original functions
createRast <- function(K, mat=NULL, endR=NA, startR=NA, nR=NA, nC=NA){
  extFull <- c(744424.4, 772144.4, 2599900.0, 2621250.0) 
  
  # NB we have to subtract 1 from ymax because we are counting from topleft corner
  extIter <- c(extFull[1], extFull[2], extFull[4]-endR*10, extFull[4]-(startR-1)*10)
  
  r <- rast(nrows=nR, ncols=nC, 
            xmin=extIter[[1]], xmax=extIter[[2]], 
            ymin=extIter[[3]], ymax=extIter[[4]],
            crs="+proj=utm +zone=46 +datum=WGS84 +units=m +no_defs", 
            resolution=10)
  
  values(r) <- mat[K, ]
  return(r)
}
writeStatsRast <- function(dat=NULL, endR=NA, startR=NA, nR=NA, nC=NA, saveFile=""){
  statsMat <- matrix(unlist(dat), nrow=length(dat)*3)
  
  b <- sapply(1:length(dat)*3, createRast, mat=statsMat, endR=endR, startR=startR, nR=nR, nC=nC)
  
  outRast <- rast(b)
  names(outRast) <- unique(names(unlist(dat)))
  
  writeRaster(outRast, saveFile, overwrite=TRUE)
}
