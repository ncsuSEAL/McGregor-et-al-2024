##########################################################
## Purpose: Spot-check different pixels in landscape maps and
##          look at individual ts plots
## Run medium:
##  - PC because already connected to SEAL
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: July 2022
##########################################################
# setwd("Z:/IanMcGregor/")
library(data.table)
library(terra)
library(MuMIn)
source("dissertation/scripts/funs/appC_analyzeLandscape.R")
source("dissertation/scripts/funs/commonB_modZ.R")
source("dissertation/scripts/funs/commonC_dupNA.R")
source("dissertation/scripts/funs/commonD_ewmaBinaryCat.R")
source("dissertation/scripts/funs/commonE_plotSummaries.R")
source("dissertation/scripts/funs/commonF_probBayes.R")

# coords must be degrees in lon-lat order or NULL
coords <- c(95.581458, 23.589703) #appRegion5
coords <- c(97.167665, 24.861578) #appRegion1
coords <- c(97.53794, 25.23482) #g0p1_7
spotCheck(coords, coordFromProbs=FALSE, appRegion=5)
spotCheck()

# different spotchecks showing the issues with each landscape
appRegion <- 1

if(appRegion==1){ #lowEver
  coords <- c(920934.8, 2756053) #day1499, i=108
  coords <- c(942835.5, 2761121) #day1653, i=119
} else if(appRegion==2){ #upEver
  coords <- c(766519.1, 2844475) #day1611, i=116, road
  coords <- c(771564.2, 2843572) #day1541, i=111
  coords <- c(783460, 2845218) #day1625, i=117
} else if(appRegion==3){ #mixDry1
  coords <- c(625763.4, 2577231) #day337, i=25
  coords <- c(619072, 2581108)
} else if(appRegion==4){ #mixDry2
  coords <- c(844815.3, 2569880) #day1331, i=96
  coords <- c(842584.9, 2566215) #day1583, i=114
} else if(appRegion==5){ #chatthin
  coords <- c(768339.3, 2606137) #day1359, i=98
  coords <- c(767702.2, 2617338) #day1667, i=120, not forest, outside chatthin
  coords <- c(754483.4, 2609110) #same as above, inside chatthin
}

coords <- c(924003.9, 2763336)
out <- spotCheck(coords, coordFromProbs=TRUE, appRegion, plots=TRUE, bayes=TRUE)

################################################################################
# plot the daily dist ratios across landscape types
library(data.table)
library(terra)

applicationStep <- 3
script <- 1
setwd("Z:/IanMcGregor/")
bayes <- TRUE
source("dissertation/scripts/args.R", local=TRUE)
source("dissertation/scripts/funs/appD_probMaps.R")

layout(matrix(1:2, ncol=2, byrow=TRUE))
sapply(c(1,5), function(X){
  appRegion <- X
  source("dissertation/scripts/argsApp.R", local=TRUE)
  dates <- getDates(saveDir, fold=paste0(appLoc, "/landscape/"))
  fileLoad <- paste0("dissertation/data/myanmar/", appLoc, 
                     "/dailyDistRatio.Rdata")
  if(bayes){
    fileLoad <- gsub(".Rdata", "_bayes.Rdata", fileLoad)
  }
  
  load(fileLoad)
  
  ## convert only the initial 0s back to NAs (representing no observations)
  nas <- c(0, diff(which(out==0)))
  if(any(which(nas > 1))){
    out[1:(first(which(nas > 1))-1)] <- NA
  } else {
    out[1:length(nas)] <- NA
  }
  
  plot(y=out*100, x=as.Date(dates, origin=as.Date("1970-01-01")), 
       xlim=c(as.numeric(as.Date("2015-06-23")), 
              as.numeric(as.Date("2020-01-31"))),
       xlab="Date", ylab="% landscape disturbed", 
       main=paste0("Region ", appRegion), type="l")
  abline(v=as.numeric(as.Date("2018-12-31")), lty=2)
})

# do the GLAD comparison
cellComp <- lapply(1:5, function(X){
  appRegion <- X
  source("dissertation/scripts/argsApp.R", local=TRUE)
  
  # Get masked NA values
  maskRegion <- createForestMask(saveDir, appRegion, nMatRow, nCol, 
                                 appLocations, perCover=30)
  pixIgnore <- which(is.na(values(maskRegion)))
  
  # Get dist 19 pixels
  dist19 <- createForestMask(saveDir, appRegion, nMatRow, nCol, 
                                 appLocations, perCover=30, returnGlad=TRUE)
  valsGlad <- values(dist19)
  distGlad <- which(!is.na(valsGlad))

  # constrain Glad values with ignore pixels
  distGlad <- distGlad[!(distGlad %in% pixIgnore)]
  
  # bring in my detections
  load(paste0("dissertation/data/myanmar/", appLoc, "/distCellIndex.Rdata"))
  
  # cells together
  both <- out[out %in% distGlad]
  onlyGlad <- setdiff(distGlad, both)
  onlyMine <- setdiff(out, both)
  
  # create empty plot
  ex <- appLocations[[appRegion]]$ext
  oneDay <- terra::rast(nrows=nMatRow, ncols=nCol, 
                        xmin=ex[[1]], xmax=ex[[2]], 
                        ymin=ex[[3]], ymax=ex[[4]],
                        crs=appLocations[[appRegion]]$crs)
  
  # set up values vector
  ## 0 = nodist, 1=my dist, 2=glad dist, 3=both
  vecFull <- rep(0, length(valsGlad))
  vecFull[both] <- 3
  vecFull[onlyGlad] <- 2
  vecFull[onlyMine] <- 1
  vecFull[pixIgnore] <- NA
  
  # plot the map
  values(oneDay) <- vecFull
  legendList <- list(legend=c("0=No detection", "1=My detection", "2=GLAD", "3=Both"),
                     col=c("grey", "blue", "red", "yellow"))
  plot(oneDay, main=paste0("2019 Disturbance Detections for Region ", appRegion),
       col=c("grey", "blue", "red", "yellow"), colNA="black", plg=legendList)
  
  print(paste0("Done with region ", appRegion))
  return(list(onlyMe = length(onlyMine), onlyGlad = length(onlyGlad), 
              both = length(both)))
})

################################################################################
# Compare raster values between training and landscape
library(data.table)
library(terra)

source("scripts/funs/commonE_plotSummaries.R")
source("scripts/funs/trainA_combineDat.R") # step 1
source("scripts/funs/commonA_processDat.R") # step 1

base <- fread("data/myanmar/trainingDataPoints.csv")
chatPts <- base[grepl("g42", pointid), pointid]
chatPts <- chatPts[!grepl("a", chatPts)]

compareTrainLand(chp="g42p2_5", base)

################################################################################
# get band values of 9-pixel neighborhood

library(terra)

# surrounding cells
nCol <- 2773

getCellNum <- function(type, nCol){
  cell0 <- ifelse(type=="chatthin", 3776126, 177086)
  cellValsnRun <- unlist(list(cell1 = cell0 - nCol - 1,
                              cell2 = cell0 - nCol,
                              cell3 = cell0 - nCol + 1,
                              cell4 = cell0 - 1,
                              cell5 = cell0 + 1,
                              cell6 = cell0 + nCol -1,
                              cell7 = cell0 + nCol,
                              cell8 = cell0 + nCol + 1))
  return(as.numeric(sort(c(cell0, cellVals))))
}
cells <- getCellNum(type="chatthin", nCol)

getCellValues <- function(path, img, bands, cells){
  f <- list.files(paste0(path, img), full.names=TRUE)
  f <- f[grepl(bands, f)]
  
  rst <- rast(f)
  
  b <- sapply(cellVals, function(X){
    coords <- as.vector(xyFromCell(r, X))
    coordVec <- vect(data.frame(lon=coords[1], lat=coords[2]), 
                     crs="EPSG: 32646")
    plot(coordVec, add=TRUE)
    return(extract(rst, coordUTM))
  })
  
  return(b)
}

path <- "Z:/IanMcGregor/"
img <- "l8c2Data/LC08_L2SP_134044_20190118_20200830_02_T1"
bands <- "B4|B5|QA_PIXEL"

img <- "s2Data/L2A/46QGM/S2A_MSIL2A_20190116T041121_N0300_R047_T46QGM_20210916T004952.SAFE/GRANULE/L2A_T46QGM_A018634_20190116T041933/IMG_DATA/R20m"
bands <- "B04|B8A|SCL"

getCellValues(path, img, bands, cells)





################################################################################
# Remove older probability files automatically
## I'm doing this the easy way by just looking at the file names, but we could
## do this comparing the modified date and memory as well
filePath <- "dissertation/data/myanmar/landscape/probs"
f <- list.files(filePath)
filesKeep <- paste0("start", matRuns$startR, "end", matRuns$endR)

library(parallel)
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("filePath", "f", "filesKeep"))
parSapply(cl, f, function(X, keep=filesKeep){
  probFiles <- list.files(file.path(filePath, X))
  filesRemove <- probFiles[!grepl(paste0(keep, collapse="|"), probFiles)]
  
  sapply(filesRemove, function(Y, day=X){
    file.remove(file.path(filePath, day, Y))
    return(print(paste0("File removed: ", day, "/", Y)))
  })
})

# for 60 rows at a time, this should be 36
b <- sapply(f, function(X) return(length(list.files(file.path(filePath, X)))))

stopCluster(cl)
