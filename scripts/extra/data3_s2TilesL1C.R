##########################################################
## Purpose: Get needed s2 tiles based on known points and automatically download
##          imagery from GEE. Then, download the s2 tiles from Google Cloud 
##          Storage
## Run medium: 
##  - Either Mac or PC for identifying the s2 tiles
##  - PC for downloading s2 L1C images bc already connected to SEAL
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 3.6.3, March-April 2021
## Last updated: July 2022
##########################################################
library(data.table)
library(terra)

s2TilesShp <- "data/spatial/s2Tiles/sentinel_2_index_shapefile.shp"
dataPath <- "data/dataMyanmar/trainingPars/trainingDataPoints.csv"

tileByPoints <- function(s2TilesShp, dataPath, chatthinOnly=FALSE, 
                         chatPath=NULL){
  # Step 1: Bring in data and make same CRS
  s2 <- vect(s2TilesShp)
  base <- fread(dataPath)
  
  v <- vect(base, geom=c("coordX", "coordY"), crs="EPSG:4326")
  if(crs(v) != crs(s2)) s2 <- project(s2, v)
  
  ##look at specific sub-section if want
  # subP <- base[coordX > 98 & coordX < 100 & coordY > 22 & coordY < 23.5, ]
  # subP <- vect(subP, geom=c("coordX", "coordY"), crs="EPSG:4326")
  
  if(chatthinOnly){
    # already EPSG: 4326
    chat <- vect(paste0("data/dataMyanmar/", chatPath))
    s2Small <- crop(s2, ext(chat))
  } else {
    s2Small <- crop(s2, ext(v))
  }
  
  ## crop s2 tile extent to only focus on smaller data
  plot(s2Small)
  text(s2Small, labels=s2Small$Name)
  plot(v, add=TRUE)
  # text(dt, labels=dt$pointid)
  
  ## relabel v if using chatthin
  if(chatthinOnly) plot(chat, add=TRUE) 
  
  # Step 2: Get unique tile names that contain points
  ## we need to use "intersects" here otherwise if we use something like
  ## "contains" we actually exclude different points. The wiki help page listed
  ## in the help docs for this function is helpful
  vInt <- as.data.table(t(relate(s2Small, v, relation="intersects")))
  colnames(vInt) <- s2Small$Name
  pointids <- v$pointid
  
  ptTile <- rbindlist(lapply(colnames(vInt), function(X){
    sub <- vInt[,..X]
    pts <- which(sub[[1]])
    if(length(pts)==0){
      pointid <- NA
    } else {
      pointid <- pointids[pts]
    }
    out <- data.table(tile = X, pointid = pointid)
    return(out)
  }))
  
  # Remove NAs, then add in lat/lon info
  ptTile <- ptTile[!is.na(pointid)]
  ptTile[, `:=` (lon = base$coordX[match(pointid, base$pointid)],
                 lat = base$coordY[match(pointid, base$pointid)])]
  
  return(list(s2Full=s2, s2Crop=s2Small, pointsDT=dt, pointsShp=v, 
              pointsIntersect=vInt, ptTile=ptTile))
}

# 1. Get the unique tile by points
tileData <- tileByPoints(s2TilesShp, dataPath)
tileData <- tileByPoints(s2TilesShp, dataPath, chatthinOnly=TRUE,
                         chatPath="spatial/CWS/CWS_Boundary.shp")

# 2. Filter the tiles
## This is the subjective part. Many tiles are overlapping with others, and
## often, points have two separate readings. This is fine and accounted for in
## the main algorithm, but processing S2 tiles ourselves is memory-intensive.
## To filter out unnecessary tiles, we only retain tiles that have >=4 points
## AND that are not clearly overlapping another tile upon inspection.

## Let's first look at the quantitative distribution.
tt <- tileData$ptTile[,.N,by=.(tile)][order(-N),]

## clearly, we have a number of tiles that just contain 1 point, but looking at 
### the plot, we can see that each of those is overlapped by another tile that 
### has more points.

## can make a bar plot if want
library(ggplot2)
ggplot(as.data.frame(tt)) +
  aes(x = reorder(tile, -N), y=N) +
  geom_bar(stat="identity") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=0.5))

## see which points are doubled with tiles
ptTile <- tileData$ptTile
tilesLow <- tt[N < 8, tile]
test <- rbindlist(lapply(tilesLow, function(X){
  focalPts <- ptTile[tile==X, pointid]
  allTiles <- unique(ptTile[pointid %in% focalPts, tile])
  allTiles <- allTiles[!(allTiles %in% tilesLow)]
  
  if(length(allTiles) > 0){
    out <- data.table(originalTile=X, coverTile=allTiles)
  } else {
    out <- data.table(originalTile=X, coverTiles=NA)
  }
  return(out)
}))

# 3. Only keep the main tiles and write to csv
## MAKE SURE TO COPY TO "SEAL/Ian/diss/data/myanmar/"
outTable <- tileData$ptTile
keep <- tt[N > 4, ]
outTable <- outTable[tile %in% keep$tile]
fwrite(outTable, "data/dataMyanmar/trainingPars/trainingDataS2Tiles.csv")

### -------- Download S2 images from Google Cloud Storage ------------------####
# Download s2 images from GEE (same as if downloading from Copernicus)
## **Faster to do this on PC**
### If this is the first time you're doing this, you need to download and
### install gsutil first (https://cloud.google.com/storage/docs/gsutil_install)


library(data.table)
# devtools::install_github("ncsuSEAL/sealR")
library(sealR)
dat <- fread("data/dataMyanmar/trainingDataS2Tiles.csv")
filePath <- "//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor/s2Data/L1C"

currentTiles <- list.files(filePath)
tiles <- setdiff(unique(dat$tile), currentTiles)

## this downloads for all images between May 2015 and Feb 2020. If it needs to
## be different, change the regExp
for(i in 1:length(tiles)){
  DownloadSentinel2(tiles[i], filePath,
                    L2Flag = FALSE,
                    regExp = "S2\\D_MSIL\\d\\D_20([1][5-9]|[2][0][0][1]).*.SAFE/$")
}
