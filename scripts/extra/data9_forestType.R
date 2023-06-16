##########################################################
## Purpose: Determine forest strata across training points in Myanmar
## Run medium:
##  - Mac, but might need PC
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, August 2022
## Last modified: August 2022
##########################################################
library(terra)
library(data.table)
source("scripts/funs/dataB_maskForest.R")
rootPath <- "data/dataMyanmar/"
maskPath <- paste0(rootPath, "spatial/forestMapMasked.tif")

################################################################################
# Step 1: Bring in forest strata map and mask to the strata we want
## note that this generally fails on Mac due to "insufficient disk space", but
### it works fine on PC.
## NOTE strata map comes from: https://smithsonian.figshare.com/articles/figure/Myanmar_Forest_Type_Map_2020/16613818/1?file=30747181
## NOTE we are outputting with EPSG: 32646 (46N) to match Chatthin
## NOTE reprojection will take a bit
rastPath <- paste0(rootPath, "spatial/Myanmar_Forest_Type_2020_utm47N_v1.tif")
outPath <- maskPath
maskNonForest(rastPath, outPath, plotRast=TRUE)

################################################################################
# Step 3: Get total number of pixels per strata in cropped, masked tif
land <- rast(maskPath)
# f <- as.data.table(freq(land))
# 
# strata <- data.table(value=0:8, landCover=c("mangrove", "lowland ever", 
#                                             "upland ever", "plantation",
#                                             "mix dec", "bamboo", "dry dec", 
#                                             "thorn", "swamp"))
# 
# totalPix <- sum(f$count)
# f[, `:=` (perc = round(count/totalPix, 2), landCover = strata$landCover)]

# Step 4: examine distribution of strata by training points
base <- fread(paste0(rootPath, "trainingPars/trainingDataPoints.csv"))
# base <- base[!grepl("a", pointid)]
train <- vect(base, geom=c("coordX", "coordY"), crs="EPSG: 4326")
train <- project(train, "EPSG: 32646")

forestPt <- as.data.table(extract(land, train, method="simple"))
colnames(forestPt) <- c("pointid", "class")
forestPt[, `:=` (pointid = base$pointid)]

## plot points over forest classes and assign nearest 
trainSub <- train[train$pointid %in% forestPt[is.na(class), pointid]]

buff <- buffer(trainSub, width=300)
n <- as.data.table(extract(land, buff))
colnames(n) <- c("ID", "class")

pts <- data.table(ID=1:nrow(trainSub), pointid=trainSub$pointid)
setkeyv(pts, "ID")
setkeyv(n, "ID")
n <- n[pts]

bl <- n[, .N, by=.(pointid, class)][!is.na(class), .SD[which.max(N)], by=pointid]

forestPt <- forestPt[!is.na(class)]
forestPt <- rbind(forestPt, bl[,N:=NULL])
forestPt <- forestPt[order(pointid), ]

b <- as.data.table(extract(land, trainSub, method="bilinear"))

totals <- forestPt[, .N, by=forestType
                   ][order(forestType),
                     ][, `:=` (perc = round(N/sum(N), 3))]
setkeyv(totals, "forestType")
setkey(strata, "value")
totals <- totals[strata]

totals[, percOriginal := f$perc]

###################################################################
library(terra)
shp <- vect("data/dataMyanmar/spatial/forestType90m/MMR_EcosystemsMap_v1_1.shp")
plot(shp)
