##########################################################
## Purpose: 1) Extract downloaded L8 images, and
##          2) Prepare text files for creating the VRTs
##          3) Actually create the vrts
## Run medium:
##  - PC because we need the connection to SEAL with the windows path, but
##    technically could do Mac though would be slower
##  - For steps 2-3, make sure wd() is SEAL/IanMcGregor/
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, September 2021
## Last updated: July 2022
##########################################################
library(parallel)
library(stringr)
source("scripts/funs/appA_createVRTs.R")

# Step 1. Crop tifs and create the VRTs (AFTER making txt files)
library(terra)
library(data.table)
library(gdalUtilities)
library(parallel)

# workflow for each sensor:
# 1: get list of files
# 2: separate out bands to be the numeric and QAQC separately
# 3: build 2 vrts (1 for numeric and 1 for QAQC)
# 4: do gdalwarp on the vrts (cubic for numeric, and nearest neighbor for QAQC)

setwd("Z:/IanMcGregor/")

appRegion <- 5
hpc <- FALSE
applicationStep <- 3
source("scripts/argsApp.R", local=TRUE)
source("scripts/funs/appA_createVRTs.R")
origVars <- ls()
s2Bands <- c("NDVI", "SCL")
epsg="32646"

# crop tifs and make VRTs in parallel
## Note that cropTIF is only for S2, whereas makeVRT is for both L8 and S2
## NOTE the first time I ran the updated version with both L8 and S2, ALL of the
## final S2 vrts had binary 0 and 1s, and I couldn't figure out why. Running the
## S2 parallelization alone seemed to fix things.

sapply(1:length(appLocations), prepLandscapeApp, appLocations, cropTIF=FALSE, 
       makeVRT=TRUE, parL=TRUE, epsg, s2Bands)

###############################################################################
# Archive
# Extract downloaded L8 images from USGS

# base <- "/share/jmgray2/yourID/s2Data" # if data are in /share directory 
# base <- "/rsstu/users/j/jmgray2/SEAL/IanMcGregor/" # if in SEAL dir #for hpc
# base <- "/Volumes/SEAL/IanMcGregor/" # if in SEAL dir #for hpc
base <- "//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor/" #for windows

# Extract L8 images from downloaded tar files
l8 <- list.files(paste0(base, 
                        "l8c2Data/Bulk_Order_L8_Test1/Landsat8_OLI_TIRS_C2_L2"))

l8Folders <- gsub(".tar", "", l8)

cl <- makeCluster(7)
clusterExport(cl, c("base", "l8", "l8Folders"))
pbapply::pbsapply(1:length(l8), function(X){
  newF <- paste0(base, "l8c2Data/", l8Folders[X])
  if(!dir.exists(newF)) dir.create(newF)
  
  untar(paste0(base,"l8c2Data/Bulk_Order_L8_Test1/Landsat8_OLI_TIRS_C2_L2/", 
               l8[X]), exdir=newF, tar="internal")
}, cl=cl)

