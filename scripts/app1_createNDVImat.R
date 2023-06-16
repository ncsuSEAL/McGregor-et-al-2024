##########################################################
## Purpose: Create landscape matrix for N rows using data from VRTs.
##          Specifically, read in N rows of data from VRTs, filter by QAQC,
##          then calculate index (here, NDVI) and return that matrix.
## Input: all necessary functions
## Variable created: matrices to be used in 9AppTask2-3.R
## Run medium:
##  - HPC submitted job via csh file of the same name in SEAL/Ian/
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last updated: June 2022
##########################################################
library(data.table)
library(parallel)
# TO RUN THIS SCRIPT - you must `cd` to the root folder, 
## e.g. cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
## this is because of how the VRT xml relates to the raster locations

# NB This code will sometimes fail with the error "Failure during Raster IO". 
## The reason for this is due to available memory, i.e. the amount of space 
## needed to have several VRTs open via GDAL at once. 
## E.g. this will work on Mac if nothing else is running, but will give error 
## otherwise (e.g. if zoom is in the background)
## This DOES run fine on HPC, but memory of the node / worker needs to be at 
## least ~16GB. Other potential workaround is simply processing fewer rows at 
## a time, but this would involve a full re-run of scripts along with redoing 
## the `matRuns` table.

# setwd("Z:/IanMcGregor/")
# source("dissertation/scripts/9App1_functions.R")

args <- commandArgs(trailingOnly=TRUE)

applyLandMatrix <- function(nRun, appRegion){

  applicationStep <- 1
  nRun <- as.numeric(nRun)
  
  # Step 1. Define and load key variables / define parameters
  scriptPath <- "scripts/"
  source(paste0(scriptPath, "args.R"), local=TRUE)
  source(paste0(scriptPath, "argsApp.R"), local=TRUE)
  source(paste0(scriptPath, "funs/appB_combineDat.R"))
  source(paste0(scriptPath, "funs/commonA_processDat.R"))
  
  # Step 2: Write out progress to text file
  write(paste0("Started task for ", appLoc, ":", task1, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), "."),
        file=txtOutFile, append=TRUE)
  
  # Step 3: Process the raw satellite data
  ## Read in data, convert to matrix, filter by QAQC and calculate index
  dataVecs <- createDataVectors(base=baseDir, saveDir, sensors, appLoc)
  allMatrices <- lapply(1:length(dataVecs$vrts), processMatrix, 
                        bandList=dataVecs$bandList,
                        vrts=dataVecs$vrts, startRow=startRow, 
                        numberRow=numberRow, nMatRow=nMatRow, nCol=nCol,
                        fileOut=txtOutFile, funs=funs, appLoc=appLoc,
                        vegIndex=vegIndex)
  names(allMatrices) <- names(dataVecs$vrts)
  
  if(saveSensorTS){
    saveRDS(allMatrices, 
            file=paste0(gsub("landscape/", 
                             paste0(appLoc, "/landscape/"), saveSensorTSDir), 
                        "mats", startRow, "-", endRow, ".RDS"))
  }
  
  # Step 4: Write out progress to text file
  write(paste0("Finished task for ", appLoc, ":", task1, " at ", 
               format(Sys.time(), "%Y-%m-%d_%H:%M:%S"), "."),
        file=txtOutFile, append=TRUE)
  return(print(paste0("Successfully done with rows ", startRow, "-", 
                      numberRow+startRow-1, 
                      " and saved the RDS file if indicated.")))
}
applyLandMatrix(nRun=args[1], appRegion=1)
# applyLandMatrix(nRun=18)
#############################################################################################
# if want to run the original script (takes longer per point, more memory), then convert the
## output of the first script to a data frame here and run the original scripts
## This can be done for spotchecking a pixel if need be

# outDT <- lapply(names(allMatrices), function(H){
#   matr <- allMatrices[[H]]$ts
#   colnames(matr) <- as.vector(sapply(vrts[[H]], function(D){
#     return(paste0("x", stringr::str_extract(D, "20[[:digit:]]{6}")))
#   }))
#   
#   rownames(matr) <- paste0("pixel_", seq(1:nrow(matr)))
#   tt <- as.data.table(as.table(matr))
#   colnames(tt) <- c("pixel", "date", "bandValue")
#   return(tt)
# })
# 
# outDT <- rbindlist(outDT)