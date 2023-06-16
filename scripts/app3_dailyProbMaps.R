##########################################################
## Purpose: Create PNGs of daily lanscape probability maps,
##          convert to a GIF, and/or calculate the ratio of 
##          disturbed pixels (prob > 0.5) for every day of the ts
## Input: all necessary functions
## Variable created: plots
## Run medium:
##  - HPC job for getDistRatio=TRUE. Csh file with same name is located
##    at SEAL/Ian/
##  - HPC interactive session for makeDailyMaps=TRUE. This is much faster than
##     Mac or PC. Inexplicably, actually submitting a job fails with no output 
##     files.
##  - Mac or PC only for makeGIF=TRUE. 
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.2, Dec 2021
## Last updated: June 2022
##########################################################
library(data.table)
library(terra)
library(parallel)

# setwd("Z:/IanMcGregor/")

runProbMaps <- function(appRegion, getDistRatio, compareGLAD, makePNGs, 
                        makeGIF, bayes, binaryDist){
  # Define key variables and functions
  applicationStep <- 3
  hpc <- TRUE
  source("scripts/argsApp.R", local=TRUE)
  source("scripts/args.R", local=TRUE)
  source("scripts/funs/appD_probMaps.R")
  
  # as a reminder, when you assign values to a raster, it starts at top left, 
  ## and goes left-right, top-bottom
  
  # Get exact dates used for probability calculation and define vars
  dates <- getDates(saveDir, fold=paste0(appLoc, "/landscape/"))
  dirPath <- paste0(saveDir, probsDir)
  if(bayes & !binaryDist) dirPath <- paste0(dirPath, "Bayes")
  if(binaryDist) dirPath <- gsub("probs", "binaryProbs", dirPath)
  
  # Here we are masking out any pixels that don't have >= perCover % forest as
  ## of the beginning of 2019
  maskRegion <- createForestMask(saveDir, appRegion, nMatRow, nCol, 
                                 appLocations, perCover=30)
  pixIgnore <- which(is.na(values(maskRegion)))
  
  if(!binaryDist){
    if(getDistRatio){
      write(paste0("Starting to calculate daily ratio for region ", appRegion, 
                   " at ", format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), 
            file=txtOutFile, append=TRUE)
      out <- dailyDistMetric(dates, dirPath, matRuns, nCol, nCuts, probThresh,
                             pixIgnore, ratio=TRUE, index=FALSE)
      if(bayes) ratioFile <- gsub(".Rdata", "_bayes.Rdata", ratioFile)
      save(out, file=ratioFile)
      print(paste0("Saved ratio of disturbance to ", ratioFile))
      
      write(paste0("Finished calculating daily ratio at ", 
                   format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
            append=TRUE)
    }
    
    if(compareGLAD){
      write(paste0("Starting to calculate GLAD comparison for region ", appRegion, 
                   " at ", format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), 
            file=txtOutFile, append=TRUE)
      out <- dailyDistMetric(dates, dirPath, matRuns, nCol, nCuts, probThresh,
                             pixIgnore, ratio=FALSE, index=TRUE)
      save(out, file=distCellFile)
      print(paste0("Saved disturbance index cells to ", distCellFile, " for 
                 GLAD comparison."))
      write(paste0("Finished calculating GLAD comparison at ", 
                   format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
            append=TRUE)
    }
  }
  
  if(makePNGs){
    write(paste0("Starting to make daily png maps at ", 
                 format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
          append=TRUE)
    
    # Create probability maps per day
    if(binaryDist) imgFolder <- paste0(imgFolder, "Binary")
    if(!binaryDist & bayes) imgFolder <- paste0(imgFolder, "Bayes")
    fullOutDir <- paste0(saveDir, appLoc, "/", imgFolder)
    if(!dir.exists(fullOutDir)) dir.create(fullOutDir)
    
    region <- appLocations[[appRegion]]
    extent <- region$ext
    
    makeDailyImgs(shapePath, plotCRS, extent, nRow=nMatRow, nCol, dirPath, 
                  timeFrame, matRuns, saveName, pal, dates, fullOutDir, region,
                  appLoc, pixIgnore, binaryDist)
    
    write(paste0("Finished creating daily png maps at ", 
                 format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
          append=TRUE)
  }
  
  ## GIF conversion can ONLY happen on Mac or PC
  if(makeGIF){
    write(paste0("Starting to convert pngs to gif at ", 
                 format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
          append=TRUE)
    
    # Convert series of images into gif
    pngFiles <- list.files(gsub("probs", "dailyMaps", dirPath), 
                           full.names=TRUE)
    pngFiles <- pngFiles[grepl(".png", pngFiles)]
    ord <- as.numeric(gsub(".*day", "", gsub(".png", "", pngFiles)))
    pngFiles <- pngFiles[order(ord)]
    
    ## smaller gif for defense
    # bookEnds <- which(grepl("1275|1429", pngFiles))
    # pngFiles <- pngFiles[bookEnds[1]:bookEnds[2]]
    # gifName <- paste0(gsub("probs", "defenseGIF", dirPath), ".gif")
    
    gifName <- paste0(gsub("probs", paste0("fullTS_", timeFrame), dirPath), 
                      ".gif")
    
    gifski::gifski(pngFiles, 
                   gif_file=gifName, 
                   width=900, height=600, delay=0.5)
    
    write(paste0("Finished creating daily gif at ", 
                 format(Sys.time(), "%Y-%m-%d_%H:%M:%S")), file=txtOutFile, 
          append=TRUE)
  }
}

runProbMaps(appRegion=1, getDistRatio=FALSE, compareGLAD=FALSE, makePNGs=TRUE, 
            makeGIF=FALSE, bayes=FALSE, binaryDist=FALSE)
