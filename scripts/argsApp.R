##########################################################
## Purpose: Apply ts models and prob model at landscape scale to get probabilities of disturbance
## Input: all necessary functions
## Variable created: matrices of daily prob of disturbance per defined chunks, saved by day
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Dec 2021
##########################################################

# All application steps
runType <- "app"
sensors <- "L8S2" #L8S2 or All

## Note that for landscape application (and in general), it doesn't matter if
## backfill is labeled as TRUE or if we do FALSE and then add do the filling 
## just prior to plotting. In other words, the plots look the same. To make
## things easier with how the code is written, I'm keeping this as TRUE.
backfill <- FALSE

### how split up landscape---------------------------------------------------###
nMatRow <- 2135                         #max row of landscape matrix
nCol <- 2773                            #max col of landscape matrix
nBreaks <- 60                           #N rows to consider at a time
baseN <- floor(nMatRow/nBreaks)

## this was with DT but kept failing in HPC (but not outside HPC)
matRuns <- data.frame(startR=seq(1,nMatRow, by=nBreaks), 
                      nR=c(rep(nBreaks, baseN), nMatRow-baseN*nBreaks))
matRuns$endR <-  matRuns$startR + matRuns$nR-1

### folders and files -------------------------------------------------------###
# baseDir="/Volumes/SEAL/IanMcGregor/"

## define folder locations and load specific Rdata for application region
saveDir <- "dissertation/data/myanmar/"
saveSensorTSDir <- paste0(saveDir, "landscape/")

load(paste0(saveDir, "appLocations.Rdata"))
appLoc <- names(appLocations)[appRegion]
substr(appLoc, 1, 1) <- toupper(substr(appLoc, 1, 1))
appLoc <- paste0("app", appLoc)

probsDir <- paste0(appLoc, "/probs")

## NB appLocations
#1 = lowEver
#2 = upEver
#3 = mixDry1
#4 = mixDry2
#5 = chatthin

# Application step-dependent variables --------------------------------------###
## specific rows to consider for this iteration
if(applicationStep %in% c(1,2,4,5)){
  startRow <- as.numeric(matRuns[nRun, 1])
  numberRow <- as.numeric(matRuns[nRun, 2])
  endRow <- as.numeric(matRuns[nRun, 3])
  
  modelFitsOut <- paste0(saveDir, appLoc, "/landscape/start", startRow,
                         "end", endRow, "_")
  outFileProbs <- paste0(saveDir, appLoc, "/probs/allProbs_start",
                         startRow, "end", endRow)
}

## Read in data, format, get index
if(applicationStep==1){
  saveSensorTS <- TRUE
  
  task1 <- paste0(" reading in VRT ts data for ", numberRow, 
                  " rows starting at Row ", startRow)
  txtOutFile <- paste0("jobs/app1_createNDVImat/appLogTask1.txt")
  if(!file.exists(txtOutFile)) file.create(txtOutFile)
}

# Calculate daily probabilities for each pixel
if(applicationStep %in% c(2,5)){
  # nStableYears <- 3.5   #length of stable period in years (3.5 = ~end of 2018)
  endStable <- as.numeric(as.Date("2018-12-31"))
  combined <- FALSE
  longTime <- FALSE
  
  outText <- paste0(" processing for nRun number ", nRun)
  txtOutFile <- paste0("jobs/app2_processLandProbs/appLogTask2.txt")
  if(!file.exists(txtOutFile)) file.create(txtOutFile)
  
  ## load the glm model generated from the training data
  logMod <- readRDS(paste0(saveDir, "probModel.RDS"))
  
  ## load the RDS file of the rows to be processed for this iteration
  # allMatrices <- readRDS(paste0(baseDir, saveDir, "landscape/mats", startRow, 
  #                               "-", endRow, ".RDS"))
  allMatrices <- readRDS(paste0(saveDir, appLoc, "/landscape/mats", startRow, 
                                "-", endRow, ".RDS"))
}

# Daily probability maps
if(applicationStep %in% c(3,4)){
  timeFrame <- 14   # number of days between each map across the full timeseries
  pal <- "half"     # either full or half
  
  txtOutFile <- paste0("jobs/appLogTask3.txt")
  if(!file.exists(txtOutFile)) file.create(txtOutFile)
  
  plotCRS <- "EPSG: 32646"
  shapePath <- paste0(saveDir, "shapefiles/CWS_Boundary.shp")
  
  saveName <- paste0(saveDir, "landscape/chatProbs20220728_")
  imgFolder <- "dailyMaps"
  ratioFile <- paste0(saveDir, appLoc, "/dailyDistRatio.Rdata")
}
