##########################################################
## Purpose: Calculate seasonality of residuals to assign a varying inflation 
##          (for training data prior to logistic model fitting; resDates=TRUE),
##          AND compare the variance in the monitoring period residuals with the
##          stable period (for landscape application static inflation; 
##          resDates=FALSE)
## Run medium:
##  - Either
## Figures for McGregor et al 2023: 
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, October 2022
## Last modified: October 2022
##########################################################
groundhog.library(data.table, groundhogDate)
source("scripts/funs/trainB_inflation.R")

## main functions
# in `inflationFactor()`, the arguments mean
## ndvi = TRUE only used for plotting figures
## resDates = TRUE & variance = FALSE to calculate training dynamic inflation
## resDates = FALSE & variance = TRUE to calculate landscape static inflation

runInflation <- function(nRun, plotSeasons, plotInflation, calcFactor, 
                         timeframe){
  variance <- FALSE
  resDates <- FALSE
  ndvi <- FALSE
  windowVec <- NULL
  returnData <- FALSE
  
  if(plotSeasons) ndvi <- TRUE
  if(plotInflation) resDates <- TRUE
  if(calcFactor) variance <- TRUE; returnData <- TRUE; windowVec <- c(45)
  
  combos <- c("sentinel1", "landsat8, sentinel2", 
              "landsat8, sentinel2, sentinel1")
  fileNames <- c("S1", "L8S2", "All")
  
  numRun <- as.numeric(nRun)
  sensKeep <- combos[numRun]
  script <- 1
  bayes <- FALSE # necessary for args.R
  source("scripts/args.R", local=TRUE)
  source("scripts/argsTrain.R", local=TRUE)
  
  # bring in ts residual data
  fileSaveLoc <- paste0(dataPath, "/trainingPars/train1_", fileNames[numRun])
  load(paste0(fileSaveLoc, ".Rdata"))
  
  ## Now run the internal functions
  resVec <- inflationFactor(oneTS, variance, resDates, ndvi, sensVec=nRun, 
                            windowVec, timeframe, returnData, dates=dates)
  
  if(plotInflation){
    # items are 
    ## 1 = original residuals
    ## 2 = doy
    ## 3 = inflation factor
    return(resVec)
  } else {
    return(print("Done"))
  }
}
makePlots <- function(plotSeasons=FALSE, plotInflation=FALSE, calcFactor=FALSE){
  if(plotSeasons) fileAdd <- "seasonality"
  if(plotInflation) fileAdd <- "residualVec"
  
  tiff(paste0("figures/figXX_", fileAdd, ".tiff"), 
      width=22, height=14, units="cm", res=350, compression="lzw")
  
  if(plotSeasons) layout(matrix(1:6, nrow=2, byrow=TRUE))
  if(plotInflation) layout(matrix(1:4, nrow=2, byrow=TRUE))
  sapply(c("stable", "monitor"), function(tf){
    sapply(c(2,3), runInflation, plotSeasons, plotInflation, calcFactor=FALSE,
           timeframe=tf)
  })
  
  dev.off()
}

##----------------------------------------------------------------------------##
# 1. Get seasonal variation in residuals for dynamic inflation factor by doy
## NOTE 1: this uses output from train1 with inflationType="static"
## NOTE 2: this function uses the original residuals and ts models for NDVI plot
## visualization. Otherwise, the seasonal inflation adjustment (by doy) AND the
## monitoring inflation factor are based on the aggregated residuals 
## (standardized and twice-spike filtered)

## plot seasonality or inflation vector
# makePlots(plotSeasons=TRUE)
# makePlots(plotInflation=TRUE)

## calculate inflation factor for dynamic inflation (ignore plots here)
### as a reminder, this is only calculated from stable period for this iteration
### correct params are calcFactor=FALSE and timeframe="stable"
if(inflationType=="static"){
  out <- lapply(c(2,3), runInflation, plotSeasons=FALSE, plotInflation=TRUE,
              calcFactor=FALSE, timeframe="stable")
  out <- data.table("L8S2" = out[[1]][[3]], "All" = out[[2]][[3]])

  ## save the stable period dynamic inflation factors
  sapply(colnames(out), function(X){
    vec <- out[, get(X)]
    save(vec, file=
          paste0("data/trainingPars/train2_seasonalAdjustment", X, ".Rdata"))
    return(print("Saved"))
  })
}

##----------------------------------------------------------------------------##
# 2. Get static inflation factor based on increased variance in 
## monitoring residuals FOR LANDSCAPE APPLICATION ONLY. 
## The calculation of this step is dependent on the seasonally-adjusted 
## residuals, so make sure script 1 has been run with the "dynamic" 
## inflationType set.

## STOP. HAVE YOU RUN SCRIPT 1 WITH dynamic INFLATION YET?
## correct parameters are nRun=2, calcFactor=TRUE, other arguments FALSE, 
## timeframe="monitor"
if(inflationType=="dynamic"){
  runInflation(nRun=2, plotSeasons=FALSE, plotInflation=FALSE, calcFactor=TRUE,
             timeframe="monitor")
}

##----------------------------------------------------------------------------##
## archive: try different windows
# layout(matrix(1:4, ncol=2, byrow=TRUE))
# inflationFactor(nodist, variance=FALSE, resDates=TRUE, ndvi=FALSE,
#                 sensVec=c(1:3), windowVec=c(15,30,45,60), timeframe="monitor",
#                 returnData=TRUE, dates=dates)

# layout(matrix(1:4, ncol=2))
# inflationFactor(nodist, variance=FALSE, resDates=FALSE, ndvi=TRUE,
#                 sensVec=c(1:3), timeframe="monitor", returnData=TRUE, dates=dates)
