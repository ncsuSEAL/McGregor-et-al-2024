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
library(data.table)
source("scripts/funs/trainB_inflation.R")

## main functions
runInflation <- function(nRun, plotSeasons, plotInflation, calcFactor){
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
  resVec <- lapply(c("stable", "monitor"), function(tf){
    inflationFactor(oneTS, variance, resDates, ndvi, sensVec=nRun, windowVec,
                    timeframe=tf, returnData, dates=dates)
  })
  
  if(plotInflation){
    ## item list: 1=original res, 2=doy, 3=std, or inflation factor
    names(resVec) <- c("stable", "monitor")
    inflation <- resVec$stable[[3]]
    return(inflation)
  } else {
    if(calcFactor){
      print("Please update `args.R` with the mean combined inflation factor")
    }
    return(print("Done"))
  }
}
makePlots <- function(plotSeasons=FALSE, plotInflation=FALSE, calcFactor=FALSE){
  if(plotSeasons) fileAdd <- "seasonality"
  if(plotInflation) fileAdd <- "residualVec"
  
  png(paste0("figures/figXX_", fileAdd, ".png"), 
      width=20, height=14, units="cm", res=350)
  layout(matrix(1:4, nrow=2))
  sapply(c(2,3), runInflation, plotSeasons, plotInflation, calcFactor=FALSE)
  dev.off()
}

# 0. Make plots showing seasonal variation for both NDVI and backscatter, and
## seasonal variation in the residuals
makePlots(plotSeasons=TRUE)
makePlots(plotInflation=TRUE)

##----------------------------------------------------------------------------##
# 1. Get seasonal variation in residuals for dynamic inflation factor by doy
## NOTE: this function uses the original residuals and ts models for NDVI plot
## visualization. Otherwise, the seasonal inflation adjustment (by doy) AND the
## monitoring inflation factor are based on the aggregated residuals 
## (standardized and twice-spike filtered)

## calculate inflation factor for dynamic inflation (ignore plots here)
### as a reminder, this is only calculated from stable period for this iteration
out <- sapply(c(2,3), runInflation, plotSeasons=FALSE, plotInflation=TRUE,
              calcFactor=FALSE)
out <- as.data.table(out)
colnames(out) <- c("L8S2", "All")

## save the stable period dynamic inflation factors
sapply(colnames(out), function(X){
  vec <- out[, get(X)]
  save(vec, file=
         paste0("data/trainingPars/train2_seasonalAdjustment", X, 
                ".Rdata"))
  return(print("Saved"))
})

##----------------------------------------------------------------------------##
# 2. Get static inflation factor based on increased variance in 
## monitoring residuals FOR LANDSCAPE APPLICATION ONLY. 
## The calculation of this step is dependent on the seasonally-adjusted 
## residuals, so make sure script 1 has been run with the "dynamic" 
## inflationType set.

## STOP. HAVE YOU RUN SCRIPT 1 WITH dynamic INFLATION YET?
runInflation(nRun=2, plotSeasons=FALSE, plotInflation=FALSE, calcFactor=TRUE)



