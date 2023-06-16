##########################################################
## Purpose: Prep accuracy metrics for either a single threshold, or run a 
##          sensitivity analysis over a sequence of thresholds
##
##          NB: This script is partitioned into sections to test different
##              run times of the workflow in order to understand the bottleneck.
##              The two iterations that are most in use are at the top.
## Run medium:
##  - Mac or PC
##
## Figures for McGregor et al 2023: None
##
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, April 2023
## Last modified: Apr 2023
##########################################################
library(data.table)
library(parallel)
source("scripts/funs/trainD_metricsPrep.R")

script <- 3
combos <- c("sentinel1", "landsat8, sentinel2", 
            "landsat8, sentinel2, sentinel1")
fileNames <- c("S1", "L8S2", "All")
sensKeep <- combos[2] # doesn't matter, both L8S2 and All are run below
bayes <- FALSE
source("scripts/args.R")
source("scripts/argsTrain.R", local=TRUE)
consecObs <- 1
filePath <- "/trainingPars/train4_fullZewmaProbs/train4_S_N.csv"
saveFile <- paste0(dataPath, filePath)

folderPath <- paste0(dataPath, "/trainingPars/train5_allMetrics")
if(!dir.exists(folderPath)) dir.create(folderPath)
fileOutSave <- paste0(folderPath, "/",
                      "train5_normal_allMetrics_conObs", consecObs, "_", 
                      windowType,".csv")
##---------------
# Ch 1
## all lambdas, all points, both sensors, single threshold = 129 s (par pts)
runAll(threshSeq=NULL, base, dates, probThresh=0.5, window, windowType, bayes, 
       saveFile, pointids, consecObs, dataPath, fileOutSave)

## all lambdas, all points, both sensors, all thresholds = 1283 s
### parallel points
threshSeq <- seq(0.1, 1, by=0.1)
runAll(threshSeq, base, dates, probThresh, window, windowType, bayes, 
       saveFile, pointids, consecObs, dataPath, fileOutSave)
