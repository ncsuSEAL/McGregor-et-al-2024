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
groundhog.library(data.table, groundhogDate)
library(parallel)
source("scripts/funs/trainD_metricsPrep.R")

script <- 3
combos <- c("sentinel1", "landsat8, sentinel2", 
            "landsat8, sentinel2, sentinel1")
fileNames <- c("S1", "L8S2", "All")
sensKeep <- combos[2] # doesn't matter, both L8S2 and All are run below
bayes <- FALSE # only for chapter 2
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
if(runAllThresholds){
       ## all lambdas, all points, both sensors, all thresholds
       ### parallel points
       threshSeq <- seq(0.1, 1, by=0.1)
       runAll(threshSeq, base, dates, probThresh, window, windowType, bayes, 
              saveFile, pointids, consecObs, dataPath, fileOutSave)
} else {
       ## all lambdas, all points, both sensors, single threshold
       
       for(obsNumber in c(1,3)){
              print(paste0("Running for ", obsNumber, 
                            " consecutive observation(s)"))
              runAll(threshSeq=NULL, base, dates, probThresh=0.5, window, 
                     windowType, bayes, saveFile, pointids, consecObs=1, 
                     dataPath, fileOutSave)
       }
}
##---------------
if(makePlots){
       ## Plot of median lags (lagFast) comparing All to L8S2
       fileLoad <- gsub("Metrics_", "MetricsThresholds_", fileOutSave)
       if(bayes) fileLoad <- gsub("normal", "bayes", fileLoad)
       test <- fread(fileLoad)
       bl <- test[, median(lagFast, na.rm=TRUE), by=.(thresh, sensor, lambda)]
       blSlow <- test[, median(lagSlow, na.rm=TRUE), by=.(thresh, sensor, lambda)]

       library(ggplot2)
       bl[, thresh := as.character(thresh)]
       blSlow[, thresh := as.character(thresh)]
       ggplot(bl, aes(x=lambda, y=V1, group=thresh)) +
       geom_path(aes(color=thresh)) +
       geom_path(data=blSlow, aes(x=lambda, y=V1, group=thresh, color=thresh), 
              linetype="dotted") +
       # scale_color_manual(values=viridis::magma(11)) +
       ylab("Median lag (days)") +
       facet_wrap(~sensor)
}

