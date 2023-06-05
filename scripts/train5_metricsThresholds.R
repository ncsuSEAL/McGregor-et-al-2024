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
bayes <- TRUE
source("scripts/args.R")
source("scripts/argsTrain.R", local=TRUE)
consecObs <- 1
filePath <- "/trainingPars/train4_fullZewmaProbs/train4_S_N.csv"
saveFile <- paste0(dataPath, filePath)
fileOutSave <- paste0(dataPath, "/trainingPars/train5_allMetrics/",
                      "train5_normal_allMetrics_conObs", consecObs, "_", 
                      windowType,".csv")
##---------------
# Ch 1
## all lambdas, all points, both sensors, single threshold = 129 s (par pts)
runAll(threshSeq=NULL, base, dates, probThresh=0.5, window, windowType, bayes, 
       saveFile, pointids, consecObs, dataPath, fileOutSave)

# Ch 1 and 2
## all lambdas, all points, both sensors, all thresholds = 1283 s
### parallel points
threshSeq <- seq(0.1, 1, by=0.1)
runAll(threshSeq, base, dates, probThresh, window, windowType, bayes, 
       saveFile, pointids, consecObs, dataPath, fileOutSave)

##---------------
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

###############################################################################
# For these iterations, need to pre-define some sample variables
dt <- fread(paste0(dataPath, "
                   trainingPars/fullZewmaProbs/train4_All_conObs1_days_5.csv"))
lambda <- 0.05
pt <- "g0p0_1"
probThresh <- 0.2
consecObs=1
bayes=TRUE
dt[, probType := NULL]

# single lambda, single point, single sensor, single threshold = 0.019 s
system.time(retrieveMetrics(pt, dt, consecObs, probThresh, window, windowType, 
                            bayes, lambda))

# single lambda, all points, single sensor, single threshold = 9.376 s
system.time(rbindlist(lapply(pointids, retrieveMetrics, dt, consecObs, 
                             probThresh, window, windowType, bayes, lambda)))

# all lambdas, all points, single sensor, single threshold = 61.599 s (parallel)
cl <- setPar(base, dates, probThresh, window, windowType, bayes, saveFile, 
             pointids, consecObs)
out <- rbindlist(lapply(l, allLamAllPt, saveFile, probThresh, par=TRUE, cl=cl))
stopCluster(cl)
