##########################################################
## Purpose: Calculate full accuracy metrics and lags, and save files
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
library(ggplot2)
script <- 3
path <- "data/trainingPars/"
fileLoc <- "train5_S_allMetricsThresholds_conObs1_days.csv"
bayes <- FALSE
source("scripts/funs/trainE_metricsCalc.R")
source("scripts/args.R")
source("scripts/argsTrain.R", local=TRUE) #error is fine, ignore it

## create and save metrics (only need to do once)
### type = bayes or normal
processMetrics(path, type="normal", fileLoc, lambdaPlot=l, window)

################################################################################
# Chapter 1, Fig 4 - how F1 scores + lags change with lambda, sensors, threshold
library(ggpubr)
library(viridis)
plotMetrics <- function(type, consecObs, windowType){
  pal <- colorRampPalette(colors = 
                            c("purple", "blue", "red", "orange", "#996633"))(10)
  for(i in c("f1", "lags")){
    if(i=="f1") f <- i else f <- "medianLags"
    
    dtFull <- lapply(type, function(X){
      dt <- fread(paste0(path, "train6_", X, "_", f, "All_conObs", consecObs, 
                         "_", windowType, ".csv"))
      dt[, type := X]
      return(dt)
    })
    
    dtFull <- rbindlist(dtFull)
    dtFull[, sensorNew := ifelse(sensor=="All", "MM", "MO")]
    
    if(i=="f1"){
      test <- dtFull[, median(f1), by=.(type, threshold, sensorNew, lambda)]
      test[, threshold := as.character(threshold)]
      gF1 <- ggplot(test, aes(x=lambda, y=V1)) +
        geom_path(aes(group=threshold, color=threshold)) +
        ylab("F1 scores") +
        xlab("") +
        ylim(0, 1) +
        # geom_hline(yintercept=0.9, linetype="dashed") +
        # geom_hline(yintercept=0.95, linetype="dashed") +
        scale_color_manual(values=pal) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
      
      
      if(length(type)==1){
        gF1 <- gF1 + 
          facet_wrap(~sensorNew)
      } else {
        gF1 <- gF1 + 
          facet_wrap(~type + sensorNew, dir="v")
      }
    } else {
      test <- dtFull
      test[, threshold := as.character(threshold)]
      gLagF <- ggplot(test, aes(x=lambda, y=medLagFast)) +
        geom_path(aes(group=threshold, color=threshold)) +
        ylim(0, 30) +
        xlab("") +
        ylab("Fastest detection (days)") +
        scale_color_manual(values=pal) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
      
      gLagS <- ggplot(test, aes(x=lambda, y=medLagSlow)) +
        geom_path(aes(group=threshold, color=threshold)) +
        ylim(0, 30) +
        ylab("Slowest detection (days)") +
        scale_color_manual(values=pal) +
        xlab(expression(lambda))
      
      if(length(type)==1){
        gLagF <- gLagF + facet_wrap(~sensorNew)
        gLagS <- gLagS + facet_wrap(~sensorNew)
      } else {
        gLagF <- gLagF + facet_wrap(~type + sensorNew, dir="v")
        gLagS <- gLagS + facet_wrap(~type + sensorNew, dir="v")
      }
    }
  }
  
  g <- ggarrange(gF1, gLagF, gLagS, labels=c("a", "b", "c"), common.legend=TRUE,
                 legend="right", nrow=3)
  
  fileOut <- "figures/figXX_thresholds.png"
  
  png(fileOut, width=20, height=18, units="cm", res=350)
  print(g)
  dev.off()
}
type <- c("normal", "bayes")
plotMetrics(type="normal", consecObs=1, windowType)
