##########################################################
## Purpose: Calculate accuracy metrics and plot
## Run medium:
##  - Mac or PC
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
##########################################################
groundhog.library(data.table, groundhogDate)
groundhog.library(ggplot2, groundhogDate)
groundhog.library(ggdark, groundhogDate)
groundhog.library(viridis, groundhogDate)

bayes <- FALSE
source("scripts/args.R")
source("scripts/funs/trainE_metricsCalc.R")
source("scripts/funs/trainF_accuracyCV.R")

script <- 3
imgFolder <- "figures/"
source("scripts/argsTrain.R", local=TRUE)

##----------------------------------------------------------------------------##
# Plot Accuracy (F1 score, PR curve, and other metrics)
## Ch 1 Figures 5, 6, S5, S6 
consecObs <- 1
foldPath <- paste0("data/trainingPars/train5_allMetrics/")
f <- list.files(foldPath)
# f <- f[grepl(paste0("conObs", consecObs), f)]
metricsFiles <- f[grepl(paste0("allMetrics_conObs1_", windowType), f)]

figPars <- data.table(img=c("f1", "pr", "accLag", "boxLag", "f1Size"),
                      width=c(18,   18,      20,        20,       18),
                      height=c(15,  14,      15,        15,       12))
figPars[, `:=` (units="cm", res=350)]
                      
# splitNum = "at least" or "at most", 2nd argument is a percentage
results <- plotAccuracy(f=metricsFiles, foldPath, randomLam=FALSE, consecObs, 
                        plotF1=TRUE, plotPR=FALSE, plotMeta=TRUE, base, 
                        splitSize=FALSE, splitNum=c("at least", 1), 
                        plotWindow=60, windowType,
                        saveFig=FALSE, figPars,
                        imgFolder=imgFolder, bayes=bayes, probType="total")
