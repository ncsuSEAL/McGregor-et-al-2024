##########################################################
## Purpose: Calculate density functions for probability estimates
## Creator: Ian McGregor, imcgreg@ncsu.edu
## Run medium:
##        - Either Mac or PC, this is run on training data
## System: R Version 4.2.2, Feb 2023
## Last modified: Mar 2023
##########################################################
library(data.table)
library(parallel)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(moments)

source("scripts/funs/trainC_byLambda.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")
source("scripts/funs/commonF_probBayes.R")
source("scripts/funs/trainG_densFuns.R")
funs <- ls()

# Likelihood (run with training data)
folderPath <- "data/trainingPars/train3_densFun"
if(!dir.exists(folderPath)) dir.create(folderPath)
filePath <- paste0(folderPath, "/fun.Rdata")

chapterPars <- 1
if(chapterPars==1){
  testNormal <- FALSE
  createFuns <- TRUE
  compFigure <- FALSE
} else {
  testNormal <- FALSE
  createFuns <- FALSE
  compFigure <- TRUE # Ch 2 Fig 9 and S2
}

## note that bayes argument is only for when compFigure=TRUE
sapply(c(2:3), getLikelihood, bayes=FALSE, ewmaDays=TRUE, filePath, funs,
       testNormal, createFuns, compFigure)


