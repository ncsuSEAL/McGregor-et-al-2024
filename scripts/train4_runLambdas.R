##########################################################
## Purpose: Run iterations of lambda to calculate zEWMA and (Bayes) probability
##          of disturbance. Output is a separate csv per lambda.
## Run medium:
##  - Mac
##
## Figures for McGregor et al 2023: None
##
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: Feb 2023
##########################################################
library(data.table)
library(parallel)
library(MuMIn)

source("scripts/funs/trainC_byLambda.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")
source("scripts/funs/commonF_probBayes.R")
funs <- ls()

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

##----------------------------------------------------------------------------##
# Step 1: Run function only for L8S2 and All for each specific window.
## NOTE - this is only for a sensitivity analysis looking over different
## monitoring windows.
modStats <- rbindlist(lapply(2:3, runCombos, funs=funs, consecObs=1, 
                             returnOnlyMod=TRUE, bayes=FALSE, 
                             probType="log", returnOnlyMetrics=FALSE))
fwrite(modStats, paste0("data/dataMyanmar/logModStats", window, "days.csv"))
##----------------------------------------------------------------------------##
# Step 2: Run over everything to move on to train5.
## Run on Mac or PC
## This creates a csv per lambda with zEWMA, probs, and probsBayes. Because it
### includes everything, bayes should always be set to TRUE here.
### ALSO - saved file does not matter if windowType=obs or days (args.R) 
### because they are handled separately in train5
modStats <- sapply(c(2:3), runCombos, funs=funs, bayes=TRUE, probType="total")

## Run in HPC
runSensCombo(nRun=args[1], funs=funs, returnOnlyMod=FALSE)

