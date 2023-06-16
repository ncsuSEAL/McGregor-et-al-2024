##########################################################
## Purpose: Extract band values at each training point per 
##          each converted Sentinel 2 L2A image and
##          write to table with format matching GEE output
## Run medium:
##  - HPC job, via csh file of the same name in SEAL/Ian/
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 3.6.3, April 2021
##########################################################
library(data.table)
library(parallel)
library(terra)

# NB PLEASE READ ALL NOTES IN THE COMPANION SCRIPT for explanations and 
## supplementary code
# ALSO note depending on the node given, this takes either >24 hours or <2 hours
## for all points. Make sure to use BSUB -R select[avx2] in the csh file.
## We could speed this up by submitting multiple jobs but we would then have to
## use rMPI, which is a hassle.

# Otherwise, this script contains the following workflow.
## 1: Run the main function to get Sentinel-2 band values.

# Step 1: Define the variables necessary for the main function
# base <- "/Volumes/SEAL/IanMcGregor"
# base <- "//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor"
# base <- "/rsstu/users/j/jmgray2/SEAL/IanMcGregor"

# STOP. Have you made sure that "data/myanmar/trainingDataPoints.csv" AND
## "data/myanmar/trainingDataS2Tiles.csv" have been updated??

source("dissertation/scripts/funs/dataA_extractS2.R")

trainPath = "dissertation/data/myanmar/trainingPars/trainingDataPoints.csv"
tilePath <- "dissertation/data/myanmar/trainingPars/trainingDataS2Tiles.csv"
xcoord <- "coordX"
ycoord <- "coordY"
prev <- NULL
tileStr <- "all" #can do a c() of specific tiles or just "all"
L2APath <- "s2Data/L2A"
outFile <- "sentinel2l2a_14Nov22.csv"
txtFile <- "dataA_extractBandsProgress.txt"
if(!file.exists(txtFile)) file.create(txtFile)

#-------------------------------------------------------------------------------

# Step 2: 
extractBands(prev, tileStr, L2APath, txtFile, 
            trainPath, tilePath, removeNArows=TRUE, outFile)
