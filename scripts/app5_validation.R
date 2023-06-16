##########################################################
## Purpose: Create table of validation pixels to look over, then calculate
##          validation accuracy metrics
## Input: all necessary functions
## Variable created: csv of validation pixels, and metrics
## Run medium:
##  - PC for creating the point tables, either for making accuracy metrics
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.3, Dec 2022
## Last updated: Dec 2022
##########################################################
library(data.table)
library(terra)
library(MuMIn)
# setwd("Z:/IanMcGregor/")

source("scripts/funs/appC_analyzeLandscape.R")
source("scripts/funs/appE_validation.R")
source("scripts/funs/commonB_modZ.R")
source("scripts/funs/commonC_dupNA.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")

# 1. Create table of validation points
allCells <- rbindlist(lapply(c(1,5), getValidationPix, numValidate=30))
fwrite(allCells, "data/validationCells1.csv") #number added to not overwrite

# 2. Analyze validation results
library(data.table)
val <- fread("data/validationCells.csv")
val[, .N, by=.(region, val)]

## true positive rate
val[region==1, .N, by=val][val==1, N] / (nrow(val) / 2)
val[region==5, .N, by=val][val==1, N] / (nrow(val) / 2)

unsure <- val[val==2, ]

## which are affected by missed cloud?
val[, missReason := ifelse(grepl("missed cloud", notes), "cloud",
                           ifelse(grepl("deciduous", notes), "dryness", 
                                  "other"))]

val[val != 1, .N, by=.(missReason, region)]
