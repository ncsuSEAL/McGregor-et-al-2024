##########################################################
## Purpose: Run Sen2Cor conversions in HPC
##          (companion script = S2_createOrders.R)
## Run medium:
##  - HPC job. See data6_runSen2Cor_batchAuto/Manual.csh
## Creators: Ian McGregor, imcgreg@ncsu.edu
## Editors: Izzi Hinks, Xiaojie Gao
## System: R Version 3.6.3, July 2021, updated Nov 22
##########################################################
library(data.table)
library(parallel)

# setwd("Z:/IanMcGregor/")

## if want to run each job separately, use
### sen2cor_rmpi/data6_runSen2Cor_batchManual.csh

## if want to run all the jobs automatically, then uncomment below
### and use sen2cor_rmpi/data6_runSen2Cor_batchAuto.csh

source("./sen2cor_rmpi/data6_funs.R")

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

convertAllS2(nRun=args[1], rmpi=FALSE)

# convertAllS2(nRun=1, rmpi=FALSE)
