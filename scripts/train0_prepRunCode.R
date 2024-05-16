# install.packages("groundhog")
# library(groundhog)
# packageList <- c("data.table", "MuMIn", "terra", "moments", "RColorBrewer",
#                 "ggplot2", "ggpointdensity", "ggpubr", "ggdark", "viridis", 
#                 "cowplot")
# groundhog.library(packageList, "2023-07-31")

##------------------------------------------------------------##
## Please see the workflow for the corresponding step descriptions
library(groundhog)

# the only necessary folders are scripts, data/sensorData, and 
# data/trainingPars/trainingDataPoints.csv

## the date used for loading the package versions
#! if this date is changed from 2023-07-31, then need to re-run the lines above
groundhogDate <- "2023-07-31"

runAnalysis <- function(groundhogDate){
    makePlots <- FALSE
    if(!dir.exists("figures")) dir.create("figures")

    # Step 3-1
    print("Step 1/: Calculating initial models and residuals.")
    inflationType <- "static"
    source("scripts/train1_getResiduals.R", local=TRUE)

    # Step 3-2
    print("Step 2/: Identifying static inflation.")
    source("scripts/train2_sdInflation.R", local=TRUE)

    # Step 3-3
    print("Step 3/: Getting new residuals with static inflation.")
    inflationType <- "dynamic"
    source("scripts/train1_getResiduals.R", local=TRUE)

    # Step 3-4
    print("Step 4/: Calculating final inflation factor for landscape run.")
    source("scripts/train2_sdInflation.R", local=TRUE)

    # Step 3-5
    print("Step 5/: Creating probability density functions for the z-scores.")
    source("scripts/train3_createProbFuns.R", local=TRUE)

    # Step 3-6
    print(paste0("Step 6/: Finding zEWMA values and probabilities of ",
                "disturbance across thresholds."))
    source("scripts/train4_runLambdas.R", local=TRUE)

    # Step 3-7
    print("Step 7/: Prepping metric calculation across all probability thresholds.")
    runAllThresholds <- TRUE
    source("scripts/train5_metricsThresholds.R", local=TRUE)

    # Step 3-8
    print(paste0("Step 8/: Calculating accuracy metrics for all probability ",
                "thresholds. Outputs are used to choose a specific threshold."))
    source("scripts/train6_accuracy.R", local=TRUE)

    # Step 3-9
    print("Step 9/: Re-running metric prep using chosen threshold.")
    runAllThresholds <- FALSE
    source("scripts/train5_metricsThresholds.R", local=TRUE)

    # Step 3-10
    print(paste0("Step 10/: Visualize F1 scores, FPR, FNR, and lag differences",
                " across lambdas for the chosen probability threshold."))
    source("scripts/train7_plotAccuracy.R", local=TRUE)

    # Step 3-11
    print("Running k-fold cross-validation for chosen lambda")
    otherPlots <- FALSE
    source("scripts/train8_crossVal.R", local=TRUE)

    print(paste0("End of initial analysis. Please refer ", 
                "back to the workflow document for other steps."))
    print(paste0("Additional plots can be made in ",
                "`train7_plotOthers.R`, and you can change the number ",
                "of consecutive observations needed to detect a ", 
                "disturbance (e.g. for the ECDF plot)."))
}

runAnalysis(groundhogDate)