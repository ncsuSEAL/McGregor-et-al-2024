##########################################################
## Purpose: Perform K-fold cross validation
## Run medium:
##  - Mac or PC
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
##########################################################
groundhog.library(data.table, groundhogDate)
groundhog.library(ggplot2, groundhogDate)

# Step 0: Load the necessary variables and data
sensKeep <- c("landsat8, sentinel2")
fileNames <- c("L8S2") # S1, L8S2, or All
script <- 1
bayes <- FALSE
source("scripts/args.R", local=TRUE)
source("scripts/argsTrain.R", local=TRUE)
source("scripts/funs/commonD_ewmaBinaryCat.R")
source("scripts/funs/commonF_probBayes.R")
source("scripts/funs/trainC_byLambda.R")
source("scripts/funs/trainD_metricsPrep.R")
source("scripts/funs/trainE_metricsCalc.R")
source("scripts/funs/trainF_accuracyCV.R")
source("scripts/funs/trainG_densFuns.R")

load(paste0(dataPath, "/trainingPars/train1_", fileNames, ".Rdata"))

runCrossVal <- function(pointids, nFolds, lambda, window, runType,
                        windowType, probType, sensKeep){
  # Step 1: Randomly split data into k-subsets
  densFun <- NULL #we make the functions new each time
  
  ## only want to define CV folds based on true disturbances
  ptDist <- pointids[!grepl("a", pointids)]
  ptNoDist <- pointids[grepl("a", pointids)]
  
  # this may produce warning if length(ptDist) is not completely
  ## divisible by nFolds. This is fine, because the function automatically
  ## accounts for this.
  folds <- split(ptDist, 1:nFolds)
  
  # Step 2: Define function: get outputs using each of the folds as a 
  ## testing data set
  fullLength <- length(folds)+1
  comp <- rbindlist(lapply(1:fullLength, runFolds, folds, lambda, window, 
                           pointids, runType, windowType, consecObs=1, 
                           probType, densFun))
  
  # Step 3 (deprecated): Calculate overall k-fold cross validation error
  # allErrors <- sapply(1:11, calcErrorCV)
  # plot(1:11, c(allErrors[1:10], ""), pch=16, xlab="Fold", ylab="Error")
  # points(11, allErrors[11], pch=16, col="red")
  
  # Step 4: Calculate metrics for plotting
  colN <- paste0("stat", 1:window)
  met <- rbindlist(lapply(unique(comp$iter), outputF1, comp, colN, window))
  return(met)
}

# Get cross validation data
## training and testing data is determined based on the number of folds. So if
## we have 10 folds, then we train on 
## 100-((length(ptDist)/10)/length(ptDist))*100 % of the data and test on the 
## other ((length(ptDist)/10)/length(ptDist))*100 %
set.seed(23469)
cv <- runCrossVal(pointids, nFolds=10, lambda, window, runType, windowType,
                  probType="total", sensKeep)

# Visualize the results
dt <- melt(cv, id.vars=c("iter", "type"))
dt <- dt[variable != "overallAcc"]
cols <- c(viridis::viridis(10))
cols <- colorRampPalette(colors = 
                          c("purple", "blue", "red", "orange", "#996633"))(10)

## Chapter 1, Figure S8: boxplots of the different variables
dt[, iter := as.numeric(iter)
][, iter := ifelse(iter <=10, paste0("Fold", iter), "Full data")
][, iter := factor(iter, levels=c(paste0("Fold", 1:10), "Full data"))]
dt[, variable := ifelse(variable=="f1", "F1",
                        ifelse(variable=="recall", "Recall", "Precision"))]

# png("figures/figXX_crossVal.png", width=20, height=14, units="cm", res=350)
ggplot(dt) +
  aes(x = iter, y = value, fill = iter) +
  geom_boxplot(alpha=0.7) +
  # scale_fill_manual(values=c(cols, "violet")) +
  scale_fill_manual(values=c(cols, "grey")) +
  theme_minimal() +
  facet_wrap(~variable) +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title=element_blank())
# dev.off()

##############################################################################
# other plots

if(otherPlots){
  ## overlay plot just showing the full data and subset
  ggplot(dt) +
    aes(rep(rep(1:window, 11), 3), value) +
    geom_point(aes(col=type)) +
    xlab("Days since disturbance") +
    facet_wrap(~variable)

  ## overlay plot of full data and subsets, different color per iteration
  dt[, iter := factor(iter, levels=as.character(1:11))]
  dtsub <- dt[iter != "11", ]
  dtfull <- dt[iter == "11", ]

  ggplot() +
    geom_path(data=dtsub, aes(rep(rep(1:window, 10), 3), value, col=iter)) +
    geom_path(data=dtfull, aes(rep(1:window, 3), value, col=iter), linewidth=2) +
    scale_color_manual(values=c(cols, "black")) +
    xlab("Lag (daily)") +
    facet_wrap(~variable)
}