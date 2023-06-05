##########################################################
## Purpose: Define key arguments for functions
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: May 2022
##########################################################
### primary run parameters -------------------------------------------------####
window <- 60             #time window for detecting disturbances
windowType <- "days"      #what kind of window to use for calculating lags? This
                          # is also used for the accuracy vs lags plots. 
                          # The options are "days" or "obs"
ewmaType <- "ewma"        #different ways to calculate ewma: ewma, chart, filter
staticInflation <- 1.2   #static inflation factor for landscape application

if(bayes){
  probThresh <- 0.10
  lambda <- 0.3
} else {
  probThresh <- 0.50     #threshold probability for indicating disturbance?
  lambda <- 0.6          #scaling factor for ewma
}

### sensor arguments -------------------------------------------------------####
vegIndex="ndvi"           # ndvi or evi2 for now (as of late Sep 2022)
pol="vh"                  #polarization for Sentinel-1

### spike filter arguments -------------------------------------------------####
spikeWin <- 5             #window (obs) for spike filter
spikeAmp <- 3             #amplitude for spike filter
spikeThresh <- 3          #threshold for spike filter
spikeTime <- 365          #max time for spike filter to consider
spike2 <- TRUE

### output arguments -------------------------------------------------------####
returnStats <- FALSE      #return the model stats of each sensor's ts model
saveStats <- FALSE        #save the model stats of each sensor's ts model?
returnHarMod <- TRUE        #return harmonic models for each sensor?
makeSumPlots <- TRUE      #retain information needed to make diagnostic plots?
saveProbModel <- FALSE    #save the overall logistic model as a model object?
probModelFile <- "data/trainingPars/probModel.RDS" 
plotProbModel <- FALSE     #plot the logistic model after calculation
pathPlanet <- "Z:/IanMcGregor/planetImagesNew/planetscope"