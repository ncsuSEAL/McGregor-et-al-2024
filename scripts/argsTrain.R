##########################################################
## Purpose: Define key arguments for functions
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: May 2022
##########################################################
### primary run parameters -------------------------------------------------####
longTime <- FALSE         # use data from L8 launch?
weights <- c(28, 0.8)     #scaling facfor for ewmaType != "ewma"
runType <- "train"

#Define the lambdas to iterate over
if(script==2){
  l <- c(0, 0.025, 0.05, 0.1, 0.2, 0.5, 1) #logistic spectrum
} else if(script==3){
  l <- seq(0.05, 1, by=0.05)  # final spectrum
}

### folders and files ------------------------------------------------------####
dataPath <- "data/"
sensorFold <- "sensorData" # folder w/SR band data, each sensor is other table
dataTrainPath <- "trainingPars/trainingDataPoints.csv" ## table of training data

### sensor arguments -------------------------------------------------------####
# sensors to use for analysis and their respective bands (here for NDVI and EVI2)
# sensKeep <- c("landsat8", "sentinel2", "sentinel1")
bandList <- list(
  "landsat8"=list("B4"=1, "B5"=2, "QA_PIXEL"=3, "QA_RADSAT"=4, "QA_AEROSOL"=5),
  "sentinel1"=list("VH"=2), #("VV" = 1, "VH" = 2)
  "sentinel2"=list("NDVI"=1, "SCL"=2)
)

################################################################################
# Format arguments

## Load metadata ---
base <- fread(file.path(dataPath, dataTrainPath))
base <- base[, `:=` (dateDist = as.Date(dateDist, format="%d/%m/%Y"),
                     datePre = as.Date(datePre, format="%d/%m/%Y"))]
pointids <- sort(unique(base[,pointid]))

if(script != 3){
  ## Format bandlist and sensor data ---
  bandList <- bandList[grepl(gsub(", ", "|", sensKeep), names(bandList))]
  sensorFiles <- list.files(file.path(dataPath, sensorFold))
  sensorFiles <- sensorFiles[grepl(paste0(names(bandList), collapse="|"), 
                                  sensorFiles)]
  sensorNames <- gsub(".csv", "", sensorFiles)

  names(bandList) <- sensorNames

  ## Define dateRange ---
  if(longTime){
    dateRange <- c(as.Date("2013-02-01"), as.Date("2020-01-31"))
  } else {
    dateRange <- c(as.Date("2015-06-23"), as.Date("2020-01-31"))
  }

  dates <- as.numeric(seq(dateRange[1], dateRange[2], by=1))
}