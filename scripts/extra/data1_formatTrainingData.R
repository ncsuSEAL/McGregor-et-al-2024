##########################################################
## Purpose: Format training data from google forms, create pointids, and
##          write out as the primary training data set used for the remainder
##          of the research workflow
## Run medium:
##  - Either PC or Mac
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 3.6.3, March-April 2021
## Last updated: July 2022
##########################################################
library(data.table)

formatData <- function(dataPath){
  # Step 1: Read in downloaded google forms csv
  pl <- fread(dataPath)
  
  # Step 2: Create location, dayDiff, pointid
  pl <- pl[!is.na(coordX), 
  ][, `:=` (dateDist = as.Date(dateDist, format="%m/%d/%Y"), 
            datePre = as.Date(datePre, format="%m/%d/%Y"),
            datePlanetDist = as.Date(dateDist, format="%m/%d/%Y"),
            datePlanetPre = as.Date(dateDist, format="%m/%d/%Y"))
  ][, `:=` (dayDiff = as.numeric(dateDist - datePre),
            location="myanmar")]
  
  #have to do this separately bc filtered out NA coordinates above
  pl <- pl[, pointid := paste0("g", nGroup, "p", nPoly, "_", nPoint)]
  pl <- pl[order(nGroup, nPoly, nPoint), ]
  
  ## remove Peru points
  pl <- pl[!grepl("200", pointid), ]
  
  ## double check the recorded dist and pre dates are correct
  testDiff <- pl[dayDiff < 0, pointid]
  if(length(testDiff) > 0){
    stop(paste0("Negative dayDiff for pointid ", testDiff))
  }
  
  #double check if you have duplicated points / coordinates at all and filter
  sameP <- pl[duplicated(pointid), ]
  sameC <- pl[duplicated(coordX), ]
  
  if(nrow(sameP) > 0 | nrow(sameC) > 0){
    stop(paste0("You have duplicated points with ", sameP$pointid, 
                " (same pointid) and ", sameC$pointid, " (same coords)"))
  }
  
  # Step 3: Give the undist points the same "disturbance date" as their counterpart
  undistPts <- pl$pointid
  undistPts <- undistPts[grepl("a", undistPts)]
  
  undistData <- rbindlist(lapply(undistPts, function(X){
    sub <- pl[pointid==gsub("a", "", X), ]
    undist <- pl[pointid==X, ]
    undist[, `:=` (dateDist=sub$dateDist, datePre=sub$datePre)]
    return(undist)
  }))
  
  full <- rbind(pl[!grepl("a", pointid)], undistData)
  full <- full[order(pointid), ]
  
  return(full)
}

# Format data, then write out to new csv
full <- formatData("data/dataMyanmar/trainingDataPointsNew.csv")
fwrite(full, "data/dataMyanmar/trainingDataPoints.csv")

## quick plot of training point stats
layout(matrix(1:4, nrow=2))
hist(full$distPercS2, main="% Dist S2")
hist(full$distPercL8, main="% Dist L8")
hist(full$dayDiff, main="Temporal res dist", breaks=30)
