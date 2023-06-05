##########################################################
## Purpose: Remove NAs and duplicates
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: May 2022
##########################################################

## -----------------------------------------------------------------------------
## removeNAs = remove NAs from both std resids and obs dates
## removeDups = remove duplicate values by date by retaining the date with the 
##              highest-anomaly observation
## -----------------------------------------------------------------------------
removeNAs <- function(listObj){
  subDates <- listObj$obsDates[!is.na(listObj$stdResid)]
  outZ <- listObj$stdResid[!is.na(listObj$stdResid)]
  return(list(stdResid=outZ, obsDates=subDates))
}
removeDups <- function(stdRes, makeSumPlots=FALSE){
  
  # remove duplicates
  dupDates <- as.vector(unlist(sapply(stdRes, function(C){
    return(C[["obsDates"]])
  })))
  dupDates <- sort(dupDates[duplicated(dupDates)])
  
  for(dup in dupDates){
    # get the zscore values from the different sensors for the duplicated date
    rg <- lapply(stdRes, function(K){
      dat <- K$stdResid
      d <- K$obsDates
      absV <- abs(dat[d==dup])
      
      #if there are duplicated dates within a sensor's own ts
      if(length(absV[!is.na(absV)]) > 1){ 
        # #if there are exactly duplicated values for the duplicated dates
        if(which(duplicated(absV)) > 1){
          absVdupVals <- absV[absV==absV[duplicated(absV)]]
          absVnoVals <- absV[absV!=absV[duplicated(absV)]]
          
          # keep only the first instance, assign -9999 to the others
          if(length(absVnoVals)==0){
            absV <- c(absVdupVals[1], rep(-9999, length(absVdupVals)-1))
          } else {
            absV <- c(absVdupVals[1], absVnoVals)
          }
        }
        
        indexNA <- which(absV!=max(absV))
        indexDates <- which(d==dup)
        dat[indexDates][indexNA] <- NA
        d[indexDates][indexNA] <- NA
        
        return(list(absV=max(absV), dat=dat, d=d))
      } else if(length(absV[!is.na(absV)]) == 1){ # for dup dates across sensors
        return(list(absV=absV[!is.na(absV)], dat=dat, d=d))
      } else {
        return(list(absV=NA, dat=dat, d=d))
      }
    })
    
    #incorporate the duplicated NAs to the z-scores
    ## these are filtered out (with the dates) with the creation of outRow
    for(K in names(stdRes)){
      stdRes[[K]]$stdResid <- rg[[K]]$dat
    }
    
    rg <- sapply(rg, function(X) return(X$absV))
    
    # re-value NAs (if a sensor didn't have an observation then)
    rg[is.na(rg)] <- -9999
    
    # get the name of the sensor that did not have the max abs zscore value
    ## it is possible that some of the duplicated dates do not have actual 
    ## values (e.g. L8 is already filtered out from before, but we have a valid 
    ## S2 reading)
    ### i.e. this assigns NA to numeric(0), so nothing happens
    turnNA <- names(rg)[which(rg!=max(rg, na.rm=TRUE))]
    
    if(length(turnNA) == 1){
      stdRes[[turnNA]][["stdResid"]][stdRes[[turnNA]][["obsDates"]]==dup] <- NA
    } else if(length(turnNA) == 2){
      stdRes[[turnNA[1]]][["stdResid"]][stdRes[[turnNA[1]]][["obsDates"]]==dup] <- NA
      stdRes[[turnNA[2]]][["stdResid"]][stdRes[[turnNA[2]]][["obsDates"]]==dup] <- NA
    } else if(length(turnNA) ==3){
      stdRes[[turnNA[1]]][["stdResid"]][stdRes[[turnNA[1]]][["obsDates"]]==dup] <- NA
      stdRes[[turnNA[2]]][["stdResid"]][stdRes[[turnNA[2]]][["obsDates"]]==dup] <- NA
      stdRes[[turnNA[3]]][["stdResid"]][stdRes[[turnNA[3]]][["obsDates"]]==dup] <- NA
    }
  }
  
  # now if we remove NAs and combine sZ, we have a cohesive timeseries
  stdRes <- lapply(stdRes, removeNAs)
  
  allStdRes <- as.vector(unlist(sapply(stdRes, function(R){
    return(R[["stdResid"]])
  })))
  
  allObsDates <- as.vector(unlist(sapply(stdRes, function(R){
    return(R[["obsDates"]])
  })))
  
  if(makeSumPlots){
    # if(combined){
    #   sensorLabs <- lapply(names(stdRes), function(m){
    #     if(grepl("sentinel1", m)){
    #       sensNames <- rep("sentinel1", length(stdRes[[m]]$obsDates))
    #     } else {
    #       sensNames <- stdRes[[m]]$spikeVals$s
    #       sensDates <- stdRes[[m]]$spikeVals$t
    #       sensNames <- sensNames[sensDates %in% stdRes[[m]]$obsDates]
    #     }
    #     return(sensNames)
    #   })
    #   names(sensorLabs) <- names(stdRes)
    #   
    #   allSensors <- as.vector(unlist(sapply(names(sensorLabs), function(R){
    #     return(sensorLabs[[R]])})))
    # }
    
    allSensors <- as.vector(unlist(sapply(names(stdRes), function(R){
      return(c(rep(R, length(stdRes[[R]]$stdResid))))})))
    
    # need to order these before ordering the dates
    allSensors <- allSensors[order(allObsDates)]
  }
  
  # this reorders allZ to be in the same order as when allDates is sorted
  allStdRes <- allStdRes[order(allObsDates)]
  allObsDates <- allObsDates[order(allObsDates)]
  
  if(makeSumPlots){
    return(data.table(stdRes=allStdRes, obsDates=allObsDates, 
                      allSensors=allSensors))
  } else {
    return(data.table(stdRes=allStdRes, obsDates=allObsDates))
  }
}
