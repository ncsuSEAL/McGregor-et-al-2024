##########################################################
## Purpose: Functions for calculating sd inflation for z-scores (and variance)
##          i.e. look at seasonality of raw residuals over the year
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, October 2022
## Last modified: October 2022
##########################################################

## -----------------------------------------------------------------------------
## Function name and purpose:
##
## calcVar = calculate variables (raw residuals, variance factor diff)
## calcVarWrap = do calcVar over the assigned sensors
## rollMeanDays = do centered rolling average over days, allowing for duplicates
##                e.g. for a window of 45, this includes the mean of days 1-45 
##                even if one of those days has more than one value with it.
## addMissingDays = fill in missing days in chronological order, add NA for 
##                  associated value
## calcPlotRollMean = wrapper for `addMissingDays` and `rollMeanDays` to plot
## boxNDVI = create boxplots of NDVI or backscatter
## applySensor = wrapper for `calcPlotRollMean` or `boxNDVI` for each sensor,
##              depending on analysis
## inflationFactor = full wrapper of everything above
## -----------------------------------------------------------------------------

calcVarNDVI <- function(sens, pt, dataList, timeframe, dates){
  ts <- dataList[[pt]]$ptData[[sens]]$ts
  mod <- dataList[[pt]]$stdRes[[sens]]$model
  spikes <- dataList[[pt]]$stdRes[[sens]]$spikes
  obsDates <- dataList[[pt]]$stdRes[[sens]]$obsDates
  
  if(is.null(mod)){
    return(data.table(val=NA, dates=NA, sensor=sens, type=NA))
  } else {
    # this is only used when looking at the stable period, because we're looking 
    ## across all training points. 
    dateMod <- obsDates[(length(mod$model$y_mod)+1)]
    
    if(grepl("sentinel1", sens)){
      d <- dataList[[pt]]$ptData[[sens]]$dir
      d <- d[!is.na(ts)]
    }
    
    predVal <- predict(mod, newdata=data.frame(t_mod=dates))
    predVal <- predVal[!is.na(ts)]
    ts <- ts[!is.na(ts)]
    
    if(grepl("sentinel1", sens)){
      if(grepl("ASC", sens)){
        dir <- "ASCENDING"
      }   else {
        dir <- "DESCENDING"
      }
      
      predVal <- predVal[ts != 0 & d == dir]
      ts <- ts[ts != 0 & d == dir]
    }
    
    # NB: ONLY obsDates is already post-spike filter
    if(any(spikes %in% c(4))){
      predVal <- predVal[spikes != 4]
      ts <- ts[spikes != 4]
    }
    
    # calc raw residuals, then calc variance of pre2019 and 2019
    ## 2019-01-01 = 17897 | 2019-12-31 = 18261
    
    val1 <- ts
    type <- ifelse(grepl("sentinel1", sens), "backscatter", "ndvi")
    
    if(timeframe=="stable"){
      datesSub <- obsDates[obsDates < dateMod]
      val1 <- val1[obsDates < dateMod]
    } else if(timeframe=="monitor"){
      datesSub <- obsDates[obsDates >= 17897]
      val1 <- val1[obsDates >= 17897]
    }
    
    return(data.table(val=val1, dates=datesSub, sensor=sens, type=type))
  }
}
calcVarRes <- function(pt, dataList, variance, resDates, timeframe){
  agg <- dataList[[pt]]$oneTS
  type <- "resid"
  
  ## When we're looking at variance or specifically the monitoring period, we 
  ## use a static (2019-01-01) because 1) we're comparing to the landscape 
  ## application and 2) we're only looking at non-disturbed training points.
  dateMod <- 17897
  
  if(variance){
    # val2 <- var(agg[obsDates < dateMod, stdRes])
    # val3 <- var(agg[obsDates >= dateMod, stdRes])
    # out <- data.table(pointid=pt, varFact=val3/val2)
    agg[, period := ifelse(obsDates < dateMod, "stable", "monitor")]
    agg[, pointid := pt]
    return(agg)
  }
  if(resDates){
    if(timeframe=="stable"){
      out <- agg[obsDates < dateMod, ]
    } else if(timeframe=="monitor"){
      out <- agg[obsDates >= dateMod, ]
    }
    colnames(out) <- c("val", "dates", "sensor")
    out[, `:=` (type = type)]
    return(out)
  }
}
calcVarWrap <- function(pt, sensNames, dataList, variance, resDates, ndvi,
                        timeframe, dates){
  if(ndvi){
    out <- lapply(sensNames, calcVarNDVI, pt, dataList, timeframe, dates)
    return(rbindlist(out))
  }
  
  if(variance | resDates){
    out <- calcVarRes(pt, dataList, variance, resDates, timeframe)
    return(out)
  }
}
rollMeanDays <- function(tab, n){
  ## add on a prefix because this is an annual count, so allow for continuity
  prefix <- tab[doy >= (365-floor(n/2)+1),][, doy := doy-365]
  suffix <- tab[doy <= (floor(n/2)+2)][, doy := doy + (365-1)]
  tab <- rbind(prefix, tab, suffix)
  tab[, doy := doy + (floor(n/2))]
  
  out <- lapply(1:max(tab$doy, na.rm=TRUE), function(doyNum){
    if(doyNum < (floor(n/2)+1)){
      valDay <- NA
    } else {
      start <- (doyNum - floor(n/2))
      end <- (doyNum + floor(n/2))
      valDay <- mean(tab[doy >= start & doy <= end, val], na.rm=TRUE)
    }
    return(data.table(valDay=valDay, doy=doyNum))
  })
  return(rbindlist(out))
}
addMissingDays <- function(subDT, mainAdd){
  ## bring in full doy sequence
  ### this shouldn't need to be in sapply but having it separate affects the 
  ### data table itself for every iteration, so putting it here to get consistent
  ### results
  subDT <- subDT[!is.na(val)][order(dates)
  ][, doy := yday(as.Date(dates, origin=as.Date("1970-01-01")))]
  doyFull <- seq(1:365)
  subDT <- rbind(subDT,
                 data.table(val=NA, dates=NA, 
                            sensor=ifelse(grepl("combined", mainAdd), "both",
                                          unique(subDT$sensor)), type=NA,
                            doy=setdiff(doyFull, unique(subDT$doy))))
  subDT <- subDT[order(doy)]
  return(subDT)
}
calcPlotRollMean <- function(X, s, subDT, mainAdd, returnData, timeframe){
  if(timeframe=="monitor"){
    main <- ""
    xaxt <- "s"
  } else {
    main <- mainAdd
    xaxt <- "n"
  }
  
  if(s==2){
    ylab <- "Mean absolute z-score"
    yaxt <- "s"
  } else {
    ylab <- ""
    yaxt="n"
  }
  
  if(s==2 & timeframe=="stable")  par(mar=c(2,4,2,0))
  if(s==3 & timeframe=="stable")  par(mar=c(2,2,2,1))
  if(s==2 & timeframe=="monitor") par(mar=c(3,4,1,0))
  if(s==3 & timeframe=="monitor") par(mar=c(3,2,1,1))
  
  fullDT <- addMissingDays(subDT, mainAdd)
  
  rollM <- rollMeanDays(fullDT, n=X)
  rollM[, std := valDay / mean(fullDT$val, na.rm=TRUE)]
  
  plotTab <- rollM[(X+1):nrow(rollM)][, doy := doy-X]
  plot(plotTab$doy, plotTab$std, xlab="", ylab=ylab, xaxt=xaxt, yaxt=yaxt,
       main=main, xlim=c(1,365), ylim=c(0.45,2))
  abline(h=1, lty=2, col="red")
  
  if(timeframe=="monitor") mtext("Day of year", side=1, line=2, cex=0.8)
  if(timeframe=="stable") axis(1, at=seq(0, 300, 100), labels = FALSE)
  
  if(s==3) axis(2, at=seq(0.5, 2.5, 0.5), labels=FALSE)
  if(s==3 & timeframe=="monitor"){
    desc <- c("Dry hot", "Monsoon", "Dry cool")
    text(x=c(98, 225, 343), y=rep(2.5, 3), labels=desc)
  }
  
  # add in seasons
  ## 32 = start of hot dry season
  ## 152 = start of wet season
  ## 274 = start of cool dry season
  
  abline(v=60, col="orange")  # 2019-03-01
  abline(v=135, col="blue")   # 2019-05-15
  abline(v=305, col="purple") # 2019-11-01
  
  if(X==60){
    legend("topleft", legend=c("Start of hot dry season", "Start of wet season",
                               "Start of cool dry season"), lty=1, 
           col=c("orange", "blue", "purple"), bty="n")
  }
  
  if(returnData) return(plotTab)
}
boxNDVI <- function(subDT, mainAdd, s, timeframe){
  mainAdd <- ifelse(s==3, "SAR", "Optical")
    
  if(s==3){
    if(timeframe=="monitor") main <- "" else main <- mainAdd
    
    s1 <- subDT[grepl("sentinel1", sensor)]
    fullDT <- addMissingDays(s1, mainAdd)
    fullDT[, weekN := week(as.Date("2018-12-31")+doy)
    ][, weekN := ifelse(weekN > 52, 52, weekN)]
    boxplot(val ~ sensor + weekN, data=fullDT, ylab="", xlab="",
            main=main, outline=FALSE, col=c("grey", "violet"), xaxt="n",
            ylim=c(-20, -5))
    axis(1, at=seq(1.5, 103.5, by=2), labels=seq(1:52))
    mtext("Backscatter", side=2, line=2.5, cex=0.8)
    
    if(timeframe=="monitor"){
      legend("topright", legend=c("ASC", "DES"), col=c("grey", "violet"), 
             bty="n", pch=16)
      mtext("Week", side=1, line=2.5, cex=0.8)
    }
    
  } else {
    if(timeframe=="monitor") main <- "" else main <- mainAdd
    fullDT <- addMissingDays(subDT, mainAdd)
    fullDT[, weekN := week(as.Date("2018-12-31")+doy)
    ][, weekN := ifelse(weekN > 52, 52, weekN)]
    boxplot(val ~ weekN, data=fullDT, ylab="", xlab="",
            main=main, outline=FALSE, ylim=c(0,1))
    mtext("NDVI", side=2, line=2.5, cex=0.8)
    
    if(timeframe=="monitor") mtext("Week", side=1, line=2.5, cex=0.8)
  }
}
applySensor <- function(s, dt, windowVec=NULL, analyzeRes=FALSE, 
                        analyzeNDVI=FALSE, timeframe, returnData){
  if(s==2){
    subDT <- dt[!grepl("sentinel1", sensor), ]
  } else if(s==3){
    subDT <- dt
  }
  
  if(s==3){
    mainAdd <- paste0("Multi-source mixed")
  } else {
    mainAdd <- paste0("Multi-source optical")
  }
  
  if(analyzeNDVI & timeframe=="stable")  par(mar=c(2,4,4,1))
  if(analyzeNDVI & timeframe=="monitor") par(mar=c(4,4,2,1))
  if(analyzeRes  & timeframe=="stable")  par(mar=c(2,4,4,1))
  if(analyzeRes  & timeframe=="monitor") par(mar=c(4,4,2,1))
  
  if(analyzeRes){
    out <- sapply(windowVec, calcPlotRollMean, s, subDT, mainAdd, returnData,
                  timeframe)
    if(returnData) return(out)
  }
  if(analyzeNDVI) boxNDVI(subDT, mainAdd, s, timeframe)
}
inflationFactor <- function(dat, variance, resDates, ndvi, sensVec=NULL, 
                            windowVec=NULL, timeframe=NULL, returnData=FALSE, 
                            dates){
  # define initial table
  sensNames <- names(dat[[1]]$ptData)
  
  ## the rationale here is that when we're looking over the stable period, we
  ## can include both the disturbed and undisturbed points, because their stable
  ## periods will be representative of a similar trend.
  ## for monitoring period, though, we want the background variance / trend that
  ## is not already exacerbated by disturbances; thus we only use the 
  ## undisturbed points.
  if(variance | timeframe=="monitor"){
    dataList <- dat[grepl("a", names(dat))]
  } else {
    dataList <- dat
  }
  
  dt <- rbindlist(lapply(names(dataList), calcVarWrap, sensNames, dataList, 
                         variance, resDates, ndvi, timeframe, dates))
  
  # Compare the pre-2019 and 2019 residual variance
  ## Note that this is already combined across the sensors because we are
  ## evaluating the factor difference of the aggregated, standardized, resids
  ## that have gone through the second spike filter.
  if(variance){
    dt <- dt[!is.na(stdRes)]
    
    # average factor difference of 2019 residuals compared to pre-2019 residuals
    ## e.g. an avg of 2 = 2019 residual variance is 2x as large as pre-2019
    stable <- mean(abs(dt[period=="stable", stdRes]))
    monitor <- mean(abs(dt[period=="monitor", stdRes]))
    comb <- round(monitor / stable, 3)
    
    # visualize the overall results for variance
    # xlab <- "Factor difference monitor:stable"
    # hist(dt$stdRes, breaks=50, xlab=xlab, main="Factor diff")
    
    return(print(c(paste0("Mean combined inflation factor = ", comb))))
  }
  
  # Look at seasonality of raw residuals by day of year
  if(resDates){
    dt[, val := abs(val)]
    out <- sapply(sensVec, applySensor, dt, windowVec=windowVec, 
                  analyzeRes=resDates, analyzeNDVI=ndvi, timeframe=timeframe, 
                  returnData)
    if(returnData) return(out)
  }
  
  # Look at NDVI seasonality by week
  if(ndvi){
    sapply(sensVec, applySensor, dt, analyzeRes=resDates, analyzeNDVI=ndvi, 
           timeframe=timeframe)
  }
}
