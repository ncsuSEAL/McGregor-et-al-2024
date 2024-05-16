##########################################################
## Purpose: Functions for k-fold cross validation and accuracy plotting
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: June 2022
##########################################################
## -----------------------------------------------------------------------------
## eqRecall = equation for recall
## eqPrecision = equation for precision
## calcF1 = calculate F1 score
## f1Vec = make vector of F1 scores (and PR) for all the points
## f1VecBySize = get F1s by size partition
## f1ByLambda = wrapper (get F1 separately for each lambda)
## calcMetrics = wrapper for getting F1s and median lags for a single threshold
## metricsBySensor = wrapper to calcMetrics for each sensor
## processMetrics = wrapper for everything else
## -----------------------------------------------------------------------------
eqRecall <- function(truePos, falseNeg){
  return(truePos / (truePos + falseNeg))
}
eqPrecision <- function(truePos, falsePos){
  return(truePos/(truePos + falsePos))
}
calcF1 <- function(beta, truePos, falseNeg, falsePos){
  valRecall <- eqRecall(truePos, falseNeg)
  valPrecision <- eqPrecision(truePos, falsePos)
  
  valF1 <- 
    ((1 + beta^2) * valRecall * valPrecision) / 
    ((beta^2*valRecall) + valPrecision)
  return(c(valF1, valRecall, valPrecision))
}
f1Vec <- function(Y, tpDT, fpDT, colN){
  
  tpVec <- tpDT[, get(colN[Y])]
  tp <- sum(tpVec==1, na.rm=TRUE)
  fn <- sum(tpVec==0, na.rm=TRUE)
  
  fpVec <- fpDT[, get(colN[Y])]
  fp <- sum(fpVec==1, na.rm=TRUE)
  tn <- sum(fpVec==0, na.rm=TRUE)
  
  oa <- (tn + tp) / (length(fpVec) + length(tpVec))
  
  # # The fp number, because we're calculating it over all the ts, is
  # ## dependent on the total number of ts considered for the TP and FN
  # ## e.g. which actually had observations
  # pointsNoObs <- which(is.na(tpVec))
  # if(length(pointsNoObs) == 0){
  #   fpTab <- subDT
  # } else {
  #   fpTab <- subDT[-pointsNoObs, ]
  # }
  # 
  # fpRate <- (sum(fpTab$nStableObs) - sum(fpTab$tn)) / 
  #   sum(fpTab$nStableObs)
  # fp <- round(fpRate * nrow(fpTab))
  
  val <- calcF1(beta=1, truePos=tp, falseNeg=fn, falsePos=fp)
  return(c(val, oa))
  # return(c(tp, fn, fp))
}
f1VecBySize <- function(Y, tpDT, fpDT, colN){
  tpVec <- data.table(vec=tpDT[, get(colN[Y])], size=tpDT$distSize, 
                      point=tpDT$point)
  
  valsF1 <- sapply(sort(unique(tpVec$size)), function(sz){
    tpTab <- tpVec[size==sz, ]
    
    fn <- tpTab[vec==0, .N]
    tp <- tpTab[vec==1, .N]
    
    fpVec <- fpDT[fpDT$point %in% paste0(tpTab$point, "a"), ][, get(colN[Y])]
    fp <- sum(fpVec==1, na.rm=TRUE)
    tn <- sum(fpVec==0, na.rm=TRUE)
    
    oa <- (tn + tp) / (length(fpVec) + length(tpVec))
    
    val <- calcF1(beta=1, truePos=tp, falseNeg=fn, falsePos=fp)
    return(c(val, oa))
  })
  valDT <- as.data.table(t(valsF1))
  colnames(valDT) <- c("f1", "recall", "precision", "overallAcc")
  valDT[, `:=` (size = sort(unique(tpVec$size)), day = Y)]
  return(valDT)
}
f1ByLambda <- function(l, subDT, sens, plotWindow, size){
  subDT <- subDT[as.character(lambda)==as.character(l), ]
  
  tpDT <- subDT[!grepl("a", point), ]
  fpDT <- subDT[grepl("a", point), ]
  
  colN <- paste0("stat", 1:plotWindow)
  
  if(size){
    valDT <- rbindlist(lapply(1:plotWindow, f1VecBySize, tpDT, fpDT, colN))
  } else {
    valsF1 <- sapply(1:plotWindow, f1Vec, tpDT, fpDT, colN)
    
    valDT <- as.data.table(t(valsF1))
    colnames(valDT) <- c("f1", "recall", "precision", "overallAcc")
  }
  
  valDT[, `:=` (lambda=l, sensor=sens)]
  
  return(valDT)
}
calcMetrics <- function(Y, dat, lambdaPlot, X, window){
  datSub <- dat[thresh==Y, ]
  outVals <- rbindlist(lapply(lambdaPlot, f1ByLambda, subDT=datSub, 
                              sens=X, plotWindow=window, size=FALSE))
  outVals <- cbind(data.table(threshold=Y), outVals)
  
  lagF <- datSub[, median(lagFast, na.rm=TRUE), by=lambda]
  lagS <- datSub[, median(lagSlow, na.rm=TRUE), by=lambda]
  lags <- cbind(data.table(X, Y), lagF, lagS$V1)
  colnames(lags) <- c("sensor", "threshold", "lambda", "medLagFast", 
                      "medLagSlow")
  
  return(list(f1=outVals, lags=lags))
}
metricsBySensor <- function(X, dt, lambdaPlot, window){
  dat <- dt[sensor==X]
  metrics <- lapply(unique(dt$thresh), calcMetrics, dat, lambdaPlot, X, window)
  
  f1 <- rbindlist(lapply(metrics, '[[', 1))
  lags <- rbindlist(lapply(metrics, '[[', 2))
  
  return(list(f1=f1, lags=lags))
}
processMetrics <- function(path, type, fileLoc, lambdaPlot, window){
  # load dt
  file <- gsub("S", type, fileLoc)
  dt <- fread(paste0(path, "train5_allMetrics/", file))
  
  ## calculate f1 per day and median lags
  allMetrics <- lapply(c("L8S2", "All"), metricsBySensor, dt, lambdaPlot, 
                       window)
  
  f1All <- rbindlist(lapply(allMetrics, '[[', 1))
  lagsAll <- rbindlist(lapply(allMetrics, '[[', 2))
  pathOut <- gsub("train5", "train6", file)
  fwrite(f1All, paste0(path, gsub("allMetricsThresholds", "f1All", pathOut))) 
  fwrite(lagsAll, paste0(path, gsub("allMetricsThresholds", 
                                    "medianLagsAll", pathOut)))
  return(print(paste0("Done processing ", type)))
}
