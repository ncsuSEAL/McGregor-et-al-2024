##########################################################
## Purpose: Functions for k-fold cross validation and accuracy plotting
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: June 2022
##########################################################
# K-fold cross validation functions
runFolds <- function(i, folds, lambda, window, pointids, runType,
                     windowType, consecObs, probType, densFun){
  print(paste0("Running fold ", i))
  
  # Calculate ewma and do binary categorization per point. It's ok to do this
  ## for both disturbed and undisturbed pts because we will use the associated
  ## undisturbed data below to help define the logistic model.
  binaryCat <- rbindlist(lapply(pointids, applyLambda, oneTS, base, lambda, 
                                         dates, window, runType, windowType,
                                         staticInflation=staticInflation))
  
  ## We're going to build the density distributions using only the training set,
  ### and then calculate the accuracy metrics for the testing set
  if(i < length(folds)+1){
    ## Identify training and testing points
    idsTest <- as.vector(unlist(folds[i]))
    idsTrain <- as.vector(unlist(folds[-i]))
    
    binaryCatTrain <- binaryCat[pointid %in% 
                                  c(idsTrain, paste0(idsTrain, "a")), ]
    binaryCatTest <- binaryCat[pointid %in% c(idsTest, paste0(idsTest, "a")), ]
  } else {
    binaryCatTrain <- binaryCat
    binaryCatTest <- binaryCat
  }
  
  funs <- createDensFuns(dat=binaryCatTrain, stable=TRUE, monitor=TRUE, 
                         saveFile=FALSE, crossVal=TRUE)
  list2env(funs, envir = .GlobalEnv)
  
  ## now generate probs for testing set
  fullProbs <- rbindlist(lapply(pointids, generateProbs, base, lambda, bayes, 
                                probType, binaryCat=binaryCatTest, logMod=NA,
                                densFun=NULL, allProbs))

  pointidsTest <- unique(binaryCatTest$pointid)
  fullProbs[, probType := NULL]
  metrics <- rbindlist(lapply(pointidsTest, retrieveMetrics, dt=fullProbs,
                              consecObs, probThresh, window, windowType, bayes,
                              lambda, dataPath))
  
  type <- ifelse(i < length(folds)+1, "subset", "full data")
  metrics[, `:=` (iter = i, type = type)]
  
  return(metrics)
}
outputF1 <- function(X, comp, colN, window){
  subDT <- comp[iter==X, ]
  
  tpDT <- subDT[!grepl("a", point), ]
  fpDT <- subDT[grepl("a", point), ]
  
  valsF1 <- sapply(1:window, f1Vec, tpDT, fpDT, colN)
  valDT <- as.data.table(t(valsF1))
  colnames(valDT) <- c("f1", "recall", "precision", "overallAcc")
  
  valDT[, iter := X][, type := ifelse(iter < 11, "subset", "full data")]
  return(valDT)
}

# Analysis functions
f1BySensor <- function(sens, dt, lambdaPlot, plotWindow, size){
  subDT <- dt[sensor==sens, ]
  outVals <- rbindlist(lapply(lambdaPlot, f1ByLambda, subDT, sens, plotWindow,
                              size))
  return(outVals)
}
formatMeta <- function(l, subDT, obsWindow, metric){
  # browser()
  subDT <- subDT[as.character(lambda)==as.character(l), ]
  
  if(metric=="fpr"){
    useDT <- subDT[grepl("a", point), ]
    bin <- 1
  } else {
    useDT <- subDT[!grepl("a", point), ]
    bin <- 0
  }
  
  metricRate <- sapply(1:obsWindow, function(column){
    metricVec <- useDT[, get(paste0("stat", column))]
    return(c(sum(metricVec==bin, na.rm=TRUE), 
             length(metricVec[!is.na(metricVec)])))
    # return(sum(fpVec==1, na.rm=TRUE) / length(fpVec[!is.na(fpVec)]))
  })
  totals <- rowSums(metricRate)
  out <- round(totals[1] / totals[2], 4)
  outDT <- data.table(out, lambda=l)
  colnames(outDT)[1] <- metric
  
  return(outDT)
  
  # fpr <- median(fpr, na.rm=TRUE)
  # numNA <- nrow(tpDT[is.na(lagFast), ])
  # lagFast <- median(tpDT$lagFast, na.rm=TRUE)
  # lagSlow <- median(tpDT$lagSlow, na.rm=TRUE)
  # lagNeg <-  median(tpDT$lagNeg, na.rm=TRUE)
  # return(data.table(lambda=l, fpr=fpr, numNA=numNA, 
  #                   medLagFast=lagFast, medLagSlow=lagSlow, medLagNeg=lagNeg))
}
metaBySensor <- function(sens, dt, lambdaPlot, obsTime, metric){
  subDT <- dt[sensor==sens, ]
  win <- obsTime[sensor==sens, wind]
  metaTab <- rbindlist(lapply(lambdaPlot, formatMeta, subDT, obsWindow=win,
                              metric=metric))
  
  metaTab[, sensor := sens]
  return(metaTab)
}

# Plot sub-functions
plotFPNR <- function(factorOrder, dt, lambdaPlot, obsTime){
  ## FPR and FNR data for plotting
  plotFPR <- rbindlist(lapply(factorOrder, metaBySensor, dt, lambdaPlot, 
                              obsTime=obsTime, metric="fpr"))
  plotFNR <- rbindlist(lapply(factorOrder, metaBySensor, dt, lambdaPlot, 
                              obsTime=obsTime, metric="fnr"))
  
  # boxplot(fpr ~ lambda + sensor, data=plotTab, 
  #         names=rep(as.character(unique(plotTab$lambda)), 2),
  #         xlab="Lambdas spanning 0.05-1",
  #         ylab="False Positive Rate",
  #         main="FPR", col=rep(c("purple", "red"),each=20))
  # legend("topleft", legend=c("With S1", "Without S1"), pch=16, 
  #        col=c("red", "purple"), bty="n")

  cols <- c("#383aaa", "#CC33FF")
  
  ## plot of FPR
  par(mar=c(1,4,2,0))
  ylim <- c(0, 0.35) # 0.3 for ch1, 0.35 for ch2
  
  plot(plotFPR$lambda, plotFPR$fpr, col=rep(cols, each=20),
       pch=c(rep(16, 20), rep(8, 20)), xlab="", ylab="Rate", xaxt="n", yaxt="n",
       ylim=ylim 
       )
  text(0.05, 0.34, labels="a", cex=1.5)
  text(0.495, 0.34, labels="FP", cex=1.5)
  axis(1, at=seq(0.2, 1, by=0.2), labels=FALSE)
  axis(2, at=seq(ylim[1], ylim[2], by=0.05), 
       labels=seq(ylim[1], ylim[2], by=0.05))
  legend("topright", legend=c("MM", "MO"), pch=c(8, 16), 
         col=c(cols[2], cols[1]), bty="n")
  
  ## plot of FNR
  par(mar=c(1,2,2,2))
  plot(plotFNR$lambda, plotFNR$fnr, col=rep(cols, each=20),
       pch=c(rep(16, 20), rep(8, 20)), xlab="", xaxt="n", yaxt="n",
       ylim=ylim
       # ylim=c(0.1,0.25)
       )
  text(0.05, 0.34, labels="b", cex=1.5)
  text(0.495, 0.34, labels="FN", cex=1.5)
  axis(1, at=seq(0.2, 1, by=0.2), labels=FALSE)
  axis(2, at=seq(ylim[1], ylim[2], by=0.05), labels=FALSE)
  
  
  return(list(plotFPR=plotFPR, plotFNR=plotFNR))
}
plotLagSensor <- function(X, medAll, quantLow, quantHigh, windowType, bayes){
  if(X=="L8S2"){
    # col <- "#383aaa"
      # main <- "Without S1"
  } else {
    # col <- "#CC33FF"
      # main <- "With S1"
  }
  col <- "black"
  
  medPlot <- medAll[sensor==X, ]
  
  quant25 <- quantLow[sensor==X, V1]
  quant75 <- quantHigh[sensor==X, V1]
  medPlot[, `:=` (quant=c(quant25, quant75),
                  type=c(rep("25%", 20), rep("75%", 20)))]
  setnames(medPlot, old="V1", new="med")
  
  if(windowType=="obs"){
    ylim <- c(0,5) 
    ylabAdd <- "observations"
  } else {
    # if(bayes) ylim <- c(0,50) else ylim <- c(0,20)
    if(bayes) ylim <- c(0,20) else ylim <- c(0,15)
    ylabAdd <- "days"
  }
  
  if(X=="L8S2"){
    par(mar=c(3,4,0,0))
    yaxt <- "s"
    ylab <- paste0("Detection Lag (", ylabAdd, ")")
    lab1 <- "c"
    lab2 <- "MO"
  } else {
    par(mar=c(3,2,0,2))
    yaxt <- "n"
    ylab <- ""
    lab1 <- "d"
    lab2 <- "MM"
  }
  
  plot(NA, NA, xlim=range(medPlot$lambda), ylim=ylim, #45 for comparison
       xlab="", main="", yaxt=yaxt, ylab=ylab)
  points(x=medPlot$lambda, y=medPlot$med, pch=17, col=col, cex=1.5)
  points(x=medPlot$lambda, y=medPlot$quant, pch=3, col=col)
  text(0.05, ylim[2]-0.5, labels=lab1, cex=1.5)
  text(0.495, ylim[2]-0.5, labels=lab2, cex=1.5)
  title(xlab=expression(lambda), line=2, cex.lab=1.5)
  
  if(X=="All"){
    axis(2, at=seq(0, 20, by=5), labels=rep("", 5))
    legend("topright", legend=c("Median value", "25th percentile"),
           pch=c(17,3), col=col, bty="n")
    # pos <- legend("topright", legend=c("Median value", "Top and \nbottom 25%"), 
    #        bty="n")
    # points(x=rep(pos$text$x, 2) - c(0.07,0.03), 
    #        y=rep(pos$text$y, each=2), 
    #        pch=c(17, 15, 3, 8), col=rep(c("#383aaa", "#CC33FF"), 2))
  }
  
  for(i in unique(medPlot$lambda)){
    subDT <- medPlot[lambda==i, ]
    lines(x=subDT$lambda, y=subDT$med, lwd=2)
  }
  
  return(medPlot)
}
plotMedianLag <- function(dt, factorOrder, windowType, bayes){
  ## Detection lag data for plotting
  tpDT <- dt[!grepl("a", point)]
  meltDT <- melt(tpDT, id.vars=c("point", "lambda", "sensor", "distSize"), 
                 measure.vars=c("lagNeg", "lagFast", "lagSlow"))
  meltDT[, col := ifelse(sensor=="All", "red", "purple")
  ][, sensor := factor(sensor, levels=factorOrder)]
  
  lags <- meltDT[variable != "lagNeg", 
  ][, variable := factor(variable, 
                         levels=c("lagFast", "lagSlow"))]
  
  ## first let's do plots of median values plus top and bottom 25th
  ## percentile
  medAll <- lags[, median(value, na.rm=TRUE), by=.(lambda, variable, sensor)]
  quantLow <- lags[, quantile(value, 0.25, na.rm=TRUE), 
                   by=.(lambda, sensor)]
  quantHigh <- lags[, quantile(value, 0.75, na.rm=TRUE), 
                    by=.(lambda, sensor)]
  
  medPlot <- lapply(c("L8S2", "All"), plotLagSensor, medAll, quantLow, 
                    quantHigh, windowType, bayes)
  medPlot <- rbindlist(medPlot)
  return(list(medPlot=medPlot, lags=lags))
}
plotAllLags <- function(bayes, windowType, lags){
  # if(bayes){
  #   ylim <- c(0,45)
  # } else {
  #   if(windowType=="obs") ylim <- c(0,10) else ylim <- c(0,40)
  # }
  ylim <- c(0,45)
  
  ticks <- c(seq(1, 40, by=1))
  ticksLab <- c(seq(1.5, 40.5, by=2))
  labs <- rep(c(unique(lags$lambda)), 2)
  
  sapply(c("L8S2", "All"), function(X){
    lagsSub <- lags[sensor==X]
    if(X=="L8S2"){
      type <- "Multi-source optical" 
      ylab <- ifelse(windowType=="obs", "Observations", "Days")
      par(mar=c(5,4,4,1))
    } else {
      type <- "Multi-source mixed"
      ylab=""
      par(mar=c(5,2,4,2))
    }
    
    boxplot(value ~ variable + lambda, data=lagsSub, outline=FALSE,
            col=rep(c("gold", "lightblue"), 40), xaxt="n", ylab=ylab,
            xlab="Lambda", main=type, ylim=ylim)
    axis(side=1, at=ticks, labels=FALSE)
    text(ticksLab, par("usr")[3]-2.5, srt = 45, xpd = TRUE, labels=labs, 
         cex=0.7)
    legend("topright", legend=c("Fastest", "Slowest"),
           col=c("gold", "lightblue"), pch=16, bty="n")
    return("yay")
  })
}
plotDistSize <- function(factorOrder, dt, lambdaPlot, plotWindow){
  par(mar=c(5,4,4,2) + 0.1)
  
  gSize <- rbindlist(lapply(factorOrder, f1BySensor, dt, lambdaPlot, 
                            plotWindow, size=TRUE))
  gSize[, `:=` (day = as.character(day))]
  gSizePlot <- ggplot(gSize) +
    aes(x = f1, fill = size) +
    geom_histogram(bins = 30L) +
    ylim(0, 550) +
    scale_fill_hue(direction = 1) +
    theme_minimal() +
    ggtitle(paste0("F1 score by size class aggregated over ",
                   "lambdas \nand post-disturbance window")) +
    facet_wrap(~sensor, labeller= 
                 labeller(sensor=c("L8S2" = "Multi-source optical",
                                   "All" = "Multi-source mixed")))
  print(gSizePlot)
  
  return(gSize)
}
plotF1PR <- function(dt, factorOrder, lambdaPlot, plotWindow, splitSize,
                     windowType, plotF1, plotPR, saveFig, figPars, plotType,
                     addedTitle){
  # Calculate F1 values per lag.
  g <- rbindlist(lapply(factorOrder, f1BySensor, dt, lambdaPlot, plotWindow,
                        size=splitSize))
  g[, `:=` (sensor = factor(sensor, levels=factorOrder))]
  
  # subset to plotWindow, then make sensor names factors
  plotDT <- rbind(
    g[sensor=="L8S2", .SD[1:plotWindow], lambda],
    g[sensor=="All", .SD[1:plotWindow], lambda]
  )
  
  plotDT[, sensorDes := ifelse(sensor=="L8S2", "Multi-source optical", 
                             "Multi-source mixed")
  ][, sensorDes := factor(sensorDes, 
                          levels=c("Multi-source optical", 
                                   "Multi-source mixed"))]
  
  # Plotting
  colLine <- viridis::viridis(length(lambdaPlot))
  # col <- c("#26E8Fc", "#CC4678FF", "#F0F921FF")
  xlab <- ifelse(windowType=="days", "Daily lag", "Observation lag")
  
  if(plotF1){
    if(splitSize){
      usedPts <- length(unique(dt[!grepl("a", point), point]))
      
      f1 <- ggplot(plotDT) +
        aes(x=c(rep(1:plotWindow, length(lambdaPlot)), 
                rep(1:plotWindow, length(lambdaPlot))), 
            y=f1) +
        geom_line(aes(color=as.character(lambda))) +
        scale_color_manual(values=colLine, name="Lambda") +
        xlab(xlab) +
        ylim(0,1) +
        ylab("F1 score") +
        facet_wrap(~sensorDes) +
        ggtitle(paste0(usedPts, "/315 points, ", addedTitle)) +
        theme_minimal()
    } else {
      # colLine <- c(rep("grey50", 7), "red", rep("grey50", 12))
      f1 <- ggplot(plotDT) +
        aes(x=c(rep(1:plotWindow, length(lambdaPlot)), 
                rep(1:plotWindow, length(lambdaPlot))), 
            y=f1) +
        geom_line(aes(color=as.character(lambda))) +
        scale_color_manual(values=colLine, name="Lambda") +
        # scale_x_continuous(breaks = function(x){
        #       unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))}) +
        xlab(xlab) +
        ylab("F1 score") +
        ylim(0,1) +
        # xlim(0,20) +
        facet_wrap(~sensorDes) +
        # facet_wrap(~sensor, labeller=facetLabs) + 
        # geom_vline(xintercept=10, linetype="dashed") +
        theme_minimal()
      # dark_theme_minimal()
      # ggtitle(paste0("F1 score", addedTitle))
    } 
    
    if(saveFig){
      pars <- figPars[img=="f1", ]
      tiff(paste0(imgFolder, "figXX_f1Scores_", windowType, "_", plotType, 
                 ".tiff"), 
          width=pars$width, height=pars$height, units=pars$units, 
          res=pars$res, compression="lzw")
      print(f1)
      dev.off()
    } else {
      print(f1)
    }
  }
  
  if(plotPR){
    rangeVec <- c(min(lambdaPlot), median(lambdaPlot), max(lambdaPlot))
    plotDT[, dayN := paste0("Day ", rep(1:60, 40))
    ][, sensor := ifelse(sensor=="L8S2", "Without S1", "With S1")]
    
    pr <- ggplot(plotDT[dayN %in% paste0("Day ", c(1, 20, 40, 60))]) +
      aes(x=recall, y=precision) +
      #geom_path(linetype=rep(1:30, 9)) +
      geom_point(aes(color=as.character(lambda))) +
      scale_color_manual(values=c(colLine), name="Lambda") +
      xlab("Recall") +
      ylab("Precision") +
      ylim(0.9,1) +
      facet_grid(sensor~dayN) +
      theme_grey() +
      ggtitle("PR curve")
    
    if(saveFig){
      pars <- figPars[img=="pr", ]
      tiff(paste0(imgFolder, "figXX_prCurve_", windowType, "_", plotType, 
                 ".tiff"), 
          width=pars$width, height=pars$height, units=pars$units, 
          res=pars$res, compression="lzw")
      print(pr)
      dev.off()
    } else {
      print(pr)
    }
  }
  
  return(g)
}

# Main plot functions
plotOtherData <- function(sens, plotTab){
  subDT <- plotTab[sensor==sens]
  
  par(mar=c(5,4,4,4)+0.3)
  plot(subDT$lambda, subDT$fpr, ylab="falsePosRate", 
       xlab="Lambda", type="l",
       lwd=2, main=sens, ylim=range(plotTab$fpr))
  abline(v=subDT[fpr==min(fpr), lambda][1], lty=2)
  par(new=TRUE)
  plot(subDT$lambda, subDT$numNA, axes=FALSE, bty="n", xlab="", ylab="", 
       col="red", type="l", lwd=2, ylim=range(plotTab$numNA, na.rm=TRUE))
  abline(v=subDT[numNA==min(numNA), lambda][1], col="red", lty=2)
  axis(side=4, at=max(plotTab$numNA)*c(0, 0.5, 1), col="red", 
       col.ticks="red")
  mtext("numNA", side=4, line=3, cex=0.8)
  par(new=TRUE)
  plot(subDT$lambda, subDT$medLagFast, axes=FALSE, bty="n", 
       xlab="", ylab="", 
       col="blue", type="l", ylim=c(0, max(plotTab$medLagSlow, 
                                           plotTab$medLagFast, na.rm=TRUE)))
  lines(subDT$lambda, subDT$medLagSlow, col="purple")
  axis(side=4, at=max(plotTab$medLagSlow, plotTab$medLagFast)*c(0.25, 0.75), 
       col="blue", col.ticks="blue")
}
plotAccuracy <- function(f, foldPath, randomLam, consecObs,
                         plotF1, plotPR, plotMeta, base, splitSize, splitNum,
                         plotWindow, windowType, saveFig, figPars, 
                         imgFolder, bayes, probType){
  if(bayes){
    f <- f[grepl("bayes", f)]
  } else {
    f <- f[!grepl("bayes", f)]
    # f <- f[grepl(probType, f)]
  }
  
  # dt <- rbindlist(lapply(f, function(X){
  #   p <- fread(paste0(foldPath, X))
  #   p[, sensor := strsplit(X, "_")[[1]][2]]
  #   return(p)
  # }))
  dt <- fread(paste0(foldPath, f))
  
  dt <- dt[, distSize := base$distPercL8[match(point, base$pointid)]
           ][, distSize := ifelse(distSize <= 0.9, "<=90%", "100%")]
  
  # dist <- dt[distType=="normal"]
  # sapply(paste0("stat", 1:60), function(q){
  #   subDT <- dist[lambda==0.4, ..q]
  #   numNA <- length(which(is.na(subDT[,1])))
  #   return(1-round(numNA/nrow(subDT), 2))
  # })
  
  if(splitSize){
    ptsDist <- base[, .(pointid, distPercL8)]
    setkeyv(ptsDist, "pointid")
    setkeyv(dt, "point")
    dt <- dt[ptsDist]
    
    threshold <- as.numeric(splitNum[2])
    
    if(splitNum[1]=="at most"){
      ptsFocus <- unique(dt[!grepl("a", point) & 
                              distPercL8 <= threshold, point])
      addedTitle <- paste0("<= ", 100*threshold, "% of Landsat-8 pixel")
    } else if(splitNum[1]=="at least"){
      ptsFocus <- unique(dt[!grepl("a", point) & 
                              distPercL8 >= threshold, point])
      addedTitle <- paste0(">= ", 100*threshold, "% of Landsat-8 pixel")
    }
    
    dt <- dt[point %in% c(ptsFocus, paste0(ptsFocus, "a"))]
    
    plotTypeAdd <- ifelse(splitNum[1]=="at least", "AL", "AM")
    plotType <- paste0("split", plotTypeAdd, splitNum[2])
  } else {
    addedTitle <- ""
    plotType <- "allData"
  }
  
  factorOrder <- c("L8S2", "All")
  
  if(randomLam){
    lambdaPlot <- seq(0.1, 1, length.out=lambdaN)
  } else {
    lambdaPlot <- unique(dt$lambda)
  }
  
  if(plotF1 | plotPR){
    outTab <- plotF1PR(dt, factorOrder, lambdaPlot, plotWindow, splitSize,
                       windowType, plotF1, plotPR, saveFig, figPars, plotType,
                       addedTitle)
  } else {
    outTab <- NA
  }
  
  if(plotMeta){
    ## We change the plot window when plotting FPR and FNR, bc when using
    ## windowType="obs", if we still use 60 as we do for "days", then
    ## we have very few points that have 60 observations prior to the end of the
    ## primary timeframe.
    if(windowType=="days") win <- plotWindow else win <- 15
    obsTime <- data.table(sensor=c("L8S2", "All"), wind=rep(win, 2))
    
    # ----------------------------------------------#
    # Plot fpr, fnr, and median lags in one figure
    if(saveFig){
      pars <- figPars[img=="accLag", ]
      tiff(paste0(imgFolder, "figXX_accLag_", windowType, "_", plotType, 
                  ".tiff"), width=pars$width, height=pars$height, 
           units=pars$units, res=pars$res, compression="lzw")
    }
    layout(matrix(1:4, nrow=2, byrow=TRUE))
    fpnr <- plotFPNR(factorOrder, dt, lambdaPlot, obsTime)
    lagData <- plotMedianLag(dt, factorOrder, windowType, bayes)
    par(mar=c(5,4,4,2)+0.1)
    
    # ----------------------------------------------#
    ## now we do boxplots of all the lag times
    if(saveFig){
      dev.off()
      
      pars <- figPars[img=="boxLag", ]
      tiff(paste0(imgFolder, "figXX_boxLag_", windowType, "_", plotType, 
                  ".tiff"),  width=pars$width, height=pars$height, 
           units=pars$units, res=pars$res, compression="lzw")
    }
    layout(matrix(1:2, nrow=1))
    plotAllLags(bayes, windowType, lags=lagData$lags)
    
    if(saveFig) dev.off()

    # ----------------------------------------------#
    ## accuracy by size class
    ### percDistL8=1 means >= 100% Landsat pixel
    if(!splitSize){
      pars <- figPars[img=="f1Size", ]
      if(saveFig){
        tiff(paste0(imgFolder, "figXX_f1Size", "_", plotType, ".tiff"), 
            width=pars$width, height=pars$height, units=pars$units, 
            res=pars$res, compression="lzw")
      }
      gSize <- plotDistSize(factorOrder, dt, lambdaPlot, plotWindow)
      if(saveFig) dev.off()
    } else {
      gSize <- NA
    }
    
    # ----------------------------------------------#
    par(mar=c(5,4,4,2) + 0.1)
    layout(matrix(1:1))
    
    fpr <- fpnr$plotFPR
    fnr <- fpnr$plotFNR
    lags <- lagData$medPlot
    lagsFull <- lagData$lags
  } else {
    fpr <- NA
    fnr <- NA
    lags <- NA
    lagsFull <- NA
    gSize <- NA
  }
  return(list(tabF1PR=outTab, fpr=fpr, fnr=fnr, lags=lags, lagsFull=lagsFull, 
              size=gSize))
  
  # Temp code to plot different lambdas from L8S2 and All on top of each other to
  ## have an easier comparison
  # layout(matrix(1:4, ncol=2, byrow=TRUE))
  # 
  # sapply(c(0.1, 0.2, unique(g$lambda)[3], 0.4), function(X){
  #   sub <- g[lambda==X, ]
  #   plot(1:30, sub[sensor=="L8S2", f1], type="l", main=paste0("lambda=", X),
  #        ylab="f1", ylim=c(0.55,0.96))
  #   lines(1:30, sub[sensor=="All", f1], col="blue")
  #   if(X==0.1) legend("bottomright", legend=c("L8S2", "All"), 
  #                     col=c("black", "blue"), lty=1, bty="n")
  #   abline(v=5, lty=2)
  # })
  
  ## old plot of just the failed detections
  # numNA <- tpDT[is.na(lagFast), .N, by=.(lambda, sensor)]
  # plot(numNA$lambda, numNA$N, pch=16, col=c(rep("red",20), rep("purple", 20)),
  #      xlab="Lambda", ylab="Number of failed detections")
  # legend("bottomleft", legend=c("With S1", "Without S1"), 
  #      col=c("red", "purple"), pch=16, bty="n")
}

# Distribution of lag times
plotLag <- function(fDates, fMetrics, l, sensorCombo, foldPath, window, 
                    windowType, obs, plotBox, plotECDF){
  
  dtList <- lapply(c("metrics", "dates"), function(Y){
    if(Y=="metrics"){
      p <- fread(fMetrics)
      p <- p[sensor == sensorCombo, ]
      p <- p[!grepl("a", point), ]
    } else if(Y=="dates"){
      files <- fDates[grepl(sensorCombo, fDates)]
      p <- fread(files)
    }
    return(p)
  })
  
  dtMet <- dtList[[1]]
  dtMet <- dtMet[as.character(lambda)==as.character(l), ]
  rangeMet <- max(dtMet$lagSlow, dtMet$lagFast, na.rm=TRUE)
  
  ## change NAs to be 70 in order to have the true % be represented on the ecdf
  dtMet[, `:=` (lagFast = ifelse(is.na(lagFast), 9999, lagFast),
                lagSlow = ifelse(is.na(lagSlow), 9999, lagSlow))]
  
  dtObs <- dtList[[2]]
  subtractDist <- rbindlist(lapply(unique(dtObs$pointid), function(pt){
    sub <- dtObs[pointid==pt, ]
    dateDist <- as.numeric(base[pointid==pt, dateDist])
    sub[, diff := datesObs - dateDist]
    return(sub)
  }))
  
  subtractDist[, obsNum := rep(seq(1,window), (nrow(subtractDist)/window))]
  if(obs==1){
    ylab1 <- "% disturbances detected"
    ylab2 <- "Days"
    yaxt <- "s"
    xlim <- c(obs-1, rangeMet+1)
    par(mar=c(3,3,2,0))
  } else if(obs==3){
    ylab1 <- ""
    ylab2 <- ""
    yaxt <- "n"
    xlim <- c(obs-1, rangeMet)
    par(mar=c(3,2,2,1))
  }
  
  if(plotECDF){
    if(windowType=="days"){
      xlab <- "Days since disturbance"
      xlim <- c(0,60)
    } else {
      xlab <- "Obs since disturbance"
      xlim <- c(0,20)
    }
    
    if(sensorCombo=="L8S2"){
      nm <- "MO"
    } else {
      nm <- "MM"
    }
    
    plot(ecdf(dtMet$lagFast), main=paste0(nm, ": nObs =", obs), 
         xlim=xlim, ylim=c(0,1), verticals=TRUE, col="#383aaa", pch=NA, lwd=2,
         xlab="", yaxt=yaxt, ylab="", col.01line="transparent")
    lines(ecdf(dtMet$lagSlow), col="#66CCFF", verticals=TRUE, pch=NA,
          col.01line="transparent")
    title(xlab=xlab, line=2)
    abline(h=0.5, lty=2)
    # abline(v=10, lty=2)
    # abline(v=20, lty=2)
    
    if(obs==1){
      title(ylab=ylab1, line=2)
    } else if(obs==3){
      axis(2, at=seq(0, 1, 0.2), labels=rep("", 6))
      legend("bottomright", col=c("#383aaa", "#66CCFF"), bty="n", lwd=c(2,1), 
             cex=0.8,
             legend=c("Fastest detection", "Slowest detection"))
    }
  }
  
  if(obs==1){
    xlim <- c(0,4)
    par(mar=c(4,4,1,1))
    xaxt <- "s"
  } else if(obs==3){
    xlim <- c(0, 10)
    par(mar=c(4,1,1,1))
    xaxt <- "n"
  }
  
  if(plotBox){
    boxplot(subtractDist[obsNum <= 30, diff] ~ 
              subtractDist[obsNum <= 30, obsNum],
            yaxt=yaxt, xaxt=xaxt, xlab="Observations since disturbance",
            ylab=ylab2)
    if(obs==3){
      axis(1, at=1:60, labels=3:62)
    }
  }
  
  return(dtMet[lambda==l,])
}

# Performance by size class
plotLagSizeClass <- function(f, l, sensorCombo){
  dtList <- lapply(1:2, function(Y){
    dataType <- ifelse(Y==1, "_metrics.csv", "_obsDates.csv")
    files <- f[grepl(dataType, f)]
    
    if(dataType=="_metrics.csv"){
      out <- rbindlist(lapply(files, function(X){
        p <- fread(paste0(foldPath, X))
        p[, sensor := strsplit(X, "_")[[1]][2]]
        p <- p[sensor==sensorCombo & !grepl("a", point), ]
        return(p)
      }))
    } else {
      p <- fread(paste0(foldPath, files[grepl(sensorCombo, files)]))
      p[, sensor := as.character(sensor)][, sensor:=sensorCombo, ]
      return(p[fill==0, ])
    }
  })
  dtMet <- dtList[[1]]
  ptsDist <- base[!grepl("a", pointid), .(pointid, distPercL8)]
  setkeyv(ptsDist, "pointid")
  setkeyv(dtMet, "point")
  dtMet <- dtMet[ptsDist]
  
  dtMet[, relSize := ifelse(distPercL8==1, "Large", "Small")]
  
  metricsBySize <- function(dtMet){
    fastMed <- dtMet[, median(lagFast, na.rm=TRUE), by=.(relSize)
    ][,`:=` (type = "median", speed="fast")]
    slowMed <- dtMet[, median(lagSlow, na.rm=TRUE), by=.(relSize)
    ][, `:=` (type = "median", speed="slow")]
    
    fastMean <- dtMet[, mean(lagFast, na.rm=TRUE), by=.(relSize)
    ][, `:=` (type = "mean", speed="fast")]
    slowMean <- dtMet[, mean(lagSlow, na.rm=TRUE), by=.(relSize)
    ][, `:=` (type = "mean", speed="slow")]
    
    return(rbind(fastMed, slowMed, fastMean, slowMean))
  }
  
  out <- metricsBySize(dtMet)
  ggplot(out) +
    aes(x = relSize, y = V1, colour = type) +
    geom_point(shape = "circle", size = 1.5) +
    ylab("Observations") + 
    xlab("Relative size") + 
    scale_color_hue(direction = 1) +
    theme_minimal() +
    facet_wrap(~speed)
}

# Mean days per sensor combo observation
avgDaysPerObs <- function(datesFiles, backfill=FALSE, window, type){
  dt <- rbindlist(lapply(datesFiles, function(X){
    p <- fread(paste0(foldPath, X))
    p[, sensor := as.character(sensor)
      ][, sensor := gsub("_obsDates.csv", "", X)]
    return(p)
  }))
  
  # get mean length of time (since disturbance) for every observation
  dt <- dt[fill==ifelse(backfill, 1, 0), ]
  
  subtractDist <- rbindlist(lapply(unique(dt$pointid), function(pt){
    sub <- dt[pointid==pt, ]
    dateDist <- as.numeric(base[pointid==pt, dateDist])
    sub[, diff := datesObs - dateDist, by = .(sensor)]
    return(sub)
  }))
  
  # create output table. Note that the "by = rep(1:30, .N/30)" can be replaced
  ## with "by = (1:.N) %% 30", but will need to replace the 0s with "30".
  if(type=="mean"){
    out <- subtractDist[, .SD[, by = rep(1:window, .N/window), 
                              mean(diff, na.rm=TRUE)], by = sensor]
  }
  
  if(type=="median"){
    out <- subtractDist[, .SD[, by = rep(1:window, .N/window), 
                              median(diff, na.rm=TRUE)], by = sensor]
  }
  setnames(out, old="V1", new="time")
  
  return(out)
}
