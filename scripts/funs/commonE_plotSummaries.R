##########################################################
## Purpose: Functions to make summary (diagnostic) plots for both training and
##          landscape applications
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.3, July 2022
## Last modified: July 2022
##########################################################

## -----------------------------------------------------------------------------
## definePoint = get exact pixel number within Chatthin for a supplied or 
##                chosen coordinate
## prepPixel = process data to format for plotting
## plotPlanet = plot Planet imagery before and after the disturbance
## plotTimeseries = plot the original sensor NDVI ts with fit model
## plotAgg = plot aggregated residuals
## plotEWMA = plot ewma calculated from the pt / pixel data
## plotProbs = plot probability timeseries
## plotDiagnostic = wrapper to plot all summary plots for training data
## spotCheck = wrapper to plot all summary plots from landscape data
## compareTrainLand = compare NDVI ts for training pt and corresponding 
##                    landscape pixel in Chatthin
## -----------------------------------------------------------------------------
getPtDes <- function(nRun, matRuns, nMatRow, nCol, extent, plotCRS, 
                     shapePath, plotRun, appRegion){
  layout(matrix(1:1))
  
  ## need to plot the underlying raster and chatthin shapefile
  if(appRegion==5 & !is.null(shapePath)){
    chat <- vect(shapePath)
    chatProj <- terra::project(chat, plotCRS)
  }
  
  # plot chatthin
  valsRun <- c(rep(1:35, each=(nCol*first(matRuns$nR))), 
               rep(36, nCol*last(matRuns$nR)))
  rRun <- terra::rast(nrows=nMatRow, ncols=nCol,
                      extent=ext(extent), crs=plotCRS,
                      vals=valsRun, names="nRun")
  
  valsCell <- c(1:(nCol*nMatRow))
  rCell <- terra::rast(nrows=nMatRow, ncols=nCol, 
                       extent=ext(extent), crs=plotCRS,
                       vals=valsCell, names="nCell")
  
  r <- c(rRun, rCell)
  
  plot(r$nRun,  plg=list(cex=1.2))
  
  if(appRegion==5 & !is.null(shapePath)) lines(chatProj, col="black")
  # lines(as.lines(buffer(chatProj, width=100)), col="black")
  
  return(r)
}
definePoint <- function(coords=NULL, coordFromProbs, nRun, matRuns, nMatRow, 
                        nCol, extent, plotCRS, shapePath, appRegion){
  rastFull <- getPtDes(nRun, matRuns, nMatRow, nCol, extent, plotCRS, 
                       shapePath, plotRun=FALSE, appRegion)
  
  if(!is.null(coords)){
    if(!coordFromProbs){
      coordVec <- vect(data.frame(lon=coords[1], lat=coords[2]),
                       crs="EPSG:4326")
      coordUTM <- project(coordVec, plotCRS)
      pointDegree <- coordVec
    } else {
      coordUTM <- vect(data.frame(lon=coords[1], lat=coords[2]), 
                       crs="EPSG:32646")
      pointDegree <- project(coordUTM, "EPSG:4326")
    }
    
    plot(coordUTM, add=TRUE)
    
    nRun <- extract(rastFull, coordUTM)$nRun
    nCell <- extract(rastFull, coordUTM)$nCell
  } else {
    point <- click(r, xy=TRUE, cell=TRUE)
    pointShp <- vect(point, geom=c("x", "y"), crs=plotCRS)
    pointDegree <- project(pointShp, "EPSG: 4326")
    
    nRun <- point$nRun
    nCell <- point$nCell
    
  }
  
  # total number of cells in groups 1 - (nRun-1)
  previousCells <- matRuns[nRun, "nR"]*nCol*(nRun-1)
  
  # number of cells in this group 
  # length((previousCells+1):(endRow*nCol*nRun))
  
  # get the exact pixel number in this group corresponding to larger raster
  pixNum <- which((previousCells+1):(matRuns[nRun, "nR"]*nCol*nRun) == nCell)
  
  return(list(cell=nCell, pixNum=pixNum, nRun=nRun, pointDegree=pointDegree))
}
prepPixel <- function(allMatrices, pixNum, runType, spikeWin, spikeThresh,
                      spikeAmp, spikeTime, returnStats, returnHarMod,
                      saveStats, makeSumPlots, window, windowType, logMod,
                      endStable, inflation, staticInflation, spike2, lambda){
  # Reformat the data
  eachPixel <- lapply(1:nrow(allMatrices$landsat8sr$ts), reformat, 
                      dat=allMatrices)
  dateSens <- lapply(allMatrices, function(X) return(X$dates))
  datesRange <- as.vector(sapply(dateSens, function(X) return(range(X))))
  dates <- seq(min(datesRange), max(datesRange), by=1)
  rm(allMatrices)
  
  ## these next steps are taken directly from `processPixel` 
  ## in `appB_analyzeLandscape.R`
  # Step 4: Fit timeseries models and calculate st. residuals
  pix <- eachPixel[[pixNum]]
  
  stablePeriod <- c(dates[1], endStable)
  sensorNames <- names(eachPixel[[1]])
  stdRes <- lapply(sensorNames, combineZ, ptData=pix, 
                   dates=dateSens, runType, spikeWin, spikeThresh, 
                   spikeAmp, spikeTime, stablePeriod, returnStats, returnHarMod,
                   saveStats, makeSumPlots, inflation, staticInflation)
  names(stdRes) <- sensorNames
  
  # Step 5: Aggregate / remove duplicates and NA
  oneTS <- removeDups(stdRes=stdRes, makeSumPlots=makeSumPlots)
  
  if(spike2){
    spikesAgg <- detectOutliers(dates=oneTS$obsDates, vals=oneTS$stdRes, 
                                spikeWin, spikeThresh, spikeAmp, spikeTime)
    
    oneTS <- oneTS[spikesAgg < 4, ]
  }
  
  # Step 6: EWMA and binary categorization
  binaryCat <- applyLambda(pt=NULL, oneTS, base=NULL, lambda, dates, window, 
                           runType, windowType, stablePeriod)
  
  return(list(stablePeriod=stablePeriod, pix=pix, dateSens=dateSens, 
              stdRes=stdRes, oneTS=oneTS, 
              binaryCat=binaryCat, dates=dates))
}
plotPlanet <- function(pt, datePre, dateDist, pathPlanet, imageTime){
  tifs <- list.files(paste0(pathPlanet, "/", pt, "/files"))
  
  ##-----------------------------------------##
  #Step 1a: Plot the pre-disturbance planet imagery
  focusRast <- tifs[grepl("SR_clip", tifs)]
  
  if(imageTime=="pre"){
    rastFile <- focusRast[grepl(gsub("-", "", datePre), focusRast)]
    main <- "Pre-Disturbance: "
  } else if(imageTime=="post"){
    rastFile <- focusRast[grepl(gsub("-", "", dateDist), focusRast)]
    main <- "Post-Disturbance: "
  }
  
  if(length(rastFile) != 0){
    r <- rast(paste0(pathPlanet, "/", pt, "/files/", rastFile))
    plotRGB(r, stretch="lin", mar=1,
            main=paste0(main, strsplit(rastFile, "_")[[1]][1]))
  } else {
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10),
         axes=FALSE)
  }
}
plotTimeseries <- function(pt, dateDist, datePre, plotDates, plotOrder, 
                           dateSens, cols, pix, stdRes, runType, dates,
                           vegIndex){
  b <- sapply(1:length(plotOrder), function(i){
    if(runType=="app" & i==3){
      plot(1, type="n", xlab="", ylab="", axes=FALSE, 
           xlim=c(0,10), ylim=c(0,10))
    } else {
      if(plotOrder[i] == "sentinel1"){
        ylim <- range(pix$sentinel1ASC$ts, pix$sentinel1DES$ts, na.rm=TRUE)
        ylab <- "Backscatter"
      } else {
        ylim <- c(0,1)
        ylab <- toupper(vegIndex)
      }
      
      ## plot the original NDVI data
      if(runType=="train"){
        d <- dates
      } else {
        d <- dateSens[[i]]
      }
      
      if(grepl("sentinel1", plotOrder[i])){
        sapply(c("ASC", "DES"), function(orb){
          colActual <- ifelse(orb=="ASC", 3, 4)
          
          dat <- pix[[paste0(plotOrder[i], orb)]]
          direct <- dat$dir[!is.na(dat$dir)]
          if(grepl("ASC", orb)){
            dir <- "ASCENDING"
          }   else {
            dir <- "DESCENDING"
          }
          dat$ts <- dat$ts[dat$dir==dir]
          d <- d[dat$dir==dir]
          
          if(orb=="ASC"){
            plot(as.Date(d, origin=as.Date("1970-01-01")), dat$ts, 
                 ylim=ylim, ylab=ylab, xlab="", col=cols[colActual], 
                 main=paste0(pt, ": ", plotOrder[i]),
                 xlim=c(as.Date(min(dates), origin=as.Date("1970-01-01")),
                        as.Date(max(dates), origin=as.Date("1970-01-01"))))
          } else {
            points(as.Date(d, origin=as.Date("1970-01-01")), dat$ts, 
                   col=cols[colActual])
          }
          
          if(!all(is.na(d))){
            ## add the harmonic model fit
            harMod <- stdRes[[paste0(plotOrder[i], orb)]]$model
            lines(dates, predict(harMod, newdata=data.frame(t_mod=dates)), 
                  col=cols[colActual])
          }
        })
        legend("bottomleft", legend=c("ASC", "DES"), pch=16, 
               col=c("blue", "violet"), bty="n")
      } else {
        dat <- pix[[plotOrder[i]]]
        
        plot(as.Date(d, origin=as.Date("1970-01-01")), dat$ts, 
             ylim=ylim, ylab=ylab, xlab="", col=cols[i], 
             main=paste0(pt, ": ", plotOrder[i]))
        
        ## add the harmonic model fit
        harMod <- stdRes[[plotOrder[i]]]$model
        lines(dates, predict(harMod, newdata=data.frame(t_mod=dates)), 
              col=cols[i])
      }
      
      abline(v=datePre, lty=2, col="grey")
      abline(v=dateDist, lty=2)
      
      # par(bg="black")
      # plot(as.Date(d, origin=as.Date("1970-01-01")), dat$ts,
      #      ylim=ylim, ylab=ylab, xlab="", col=cols[i],
      #      main="Landsat-8 SR",
      #      col.lab = "white", col.main = "white", col.axis="white", fg="white")
      # abline(v=dateDist, lty=2, col="white")
      
      if(i==1){
        legend("bottomleft", legend=c("End of stable period", "Disturbance date"),
               col=c("grey", "black"), lty=2, bty="n")
        # legend("bottomleft", legend=c("Disturbance date"),
        #        col=c("white"), lty=2, bty="n", text.col="white")
      }
    }
  })
}
plotAgg <- function(pt, dateDist, datePre, xlim, aggRes){
  aggRes[, plotCol := ifelse(allSensors=="landsat8sr", "black",
                             ifelse(allSensors=="sentinel2l2a", "orange",
                                    ifelse(allSensors=="sentinel1ASC", "blue",
                                           "violet")))]
  plot(as.Date(aggRes$obsDates, origin=as.Date("1970-01-01")), aggRes$stdRes,
       xlim=xlim, col=aggRes$plotCol, ylab="Anomaly", xlab="",
       main="Combined standardized residuals")
  abline(v=datePre, lty=2, col="grey")
  abline(v=dateDist, lty=2, col="black")
  
  # aggRes[, plotCol := ifelse(allSensors=="landsat8sr", "#CCCCCC",
  #                            ifelse(allSensors=="sentinel1", "blue", 
  #                                   "#FF9933"))]
  # aggRes[, plotCol := ifelse(allSensors=="landsat8sr", cols[1],
  #                            ifelse(allSensors=="sentinel1", cols[3], 
  #                                   cols[2]))]
  # par(bg="black")
  # plot(as.Date(aggRes$obsDates, origin=as.Date("1970-01-01")), aggRes$stdRes,
  #      xlim=xlim, col=aggRes$plotCol, ylab="Anomaly", xlab="",
  #      main="Combined standardized residuals", pch=16,
  #      col.lab = "white", col.main = "white", col.axis="white", fg="white")
  # abline(v=dateDist, lty=2, col="white")
  # legend("bottomleft", legend=c("Landsat-8", "Sentinel-2", "Sentinel-1"),
  #        pch=c(16,16,16), col=c(cols), bty="n", text.col="white")
  # par(bg="white")
}
plotEWMA <- function(pt, dateDist, datePre, xlim, ewma, ylim,
                     plotLamOrder, lineCols, lambdaVec, runType){
  sapply(1:length(plotLamOrder), function(X){
    dat <- ewma[[plotLamOrder[X]]]
    if(runType=="train") dat <- dat[pointid==pt, ]
    if(X==1){
      plot(as.Date(dat$date, origin=as.Date("1970-01-01")), dat$ewma, 
           xlim=xlim, ylim=ylim, xlab="", ylab="Anomaly", type="l",
           col=lineCols[X],
           main="Weighted standardized residuals (ewma)")
    } else {
      lines(dat$date, dat$ewma, col=lineCols[X])
    }
  })
  abline(v=datePre, lty=2, col="grey")
  abline(v=dateDist, lty=2)
  legend("bottomleft", legend=paste0("Lambda=", lambdaVec), 
         col=rev(lineCols[1:length(lambdaVec)]), bty="n", lty=1)
  
  
  # plot(as.Date(dat$date, origin=as.Date("1970-01-01")), dat$ewma, 
  #      xlim=xlim, xlab="", ylab="Anomaly", type="l",
  #      col="#FF99CC",
  #      main="Weighted standardized residuals (ewma)",
  #      col.lab = "white", col.main = "white", col.axis="white", fg="white")
  # abline(v=dateDist, lty=2, col="white")
  # legend("bottomleft", legend=paste0("Lambda=", rev(lambdaVec)[X]), 
  #        col="#FF99CC", bty="n", lty=1, text.col="white")
}
plotProbs <- function(pt, dateDist, datePre, xlim, ewma, runType, logMod=NULL,
                      plotLamOrder=NULL, lineCols=NULL, dates=NULL, 
                      plotOrder=NULL, densFun=NULL, lambdaVec, bayes=FALSE,
                      allProbs=NULL, plots){
  
  if(runType=="app"){
    binaryCat <- ewma[[1]]
    # probs <- predict(logMod, data.frame(ewma=binaryCat$ewma), type="response")
    if(!is.null(densFun)){
      load(densFun[1], envir=.GlobalEnv)
      load(densFun[2], envir=.GlobalEnv)
    }
    
    if(bayes){
      probs <- applyBayes(allProbs, pt, binPt=binaryCat, runType="app")
    } else {
      probs <- sapply(binaryCat$ewma, calcProb, prior=NULL)
    }
    
    ## expand probs to fill in for all dates (i.e. if we didn't have an obs for
    ## a certain date, we fill it by bringing in the last non-NA obs)
    probsFull <- dates
    probsFull[!(dates %in% binaryCat$date)] <- NA
    probsFull[!is.na(probsFull)] <- probs
    probsFull <- nafill(probsFull, type="locf")
    
    if(plots){
      plot(as.Date(dates, origin=as.Date("1970-01-01")), probsFull, 
           xlim=xlim, ylim=c(0,1), xlab="", ylab="Probability", type="l", 
           col=lineCols, main="Probability of disturbance")
      
    }
    return(probsFull)
    
  } else if(runType=="train"){
    sapply(1:length(plotLamOrder), function(X){
      dat <- ewma[[plotLamOrder[X]]]
      dat <- dat[pointid==pt, ]
      # if(is.null(logMod)){
      #   logMod <- fitLogMod(dat, plotProbModel=FALSE, saveProbModel=FALSE,
      #                       probModelFile="", pointids, windowType)
      # }
      # probsFull <- predict(logMod, data.frame(ewma=dat$ewma), type="response")
      
      load(gsub("N", lambdaVec[X]*100, densFun[1]), envir=.GlobalEnv)
      load(gsub("N", lambdaVec[X]*100, densFun[2]), envir=.GlobalEnv)
      
      if(bayes){
        probsFull <- applyBayes(allProbs, pt, binPt=dat)
      } else {
        probsFull <- sapply(dat$ewma, calcProb, prior=NULL)
      }
      
      if(X==1){
        plot(as.Date(dat$date, origin=as.Date("1970-01-01")), probsFull,
             xlim=xlim, ylim=c(0,1), xlab="", ylab="Probability", type="l",
             col=lineCols[X], main="Probability of disturbance")
      } else {
        lines(as.Date(dat$date, origin=as.Date("1970-01-01")),
              probsFull, col=lineCols[X])
      }
    })
    abline(v=datePre, lty=2, col="grey") 
    abline(v=dateDist, lty=2, col="black")
    
    # for EDFM poster
    # par(bg="black")
    # sapply(c(1, 3), function(X){
    #   dat <- ewma[[plotLamOrder[X]]]
    #   dat <- dat[pointid==pt, ]
    #   probsFull <- predict(logMod, data.frame(ewma=dat$ewma), type="response")
    # 
    #   if(X==1){
    #     plot(as.Date(dat$date, origin=as.Date("1970-01-01")), probsFull,
    #          xlim=xlim, ylim=c(0,1), xlab="", ylab="Probability", type="l",
    #          col=col[1], main="Probability of disturbance", lwd=2,
    #          col.lab = "white", col.main = "white", col.axis="white", fg="white")
    #   } else {
    #     lines(as.Date(dat$date, origin=as.Date("1970-01-01")),
    #           probsFull, col=col[1], lwd=2)
    #   }
    # })
    # abline(v=dateDist, lty=2, col="white")
    
    # par(bg="black")
    # lineCols <- rev(c("#26E8Fc", "#CC4678FF", "#F0F921FF"))
    # sapply(1:length(plotLamOrder), function(X){
    #   dat <- ewma[[plotLamOrder[X]]]
    #   # logMod <- fitLogMod(dat, plotProbModel, saveProbModel, probModelFile)
    #   dat <- dat[pointid==pt, ]
    #   probsFull <- predict(logMod, data.frame(ewma=dat$ewma), type="response")
    # 
    #   if(X==1){
    #     plot(as.Date(dat$date, origin=as.Date("1970-01-01")), probsFull,
    #          xlim=xlim, ylim=c(0,1), xlab="", ylab="Probability", type="l",
    #          col=lineCols[X], main="Probability of disturbance", lwd=2,
    #          col.lab = "white", col.main = "white", col.axis="white", fg="white")
    #   } else {
    #     lines(as.Date(dat$date, origin=as.Date("1970-01-01")),
    #           probsFull, col=lineCols[X], lwd=2)
    #   }
    # })
    # abline(v=datePre, lty=2, col="grey")
    # abline(v=dateDist, lty=2, col="white")
    # par(bg="white")
  }
}
plotDiagnostic <- function(pt, oneTS, lambdaVec, window, runType, 
                           base, planet, vegIndex, sensKeep, logMod=NULL,
                           bayes=FALSE, allProbs=NULL){
  dateDist <- base[pointid==pt, dateDist]
  datePre <- base[pointid==pt, datePre]
  
  #Step 1. Get planet images before and after
  if(planet) layout(matrix(1:8, ncol=2)) else layout(matrix(1:6, ncol=2))
  
  ##-----------------------------------------##
  if(planet){
    #Step 1a: Plot the pre-disturbance planet imagery
    pathPlanet <- "Z:/IanMcGregor/planetImagesNew/planetscope"
    plotPlanet(pt, datePre, dateDist, pathPlanet, imageTime="pre")
    par(mar=c(2,4,4,1))
  }
  
  ##-----------------------------------------##
  #Define variables for remaining plots
  ptData <- oneTS[[pt]]$ptData
  stdRes <- oneTS[[pt]]$stdRes
  aggRes <- oneTS[[pt]]$oneTS
  plotDates <- as.Date(dates, origin=as.Date("1970-01-01"))
  plotOrder <- c("landsat8sr", "sentinel2l2a", "sentinel1")
  plotOrder <- plotOrder[grepl(gsub(", ", "|", sensKeep), plotOrder)]
  cols <- c("black", "orange", "blue", "violet")
  xlim <- c(min(plotDates), max(plotDates))
  dateSens <- lapply(stdRes, function(X){return(X$obsDates)})
  
  if(length(plotOrder==3)) nRun <- 3 else nRun <- 2
  dens <- paste0("data/dataMyanmar/train3_densFun/funStable_N_run", nRun, ".Rdata")
  densFun <- c(dens, gsub("Stable", "Dist", dens))
  
  ##-----------------------------------------##
  #Step 2. Plot original ts data along with the chosen "best" harmonic model
  plotTimeseries(pt, dateDist, datePre, plotDates, plotOrder, dateSens, cols, 
                 pix=ptData, stdRes, runType, dates, vegIndex)
  
  if(length(plotOrder) < 3){
    plot(NULL, axes=FALSE, xlim=c(0,10), ylim=c(0,10), xlab="", ylab="")
  }
  
  ##-----------------------------------------##
  if(planet){
    #Step 2b: Plot the post-disturbance tif
    plotPlanet(pt, datePre, dateDist, pathPlanet, imageTime="post")
    par(mar=c(2,4,4,1))
  }
  
  ##-----------------------------------------##
  #Step 3. Plot aggregated residuals colored by sensor source
  plotAgg(pt, dateDist, datePre, xlim, aggRes)
  
  ##-----------------------------------------##
  # Step 4. Apply lambdas and get ewma
  ewma <- lapply(lambdaVec, function(lam){
    binaryCat <- rbindlist(lapply(pointids, applyLambda, oneTS, base, 
                                  lambda=lam, dates, window, runType,
                                  windowType, staticInflation=staticInflation))
    # binaryCat[, dist := ifelse(grepl("a", pointid), 0, dist)]
    return(binaryCat)
  })
  ylim <- range(sapply(ewma, function(X){
    return(range(X[pointid==pt, ewma], na.rm=TRUE))
  }))
  
  ## we change the plot order because higher lambdas = larger data ranges
  plotLamOrder <- order(lambdaVec, decreasing=TRUE)
  lineCols <- c("#663399", "#CCCC00", "#336699")
  
  plotEWMA(pt, dateDist, datePre, xlim, ewma, ylim, plotLamOrder, lineCols, 
           lambdaVec, runType)
  
  ##-----------------------------------------##
  # Step 5. Apply lambdas and get probabilities
  plotProbs(pt, dateDist, datePre, xlim, ewma, runType, logMod,
            plotLamOrder, lineCols, dates=NULL, plotOrder, densFun, lambdaVec,
            bayes=bayes, allProbs=allProbs)
  
}
spotCheck <- function(coords=NULL, coordFromProbs, appRegion, plots, bayes){
  # Step 1: Determine which pixel to inspect
  applicationStep <- 4
  nRun <- 1
  source("dissertation/scripts/argsApp.R", local=TRUE)
  source("dissertation/scripts/args.R", local=TRUE)
  
  point <- definePoint(coords, coordFromProbs, nRun, matRuns, nMatRow, 
                       nCol, extent=appLocations[[appRegion]]$ext, plotCRS, 
                       shapePath, appRegion)
  nRun <- point$nRun
  cell <- point$cell
  
  # Step 2: Process individual pixel data
  ## now reload the args again with the correct data
  applicationStep <- 5
  source("dissertation/scripts/argsApp.R", local=TRUE)
  
  pixNum <- point$pixNum
  pixData <- prepPixel(allMatrices, pixNum, runType, spikeWin, spikeThresh,
                       spikeAmp, spikeTime, returnStats, returnHarMod,
                       saveStats, makeSumPlots, window, windowType, logMod,
                       endStable, inflation, staticInflation, spike2, lambda)
  
  # Step 3: Make the summary plots for the specific pixel
  pix <- pixData$pix
  dateSens <- pixData$dateSens
  stdRes <- pixData$stdRes
  oneTS <- pixData$oneTS
  binaryCat <- pixData$binaryCat
  dates <- pixData$dates
  
  if(sensors=="L8S2") run <- 2
  dens <- paste0(saveDir, "trainingPars/train3_densFun/funStable_", lambda*100, "_run", 
                 run, ".Rdata")
  densFun <- c(dens, gsub("Stable", "Dist", dens))
  
  if(bayes){
    load(paste0(saveDir, appLoc, "/allStaticProbs.Rdata"))
    allProbs <- allProbs[cell]/10000
  } else {
    allProbs <- NULL
  }
  
  if(plots){
    ## define necessary variables
    pt <- paste0("Pixel ", cell)
    datePre <- max(pixData$stablePeriod)
    dateDist <- max(pixData$stablePeriod) + 1
    plotDates <- as.Date(pixData$dates, origin=as.Date("1970-01-01"))
    plotOrder <- c("landsat8sr", "sentinel2l2a", "")
    cols <- c("black", "orange", "blue")
    xlim <- c(min(plotDates), max(plotDates))
    lambdaVec <- lambda
    plotLamOrder <- 1
    lineCols <- c("#663399")
    
    ## actually make the plots
    layout(matrix(1:6, ncol=2))
    plotTimeseries(pt, dateDist, datePre, plotDates, plotOrder, dateSens, 
                   cols, pix, stdRes, runType, dates, vegIndex)
    
    plotAgg(pt, dateDist, datePre, xlim, aggRes=oneTS)
    
    plotEWMA(pt, dateDist, datePre, xlim, ewma=list(binaryCat), 
             ylim=range(binaryCat$ewma, na.rm=TRUE),
             plotLamOrder, lineCols, lambdaVec, runType)
  }
  
  pr <- plotProbs(pt, dateDist, datePre, xlim, ewma=list(binaryCat), runType, 
                  logMod=NULL, plotLamOrder=NULL, lineCols, dates, plotOrder,
                  densFun, lambdaVec, bayes, allProbs, plots)
  pixData <- append(pixData, list("probs"=pr))
  
  return(list(point=point, pixData=pixData))
}
compareTrainLand <- function(chp, base){
  applicationStep <- 4
  nRun <- 1
  hpc <- FALSE
  source("scripts/argsApp.R", local=TRUE)
  source("scripts/args.R")
  
  coords <- as.numeric(unlist(base[pointid==chp, .(coordX, coordY)]))
  point <- definePoint(coords, nRun, matRuns, nMatRow, nCol, extent, 
                       plotCRS, shapePath)
  
  nRun <- point$nRun
  
  # Step 2: Process individual pixel data
  ## now reload the args again with the correct data
  applicationStep <- 5
  source("scripts/argsApp.R", local=TRUE)
  
  # Get training data for same point
  sensKeep <- "landsat8, sentinel2"
  source("scripts/argsTrain.R", local=TRUE)
  sensorNames <- names(bandList)
  ptData <- lapply(sensorNames, getIndex, dataPath=dataPath, 
                   sensorFold=sensorFold, point=chp, 
                   sensorFiles=sensorFiles, bandList=bandList, 
                   dates=dates)
  names(ptData) <- sensorNames
  
  # Plot comparisons
  
  ptDataLand <- lapply(allMatrices, function(X){
    m <- X$ts
    return(list(ts=m[point$pixNum, ]))
  })
  
  allDates <- as.Date(dates, origin=as.Date("1970-01-01"))
  
  layout(matrix(1:4, ncol=2))
  sapply(names(ptData), function(X){
    plotDates <- as.Date(allMatrices[[X]]$dates, origin=as.Date("1970-01-01"))
    plotDates <- unique(plotDates)
    
    train <- ptData[[X]]$ts
    allDates <- allDates[!is.na(train)]
    train <- train[!is.na(train)]
    
    land <- ptDataLand[[X]]$ts
    plotDates <- plotDates[!is.na(land)]
    land <- land[!is.na(land)]
    
    plot(allDates, train, xlab="", ylab="ndvi", 
         xlim=range(dates), ylim=c(0,1),
         main=paste0(chp, ": ", X, " training"))
    plot(plotDates, land, 
         xlim=range(dates), xlab="", ylim=c(0,1), ylab="ndvi", 
         main=paste0(chp, ": ", X, " landscape"))
    
    # Visualize date overlapping
    # layout(matrix(1:1))
    # full <- unique(sort(c(plotDates, allDates)))
    # fullTrain <- ifelse(full %in% allDates, 1, 0)
    # fullLand <- ifelse(full %in% plotDates, 1, 0)
    # plot(1:length(full), fullTrain, pch=16, xlab="Observation number")
    # points(1:length(full), fullLand, col="red", pch=16)
    # legend("right", legend=c("Training", "Landscape"), 
    #        col=c("black", "red"), pch=16, bty="n")
  })
}
