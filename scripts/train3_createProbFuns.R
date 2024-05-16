##########################################################
## Purpose: Calculate density functions for probability estimates
## Creator: Ian McGregor, imcgreg@ncsu.edu
## Run medium:
##        - Either Mac or PC, this is run on training data
## System: R Version 4.2.2, Feb 2023
## Last modified: Mar 2023
##########################################################
groundhog.library(data.table, groundhogDate)
groundhog.library(parallel, groundhogDate)
groundhog.library(ggplot2, groundhogDate)
groundhog.library(ggpointdensity, groundhogDate)
groundhog.library(viridis, groundhogDate)
groundhog.library(RColorBrewer, groundhogDate)
groundhog.library(cowplot, groundhogDate)
groundhog.library(ggpubr, groundhogDate)
groundhog.library(moments, groundhogDate)

source("scripts/funs/trainC_byLambda.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")
source("scripts/funs/commonF_probBayes.R")
source("scripts/funs/trainG_densFuns.R")
funs <- ls()

# Likelihood (run with training data)
probHist <- function(outData, plotHist, histInterval){
  probDT <- data.table(breaks=c(seq(histInterval, 1, by=histInterval)))
  
  if(plotHist) layout(matrix(1:2))
  
  for(i in 0:1){
    if(i==0){
      sl <- outData[grepl("a", names(outData))]
    } else {
      sl <- outData[!grepl("a", names(outData))]
    }
    
    h <- hist(unlist(sl), breaks=1/histInterval, plot=FALSE)
    h$counts <- h$counts / sum(h$counts)
    
    if(plotHist){
      main <- paste0("Probability of getting a value on x-axis \ngiven", 
                     "pixel is not disturbed")
      if(i==1) main <- gsub("not ", "", main)
      
      plot(h, freq=TRUE, yaxt="n", ylim=c(0,1), ylab="Probability", 
           main=main, xlab="Algorithm output value")
      axis(2, at=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))
    }
    
    probDT[, `:=` (prob = h$counts)]
    
    if(i==0) setnames(probDT, old="prob", new="probNoDist")
    if(i==1) setnames(probDT, old="prob", new="probDist")
  }
  return(probDT)
}
testNorm <- function(X, outData, l){
  dat <- outData[[X]]
  dat <- rbindlist(dat)
  d0 <- c(dat[dist==0, ewma], dat[grepl("a", pointid) & dist==1, ewma])
  hist(d0, prob=TRUE, breaks=50, main=paste0("Lambda = ", l[X]))
  
  return(c(skewness(d0), kurtosis(d0)))
  
  # plot(d0$ewma, main="Undisturbed ewma values")
  # hist(d0$ewma, breaks=50)
  # 
  # # We remove g42p2_7a because it has terrible model fit
  # d0 <- d0[!grepl("g42p2_7a", pointid), ]
  # test <- d0[ewma <= 2 & ewma >= -2, ]
  # 
  # plot(d0$ewma, main="Large outlier removed")
  # abline(h=2, lty=2, col="white")
  # abline(h=-2, lty=2, col="white")
  # hist(test$ewma, breaks=50, main="Only data btwn 2 and -2 (98% of data)")
  # 
  # # We can't rely on KS-test here because it is highly sensitive to any 
  # ## deviation from normality given the large sample size. E.g. here, a ks-test 
  # ## gives a p-value that says our data is strongly not normal
  # 
  # # Similarly, the p-value from a t-test is sensitive to the seed used
  # ## when generating rnorm to compare against. However, most of the time in this
  # ## case, the p-value is insignificant, indicating our data has a mean that is
  # ## considered equivalent to an ideal normal distribution.
  # print("Check the previous plots as well")
  # n <- 60
  # set.seed(34)
  # p <- sapply(sample(1:5000, n), function(X){
  #   set.seed(X)
  #   normIdeal <- rnorm(nrow(d0))
  #   pv <- t.test(d0$ewma, normIdeal)
  #   return(pv$p.value)
  # })
  # 
  # layout(matrix(1:1))
  # boxplot(p, main=paste0("Distribution of p-values over ", n, " t-tests"))
}
prepPlotData <- function(fileStable, fileMonitor, g){
  ## load the density estimators
  load(fileStable, envir=.GlobalEnv)
  load(fileMonitor, envir=.GlobalEnv)
  
  ##-------------------------------------------------------##
  ## Panels A+B Prep data for histograms
  gNoDist <- g[grepl("a", pointid) & dist==0] #undisturbed points, stable period
  gDist <- g[!grepl("a", pointid) & dist==1, ] #disturbed points, dist window
  
  histTab <- rbind(gNoDist, gDist)
  
  ## Panel A (ggHist)
  ### histograms showing distr. of zEWMA values for stable and monitoring
  ### points
  histTab[, dist := as.character(dist)]
  q <- density(histTab$ewma)
  
  r <- range(histTab$ewma)
  zSeq <- seq(r[1], r[2], 0.1)
  probSeq <- sapply(zSeq, calcProb, prior=NULL)
  seqVals <- data.table(z=zSeq, prob=probSeq)
  seqVals[, y := probSeq*max(q$y)]
  
  ## Panel B (ggHistRS)
  ### histograms showing distributions of zEWMA values for stable and 
  ### monitoring points but converted to prob
  testDist <- sapply(gDist$ewma, calcProb, prior=NULL)
  testNoDist <- sapply(gNoDist$ewma, calcProb, prior=NULL)
  histTabRS <- rbind(data.table(prob=testDist, dist="1"), 
                  data.table(prob=testNoDist, dist="0"))
  
  ##-------------------------------------------------------##
  ## Panels C + D Contour and density plots
  ### use gDist as created above (disturbed points only in dist window)
  
  ### calculate simulation of posteriors. Note we bound the zEWMA to max 0 
  ### for visualization purposes because the higher lambdas do extend > 0.
  prior <- seq(0.01,1,by=0.01)
  zewma <- rev(seq(min(gDist$ewma), 0, by=0.1))
  # zewma <- rev(seq(min(gDist$ewma), max(gDist$ewma), by=0.5))
  
  ### calculate the labels for the plot using no prior
  plotLabs <- sapply(zewma, calcProb, prior=NULL)
  plotLabs <- round(plotLabs, 2)
  
  ### calculate actual posteriors for the contour using the priors
  #### We also bound this to 1 for visualization purposes (i.e. there are 
  #### multiple instances of "1", so we only want to show the first)
  zewma <- zewma[which(plotLabs < 1)]
  datDist <- gDist[ewma >= min(zewma) & ewma <= max(zewma)]
  
  posts <- sapply(prior, function(p){
    out <- sapply(zewma, calcProb, prior=p)
    return(out)
  })
  
  ### calculate breaks for the plots
  breaksN <- seq(0,1, length.out=5)
  breaks <- seq(max(zewma), min(zewma), length.out=5)
  # breaks <- zewma[breaksN]
  # plotLabs <- plotLabs[as.vector(breaksN)]
  
  ### prepare table for contour plotting
  bl <- melt(data.table(posts))
  bl[, `:=` (zewma = rep(zewma, length(prior)),
             prior = rep(prior, each=length(zewma)),
             posterior = value)]
  bl <- bl[!is.na(value), ]
  
  ##-------------------------------------------------------##
  return(list(histTab=histTab, seqVals=seqVals, histTabRS=histTabRS, bl=bl, 
              datDist=datDist, breaks=breaks, breaksN=breaksN))
}
saveDensityEst <- function(X, l, outData, nRun, filePath, stable, monitor){
  dat <- outData[[X]]
  # dat <- rbindlist(dat)
  
  fileStable <- gsub(".Rdata", 
               paste0("Stable_", l[X]*100, "_run", nRun, ".Rdata"), 
               filePath)
  fileMonitor <- gsub("Stable", "Dist", fileStable)
  
  createDensFuns(dat, stable, monitor, saveFile=TRUE, fileStable, fileMonitor)
  
  return(c(fileStable, fileMonitor))
}
defineDens <- function(n, l, outData, fileNames, nRun, allProbs, testNormal, 
                       plotHist, plotContours){
  
  g <- outData[[n]]
  setkeyv(g, "pointid")
  setkeyv(allProbs, "pointid")
  g <- g[allProbs]
  
  if(plotHist){
    gDist <- g[!grepl("a", pointid)]
    gDist <- gDist[dist==1, ewma]
    hist(gDist, prob=TRUE, breaks=50, xlab="zEWMA", 
         main=paste0("Prob Hist of zEWMA for lambda = ", l[n]))
    lines(density(gDist))
  }
  
  if(plotContours){
    plotPars <- prepPlotData(fileStable=fileNames["stable",n],
                             fileMonitor=fileNames["monitor", n], 
                             g)
    ## Distribution of priors comparing dist and undist points
    # allProbs[, dist := ifelse(grepl("a", pointid), "0", "1")]
    # ggplot(allProbs, aes(x=prior)) + 
    #   geom_density(aes(y=after_stat(density), fill=dist), alpha=.2) +
    #   geom_density(aes(y=after_stat(density), fill=dist), alpha=.2) +
    #   scale_fill_manual(values = c("blue", "red"), labels=c("D=0", "D=1"),
    #                     name="")
    
    ## Panel A
    if(n==20) yax2 <- "Posterior probability" else yax2 <- ""
    ggHist <- ggplot(plotPars$histTab, aes(x=ewma)) + 
      geom_density(aes(y=after_stat(density), fill=dist), alpha=.2) +
      geom_density(aes(y=after_stat(density), fill=dist), alpha=.2) +
      # geom_vline(xintercept = 0, linetype="dashed", color = "white") +
      scale_fill_manual(values = c("blue", "red"), labels=c("D=0", "D=1"),
                        name="") +
      geom_path(data=plotPars$seqVals, aes(x=z, y=y), col="black", 
                linewidth=1.5) +
      scale_y_continuous(sec.axis=sec_axis(~. / max(plotPars$seqVals$y), 
               name = "Posterior probability")) +
      # geom_hline(yintercept=max(plotPars$seqVals$y)*0.5, linetype="dashed") +
      ylab("Density") +
      xlab("zEWMA") +
      ggtitle(paste0("Lambda = ", l[n])) +
      theme_minimal()
    
    ## Panel B
    ggHistRS <- ggplot(plotPars$histTabRS, aes(x=prob)) + 
      geom_density(aes(y=after_stat(density), fill=dist), alpha=.2) +
      geom_density(aes(y=after_stat(density), fill=dist), alpha=.2) +
      # geom_vline(xintercept = 0, linetype="dashed", color = "white") +
      scale_fill_manual(values = c("blue", "red"), labels=c("D=0", "D=1"),
                        name="") +
      ylab("Density") +
      # scale_x_continuous(breaks=breaksN, labels=breaksN) +
      xlab("") +
      theme_minimal()
    
    ## Panel C
    con <- ggplot(plotPars$bl, aes(x=zewma, y=prior, z=value)) +
      geom_contour_filled() + 
      labs(fill="Posterior") +
      scale_fill_manual(values=c("black", rev(brewer.pal(9, "Blues"))),
                        labels=as.character(seq(0.1,1,by=0.1))) +
      # xlim(xlim[1], xlim[2]) +
      # scale_x_continuous(breaks=breaks,labels=breaksN) +
      ylab("Prior probability") +
      scale_x_reverse(breaks=plotPars$breaks,labels=plotPars$breaksN) +
      xlab("") +
      theme_minimal()
    
    ## Panel D
    dens <- ggplot(plotPars$datDist, aes(x = ewma, y = prior)) +
      stat_density2d(aes(fill = after_stat(density)^0.25), geom = "tile", 
                     contour = FALSE, n = 200) +
      labs(fill="Density") +
      ylab("Prior probability") +
      ylim(0, 1) +
      scale_fill_continuous(type="viridis", breaks = c(0.1, 1.65), 
                            labels = c("Low", "High")) +
      # scale_x_continuous(breaks=breaks,labels=breaksN) +
      xlab("RS probability of disturbance") +
      scale_x_reverse(breaks=plotPars$breaks,labels=plotPars$breaksN) +
      theme_minimal()
    
    if(n != 1){
      ggHist <- ggHist +     ylab("")
      ggHistRS <- ggHistRS + ylab("")
      con <- con +           ylab("")
      dens <- dens +         ylab("")
    }
    
    return(list(ggHist=ggHist, ggHistRS=ggHistRS, con=con, dens=dens))
  }
}
getLikelihood <- function(nRun, ewmaDays, filePath, funs, testNormal, createFuns){
  combos <- c("sentinel1", "landsat8, sentinel2", 
              "landsat8, sentinel2, sentinel1")
  fileNames <- c("S1", "L8S2", "All")
  
  numRun <- as.numeric(nRun)
  sensKeep <- combos[numRun]
  returnOnlyMod <- FALSE
  bayes <- FALSE
  if(returnOnlyMod) script <- 2 else script <- 3
  source("scripts/args.R", local=TRUE)
  source("scripts/argsTrain.R", local=TRUE)
  
  #Bring in the residuals
  fileLoadLoc <- paste0(dataPath, "/trainingPars/train1_", fileNames[numRun])
  load(paste0(fileLoadLoc, ".Rdata"))
  
  # get the training data ts probabilities for all points for just the 
  ## post-disturbance window
  if(ewmaDays){
    ewmaOnly <- TRUE
    returnProbs <- FALSE
  } else {
    returnProbs <- TRUE
    ewmaOnly <- FALSE
  }
  
  ## run the main code
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("oneTS", "base", "dates", "window", "probThresh", 
                      "runType", "windowType", "pointids", "dataPath",
                      "staticInflation", "returnProbs"), envir=environment())
  clusterExport(cl, funs, envir=.GlobalEnv)
  
  ## get observation metrics or logmod only depending on context
  ### note: we set bayes as FALSE here because we calculate the bayesian probs
  ### later for the multi-panel plot
  outData <- parallel::parLapply(l, processLambda, oneTS, base, dates,
                                 window, runType, windowType,
                                 pointids, dataPath, 
                                 staticInflation, bayes=FALSE, allProbs=NA, 
                                 probType="total", ewmaOnly=TRUE, cl=cl)
  stopCluster(cl)
  
  # using histogram of the prob values per disturbed and undisturbed points,
  ## calculate the chance of getting a certain prob given there has or has not
  ## been a disturbance
  # probBreaks <- probHist(outData, plotHist, histInterval)
  # fwrite(probBreaks, filePath)
  
  if(testNormal){
    layout(matrix(1:20, ncol=5, byrow=TRUE))
    sapply(1:20, testNorm, outData, l)
    layout(matrix(1:1))
  }
  
  if(createFuns) saveF <- TRUE else saveF <- FALSE
  fileNames <- sapply(1:length(outData), saveDensityEst, l, outData, nRun, 
                      filePath, stable=saveF, monitor=saveF)
  rownames(fileNames) <- c("stable", "monitor")
  
  return(print(paste0("Done with run ", nRun, ". Please copy .Rdata file for",
                      " chosen landscape lambda to ",
                      "SEAL/data/myanmar/trainingPars/")))
}

dirPath <- "data/trainingPars/train3_densFun"
if(!dir.exists(dirPath)) dir.create(dirPath)
filePath <- paste0(dirPath, "/fun.Rdata")

sapply(c(2:3), getLikelihood, ewmaDays=TRUE, filePath, funs, testNormal=FALSE, 
      createFuns=TRUE)