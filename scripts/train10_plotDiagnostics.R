##########################################################
## Purpose: Make multi-panel summary plots (ts models, probs, etc) along with
##          PLANET images for the disturbed points. Also compare distributions
##          of zscores between disturbed and undisturbed points.
## Run medium:
##  - Mac
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
##########################################################
library(groundhog)
groundhogDate <- "2023-07-31"
groundhog.library(data.table, groundhogDate)
groundhog.library(terra, groundhogDate)
groundhog.library(parallel)
groundhog.library(ggplot, groundhogDate)
groundhog.library(scales, groundhogDate) #for histogram of disturbance dates

# NOTE
## In order to run this successfully, you MUST have run script 
## train1_getResiduals.R with the arguments `returnHarMod` and `makeSumPlots`
## set to TRUE. These are found in args.R

# Basic variables needed
sensKeep <- c("landsat8, sentinel2")
script <- 1
bayes <- FALSE
source("scripts/args.R")
source("scripts/argsTrain.R")
source("scripts/funs/trainC_byLambda.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")
source("scripts/funs/commonE_plotSummaries.R")
source("scripts/funs/commonF_probBayes.R")
load("data/trainingPars/train1_L8S2.Rdata")

lambdaVec <- c(0.1, 0.5, 1) # which lambdas do we want to plot?

ptsNoDist <- pointids[grepl("a", pointids)]
ptsDist <- pointids[!grepl("a", pointids)]

load(paste0(dataPath, "/trainingPars/allStaticProbs.Rdata"))
allProbs <- cbind(base[, .(pointid)], 
                  data.table(prior=allProbs[,1]/10000))

##--------------------------------------------##
# Actually run the plotting code

## if planet=TRUE, then pathPlanet is set in `args.R`

# individual plot
plotDiagnostic(pt="g0p0_1", oneTS, lambdaVec, window, runType, 
               base, planet=FALSE, vegIndex, sensKeep, logMod=NULL,
               bayes=TRUE, allProbs)

## Writing everything to pdf will take several minutes
pdf(paste0("diagnosticPlotsChat_", fill, ".pdf"), height=10, width=12)
sapply(pts, plotDiagnostic, oneTS, lambdaVec, window,  
       runType, base, planet=FALSE, vegIndex, sensKeep)
dev.off()

pdf(paste0("diagnosticPlotsG10end_", fill, ".pdf"), height=10, width=12)
pbapply::pbsapply(pts, plotDiagnostic, oneTS, lambdaVec, window, 
       runType, base, planet=FALSE, vegIndex, sensKeep)
dev.off()

###############################################################################
# Histogram of disturbance dates
layout(matrix(1:1))
breakpoints <- sapply(1:12, function(X){
  return(as.Date(paste0("2019-", X, "-01")))
})

png("figures/figXX_distByMonth.png", width=20, height=16,
    units="cm", res=350)
par(mar=c(5,4,2,1))
hist(base[!grepl("a", pointids), dateDist], freq=FALSE,
     breaks=c(as.Date(breakpoints, origin=as.Date("1970-01-01")), 
              as.Date("2020-01-01")),
     col=c(alpha("orange", 0.4), 
           rep(alpha("blue", 0.4), 4), 
           rep(alpha("purple", 0.4), 4), 
           rep(alpha("orange", 0.4), 3)),
     xlab="Month from Jan 2019 - Jan 2020", main=""
     # main="Distribution of disturbances by month"
)
legend("topright", legend=c("Cool dry season", "Hot dry season", "Wet season"),
       col=c(alpha("orange", 0.4), alpha("blue", 0.4), alpha("purple", 0.4)),
       pch=15, bty="n", pt.cex=2)
dev.off()
par(mar=c(5,4,4,2) + 0.1)
###############################################################################
# Met data
## Plot met data from GEE to support evidence of having a short stable period
# precip, pdsi, soil moisture, water deficit
plotMet <- function(m, metDT, metVars, mets, units, labs, labsY){
  # note that mean doesn't actually do anything here, but gets us the 3 columns
  ## we need
  test <- metDT[, mean(get(metVars[m])), by=.(date, region)]
  colnames(test)[3] <- metVars[m]
  
  if(metVars[m]=="pdsi"){
    ylab=mets[m]
  } else {
    ylab=paste0(mets[m], " [", units[m], "]")
  }
  
  plot(NA, xlim=range(test$date), ylim=range(test[, get(metVars[m])]), xaxt="n",
       xlab="", ylab="")
  title(ylab=ylab, line=2)
  if(m %in% c(2,4)){
    axis(1, at=c(as.Date("2016-01-01"), as.Date("2017-01-01"), 
                 as.Date("2018-01-01"), as.Date("2019-01-01"), 
                 as.Date("2020-01-01")),
         labels=seq(2016,2020, 1))
  } else {
    axis(1, at=c(as.Date("2016-01-01"), as.Date("2017-01-01"), 
                 as.Date("2018-01-01"), as.Date("2019-01-01"), 
                 as.Date("2020-01-01")),
         labels=FALSE)
  }
  
  if(m==1){
    legend("bottomleft", legend=c("Region 1", "Region 2"), 
           col=c("black", "blue"), lty=c(1,3), bty="n")
  }
  
  sapply(unique(test$region), function(X){
    if(X==1){
      col <- "black" 
      lty <- 1
    } else {
      col <- "blue"
      lty <- 3
    }
        
    dtSub <- test[region==X, ]
    lines(dtSub$date, dtSub[, get(metVars[m])], col=col, lty=lty)
  })
  abline(v=as.Date("2019-01-01"), lty=2)
  text(as.Date("2019-12-01"), labsY[m], labels=labs[m], cex=2)
}

metDT <- rbindlist(lapply(c(1,5), function(r){
  dt <- fread(paste0("data/sensorData/metData/metDataRegion",
                        r, ".csv"))
  dt[, `:=` (grp = gsub("[[:digit:]]*_", "", `system:index`),
                date = as.Date(date, format="%m/%d/%Y"),
                region = r)]
  return(dt)
}))

png("figures/figXX_metData.png", width=16, height=12, res=350,
    units="cm")
layout(matrix(1:4, ncol=2))
metVars <- c("pdsi", "pr", "soil", "def")
labs <- c("a", "c", "b", "d")
labsY <- c(275, 440, 2950, 1175)
units <- c("", "mm", "mm", "mm")
mets <- c("PDSI", "Precipitation", "Soil moisture", "Water deficit")
sapply(1:length(metVars), function(m){
  if(m %in% c(1,3)){
    par(mar=c(1,3,2,1))
  } else {
    par(mar=c(3,3,1,1))
  }
  
  plotMet(m, metDT, metVars, mets, units, labs, labsY)
})
dev.off()

###############################################################################
# Size classes of training data
distPts <- base[pointid %in% ptsDist]

## simple ECDF plot
plot(ecdf(distPts$distPercL8), verticals=TRUE, xaxt="n",
     ylab="Cumulative % of training data",
     xlab="Size of disturbance (m^2)",
     main="Distribution of disturbance sizes")
axis(1, at=seq(0.2, 1, by=0.2), 
     labels=c(as.character(seq(0.2, 0.8, by=0.2)*900), ">=900"))

## complicated ECDF with barplot included
subDT <- distPts[, .N, by=distPercL8
               ][order(distPercL8, decreasing=TRUE)
                 ][, percent := round(prop.table(N), 3)
                   ][, fullPerc := cumsum(percent)]

par(mar=c(5,4,4,4))
ylim <- c(0,200)
barplot(subDT$N, ylim=ylim, space=0, yaxt="n", 
        main="Distribution of disturbance sizes", ylab="Number of training points",
        xlab="Size of disturbance (m^2)")
lines(seq(0.5, 8.5, by=1), subDT$fullPerc*ylim[2], col="red", lwd=2)
abline(h=0.75*ylim[2], lty=2, col="red")
axis(1, at=seq(0.5, 8.5,by=1), labels=c(">= 900", seq(0.9, 0.2, by=-0.1)*900), 
     line=0)
axis(2, at=seq(0, ylim[2], by=20), labels=seq(0, ylim[2], by=20), line=0)
axis(4, at=seq(0, ylim[2], by=40), labels=seq(0, 1, length.out=6), col="red", 
     col.axis="red")
mtext("Cumulative % of training data", side = 4, line = 3, col="red")

###############################################################################
# Chapter 1, Figure S4
## Distribution of dist and undist zEWMA values by sample lambda
lam <- c(0.05, 0.5, 1)
cols <- c("purple", "#00CCFF", "#339900")

png("figures/figXX_distributions.png", width=20, height=14,
    res=350, units="cm")
layout(matrix(1:2, ncol=2))
for(q in c("L8S2", "All")){
  if(q=="All"){
    par(mar=c(5,2,2,1))
    ylim <- ""
    yaxt="n"
    main <- "Multi-source mixed"
  } else {
    par(mar=c(5,4,2,0))
    ylab <- "Density"
    yaxt="s"
    main <- "Multi-source optical"
  }
  
  # xlim <- c(-41, 28)
  plot(NA, xlim=c(-20, 5), ylim=c(0, 2.3), xlab="zEWMA", yaxt=yaxt,
       ylab=ylab, main=main)
  
  if(q=="All"){
    axis(2, at=seq(0,2,0.5), labels=FALSE)
  } else {
    legend("topleft", legend=c(as.character(lam), "Disturbed", "Undisturbed"),
           lty=c(1,1,1,1,2), col=c(cols, "black", "black"), bty="n", 
           lwd=rep(2, 5))
  }
  
  sapply(1:length(lam), function(X){
    path <- "data/trainingPars/train4_fullZewmaProbs/"
    dt <- fread(paste0(path, "train4_", q, "_", lam[X]*100, ".csv"))
    nodist <- rbind(dt[grepl("a", pointid)], 
                    dt[!grepl("a", pointid) & dist==0])
    dist <- dt[dist==1 & !grepl("a", pointid)]
    
    lines(density(nodist$ewma), col=cols[X], lty=2, lwd=2)
    lines(density(dist$ewma), col=cols[X], lty=1, lwd=2)
  })
}
dev.off()

###############################################################################
# Compare min and max z-score distributions between dist and undis points
test <- rbindlist(lapply(pointids, function(pt){
  aggRes <- oneTS[[pt]]$oneTS
  return(data.table(pointid=pt, minZ=min(aggRes$stdRes, na.rm=TRUE), 
                    maxZ=max(aggRes$stdRes, na.rm=TRUE)))
}))
layout(matrix(1:2, ncol=2))
plot(density(test[grepl("a", pointid), minZ]), col="blue", xlim=c(-50,0), 
     main="Distribution of minZ")
lines(density(test[!grepl("a", pointid), minZ]))
plot(density(test[grepl("a", pointid), maxZ]), col="blue", ylim=c(0,1),
     main="Distribution of maxZ")
lines(density(test[!grepl("a", pointid), maxZ]))
legend("topright", legend=c("disturbed", "undisturbed"), col=c("black", "blue"),
       lty=1, bty="n")
###############################################################################
# Plots showing number of observations by sensor, and per date
## (aggregated over all the training data)

fileNames <- c("S1", "L8S2", "All")
numRun <- 3

dataPath <- "data/"
fileLoadLoc <- paste0(dataPath, "/trainingPars/train1_", fileNames[numRun])
load(paste0(fileLoadLoc, ".Rdata"))

allDates <- lapply(names(oneTS), function(X){
  pt <- oneTS[[X]]$stdRes
  out <- rbindlist(lapply(names(pt), function(s){
    return(data.table(date=pt[[s]]$obsDates, sensor=s))
  }))
  out[, pt := X]
})
allDates <- rbindlist(allDates)

## remove those which had very low numbers, since models weren't fit to them
lose <- c(which(allDates$sensor=="sentinel1DES" & allDates$pt=="g0p0_8a"),
          which(allDates$sensor=="sentinel1ASC" & allDates$pt=="g10p0_68a"))
allDates <- allDates[-lose, ]

# Consolidate the data for plots
bl <- allDates[, .N, by=.(date, sensor)]

## add qualifiers
bl[, `:=` (dateQ = as.Date(date, origin=as.Date("1970-01-01")),
           nq = round(N/630, 2))]
bl[, sensorNew := ifelse(sensor=="landsat8sr", "Landsat-8",
                         ifelse(sensor=="sentinel2l2a", "Sentinel-2",
                                ifelse(sensor=="sentinel1ASC", 
                                       "Sentinel-1 Asc",
                                       "Sentinel-1 Des")))]
bl[, sensorNew := factor(sensorNew, levels=c("Sentinel-1 Des", "Sentinel-1 Asc",
                                             "Sentinel-2", "Landsat-8"))]
bl <- bl[order(sensorNew)]

# Boxplot of number of observations per sensor (full timeframe)
modStart <- as.numeric(as.Date("2019-01-01"))

## pars
cols <- viridis::viridis(4)
dateAx <- c(as.Date("2016-01-01"), as.Date("2017-01-01"), 
            as.Date("2018-01-01"), as.Date("2019-01-01"), 
            as.Date("2020-01-01"))

png("figures/figXX_sensorObsN.png", width=14, height=18,
    units="cm", res=350)
par(mar=c(4,2,2,1))
boxplot(N ~ sensorNew, data=bl, xlab="Number of observations", ylim=c(0,360),
        ylab="", horizontal=TRUE, col=alpha(cols, 0.8), boxwex=0.5)
dev.off()

# Now get the dates representing the median number of observations
png("figures/figXX_sensorObsTime.png", width=20, height=14,
    units="cm", res=350)
par(mar=c(4,2,2,1))
plot(NA, ylim=c(1,4), xlim=c(as.Date("2015-06-15"), as.Date("2020-01-31")),
     xlab="", yaxt="n", xaxt="n", ylab="")
axis(1, at=dateAx, labels=c("2016", "2017", "2018", "2019", "2020"))

sens <- unique(bl$sensorNew)
sapply(1:length(sens), function(X){
  most <- bl[sensorNew==sens[X] & N >= 100]
  points(most$dateQ, rep(X, nrow(most)), pch=15, cex=0.8, col=cols[X])
})

axis(2, at=c(1,2,3,4), labels=sens, las=0)
dev.off()

## number of points that had observation for each sensor over time
bl[, col := ifelse(sensor=="landsat8sr", cols[4],
                   ifelse(sensor=="sentinel2l2a", cols[3],
                          ifelse(sensor=="sentinel1ASC", cols[2], cols[1])))]
bl[, sensorNew := factor(sensorNew, levels=rev(levels(sensorNew)))]
g <- ggplot(bl, aes(x=date, y=nq, color=sensorNew)) +
  geom_bar(stat="identity") + 
  ylim(0,0.6) +
  scale_color_manual(values=rev(cols), guide="none") +
  xlab("") +
  ylab("Number of points that had observation post-filtering (%)") +
  scale_x_continuous(breaks=as.numeric(dateAx), 
                     labels=c("2016", "2017", "2018", "2019", "2020")) +
  theme_minimal() +
  facet_wrap(~sensorNew, ncol=1)

png("figures/figXX_sensorObsTimeAll.png", width=18, height=18,
    units="cm", res=350)
print(g)
dev.off()

###############################################################################
# Density of FP and TP for each lambda by ewma value
layout(matrix(1:1))
source("scripts/funs/trainC_byLambda.R")
source("scripts/funs/commonD_ewmaBinaryCat.R")
funs <- ls()

script <- 1
sensKeep <- c("landsat8, sentinel2")
source("scripts/args.R", local=TRUE)
source("scripts/argsTrain.R", local=TRUE)
fileNames <- "L8S2"

#Bring in the residuals
fileLoadLoc <- paste0(dataPath, "/gparetoSim/threshold", probThresh*100, 
                      "/", fileNames)
load(paste0(fileLoadLoc, ".Rdata"))

l <- c(0.1, 0.5, 1)
comb <- TRUE

## run the main code
parseVals <- function(X, probsList, iter){
  subList <- probsList[[paste0("lambda_", X)]]
  
  tabOrig <- rbindlist(subList)
  tabOrig <- tabOrig[!grepl("a", pointid)]
  
  if(iter=="density"){
    densDist <- density(tabOrig[dist==1, ewma])
    densNone <- density(tabOrig[dist==0, ewma])
    return(list(densDist=densDist, densNone=densNone))
  } else if(iter=="vals"){
    nodist <- as.numeric(tabAll[dist==0, probs])
    dist <- as.numeric(tabAll[dist==1, probs])
    return(list(nodist=nodist, dist=dist))
  }
}
plotDens <- function(lam, type, xNone, yNone, xDist, yDist, l, cols, dens, obs){
  if(type=="none"){
    xRange <- xNone
    yRange <- yNone
    main <- "Undisturbed"
  } else if(type=="dist"){
    xRange <- xDist
    yRange <- yDist
    main <- "Disturbed"
  } else if(type=="comb"){
    xRangeFull <- range(xNone, xDist)
    xRange <- c(-25, max(xRangeFull))
    yRange <- range(yNone, yDist)
    main <- "Distrib. of disturbed and undist ewma values per lambda"
  }
  
  if(lam==1){
    plot(NULL, xlim=xRange, ylim=yRange,
         xlab="EWMA value", ylab="Density",
         main=main)
    
    if(type=="comb"){
      legend("topleft", legend=c(l, "Disturbed", "Undisturbed"), 
             col=c(cols, "black", "black"), bty="n", lty=c(1,1,1,1,2))
    } else {
      legend("topleft", legend=l, col=cols, bty="n", lty=1)
    }
  }
  
  if(type=="none"){
    lines(dens[[lam]]$densNone$x, dens[[lam]]$densNone$y, col=cols[lam])
  }
  
  if(type=="dist"){
    lines(dens[[lam]]$densDist$x, dens[[lam]]$densDist$y, col=cols[lam])
  }
  
  if(type=="comb"){
    lines(dens[[lam]]$densNone$x, dens[[lam]]$densNone$y, col=cols[lam], lty=2,
          lwd=2)
    lines(dens[[lam]]$densDist$x, dens[[lam]]$densDist$y, col=cols[lam], lwd=2)
  }
  
  if(type=="comb") return(xRangeFull)
}
prepDensPlot <- function(obsNum, probsList, l, comb){
  le <- length(l)
  dens <- lapply(l, parseVals, probsList, iter="density")
  yNone <- c(0, max(sapply(1:le, function(X) return(max(dens[[X]]$densNone$y)))))
  xNone <- range(sapply(1:le, function(X) return(range(dens[[X]]$densNone$x))))
  
  yDist <- c(0, max(sapply(1:le, function(X) return(max(dens[[X]]$densDist$y)))))
  xDist <- range(sapply(1:le, function(X) return(range(dens[[X]]$densDist$x))))
  
  cols <- viridis::viridis(le)
  
  if(comb){
    sapply(1:length(l), plotDens, type="comb", xNone, yNone, xDist, yDist, l, 
           cols, dens, obsNum)
  } else {
    sapply(1:length(l), plotDens, type="none", xNone, yNone, xDist, yDist, l, 
           cols, dens, obsNum)
    sapply(1:length(l), plotDens, type="dist", xNone, yNone, xDist, yDist, l, 
           cols, dens, obsNum)
  }
}

if(!comb) layout(matrix(1:2, ncol=2))

## note that the 1 vs 3 obs only changes the first few ewma / prob values of each
## ts but that's it (and these are small changes), so overall the range and 
## density of things doesn't change.
visDens <- function(l){
  cl <- makeCluster(detectCores()-1)
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, ls(), envir=environment())
  clusterExport(cl, c(funs), envir=.GlobalEnv)
  probsAll <- parallel::parLapply(l, processLambda, oneTS, base, dates,
                                  window, plotProbModel=FALSE, 
                                  saveProbModel, probModelFile, 
                                  probThresh, runType, windowType,
                                  pointids, consecObs=1, ewmaOnly=TRUE, 
                                  staticInflation=staticInflation,
                                  cl=cl)
  stopCluster(cl)
  names(probsAll) <- paste0("lambda_", l)
  
  prepDensPlot(obsNum, probsList=probsAll, l, comb=TRUE)
}

visDens(l)
layout(matrix(1:1))
