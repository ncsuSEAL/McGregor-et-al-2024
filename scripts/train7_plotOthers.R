##########################################################
## Purpose: Calculate accuracy metrics and plot
## Run medium:
##  - Mac or PC
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
##########################################################
groundhog.library(data.table, groundhogDate)
groundhog.library(ggplot2, groundhogDate)
groundhog.library(ggdark, groundhogDate)
groundhog.library(viridis, groundhogDate)

bayes <- FALSE
source("scripts/args.R")
source("scripts/funs/trainE_metricsCalc.R")
source("scripts/funs/trainF_accuracyCV.R")

script <- 3
imgFolder <- "figures/"
source("scripts/argsTrain.R", local=TRUE) #the error is ok; keep moving forward

##----------------------------------------------------------------------------##
# Summary numbers from results
## NOTE: need to run `train7_plotAccuracy.R` first

b <- results$tabF1PR
b[lambda %in% c(0.05, 1), .(f1, lambda, sensor)][, min(f1), by=.(lambda, sensor)]
b[lambda %in% c(0.05, 1), .(f1, lambda, sensor)][, max(f1), by=.(lambda, sensor)]

b <- results$lags
b[lambda %in% c(0.05, 1), ][, min(value, na.rm=TRUE), by=.(lambda, sensor)]
b[lambda %in% c(0.05, 1), ][, max(value, na.rm=TRUE), by=.(lambda, sensor)]

b <- results$fpr
b[lambda %in% c(0.05, 1), .(fpr, lambda, sensor)]

b <- results$fnr
b[lambda %in% c(0.05, 1), .(fnr, lambda, sensor)]
##----------------------------------------------------------------------------##
# Make split size plots for FPR, FNR, and lags
plots <- c("fpr", "fnr", "lags")
plotData <- lapply(1:2, function(sn){
  if(sn==1){
    splitNum <- c("at most", 0.9)
  } else {
    splitNum <- c("at least", 1)
  }
  
  # splitNum = "at least" or "at most", 2nd argument is a percentage
  results <- plotAccuracy(f=metricsFiles, foldPath, randomLam=FALSE, consecObs, 
                          plotF1=FALSE, plotPR=FALSE, plotMeta=TRUE, base, 
                          splitSize=TRUE, splitNum=splitNum, 
                          plotWindow=window, windowType,
                          savePNG=FALSE, pngPars,
                          imgFolder="figures/", bayes=FALSE,
                          probType="total")
  
  out <- lapply(plots, function(X){
    dt <- results[[X]]
    return(dt[, size := splitNum[2]])
  })
  names(out) <- plots
  
  return(out)
})

imgFolder <- "figures/"
png(paste0(imgFolder, "figXX_accLag_", windowType, "_splitSize.png"),
    width=pngPars[img=="accLag", width], height=pngPars[img=="accLag", height],
    units=unique(pngPars$units), res=unique(pngPars$res))

layout(matrix(1:4, ncol=2, byrow=TRUE))
sapply(plots, function(X){
  if(X=="fpr"){
    ylab <- "Overall FPR"
  } else if(X=="fnr"){
    ylab <- "Overall FNR"
  } else if(X=="lags"){
    ylab <- "Detection lag (days)"
  }
  
  dt <- rbind(plotData[[1]][[X]], plotData[[2]][[X]])
  
  if(X != "lags"){
    par(mar=c(2,4,2,2))
    dt[, `:=` (col = ifelse(sensor=="L8S2", "purple", "red"),
               pch = ifelse(size==1, 16, 1),
               lty = ifelse(size==1, 1, 2))]
    
    plot(dt$lambda, unlist(dt[, ..X]), col=dt$col, pch=dt$pch, cex=0.7, xlab="", 
         ylab=ylab)
    
    sapply(unique(dt$size), function(s){
      subdt <- dt[size==s]
      for(i in unique(subdt$sensor)){
        lines(subdt[sensor==i, lambda], unlist(subdt[sensor==i, ..X]), 
              col=subdt[sensor==i, col], lty=subdt[sensor==i, lty])
      }
    })
    
    if(X=="fnr"){
      legend("topleft", 
             legend=c("MM >= 100%", "MM <= 90%", 
                      "MO >= 100%", "MO <= 90%"), 
             col=c("red", "red", "purple", "purple"),
             pch=c(16, 1, 16, 1),
             lty=c(1, 2, 1, 2), bty="n")
    }
  } else {
    par(mar=c(4,4,2,2))
    dt[, `:=` (col = ifelse(sensor=="L8S2", "purple", "red"),
               pch = ifelse(size==1, 4, 3),
               lty = ifelse(size==1, 1, 2))]
    
    sapply(unique(dt$sensor), function(sens){
      subdt <- dt[sensor==sens, ]
      
      plot(NA, xlim=range(dt$lambda), ylim=c(0,15),
           xlab="Lambda", main="", ylab=ylab)
      
      sapply(unique(dt$size), function(s){
        subdt <- subdt[size==s]
        for(i in unique(subdt$sensor)){
          for(j in unique(subdt$variable)){
            b <- subdt[sensor==i & variable==j, ]
            lines(b$lambda, b$med, col=b$col, lty=b$lty)
            points(b$lambda, b$quant, col="black", pch=b$pch)
          }
        }
      })
    })
    legend("topright", 
           legend=c("Median value", "", "25th percentile:", 
                    ">= 100%", "<= 90%"),
           lty=c(1, NA, NA, NA, NA),
           pch=c(NA, NA, NA, 4, 3),
           col=c("black", NA, NA, "black", "black"), bty="n")
  }
})
dev.off()

par(mar=c(5,4,4,2) + 0.1)
layout(matrix(1:1))

##----------------------------------------------------------------------------##
# ECDF plot of disturbance sizes
## Chapter 1 Figures 7 and S7
b <- base[!grepl("a", pointid)]
plot(ecdf(b$distPercL8), verticals=TRUE, xaxt="n", 
     xlab="Size of disturbance (m^2)", ylab="% of training disturbances",
     main="Distribution of disturbance sizes")
axis(1, at=seq(0.2,1,0.2), labels=c(900*seq(0.2,0.8,by=0.2), ">=900"))
##----------------------------------------------------------------------------##
# Plot ECDF of lag times comparing the reliance on only 1 observation vs
## 3 observations
consecObs <- c(1,3)

layout(matrix(1:2, ncol=2))

### l = lambda
### sensorCombo = "All" or "L8S2"
plotMedianLags <- function(obs, l, sensorCombo, window, windowType, 
                           plotBox, plotECDF, probType){
  foldPath <- paste0("data/trainingPars")
  
  ## obsDates
  fDates <- list.files(foldPath, pattern="obsDates", full.names=TRUE)
  
  ## metrics
  fMetrics <- list.files(paste0(foldPath, "/train5_allMetrics"), 
                         pattern="Metrics_", full.names = TRUE)
  fMetrics <- fMetrics[!grepl("bayes", fMetrics) & 
                         grepl(paste0("conObs", obs), fMetrics) & 
                         grepl(windowType, fMetrics)]
  
  out <- plotLag(fDates, fMetrics, l, sensorCombo, foldPath, window, windowType, 
                 obs, plotBox, plotECDF)
  return(out)
}

## ecdf plot just for L8S2
png(paste0(imgFolder, "figXX_ecdfBasic.png"), width=16, height=10, units="cm",
    res=350)
layout(matrix(1:2, ncol=2))
metrics <- lapply(consecObs, plotMedianLags, l=0.6, sensorCombo="L8S2", 
                  window=60, windowType, plotBox=FALSE, plotECDF=TRUE,
                  probType="total")
par(mar=c(5,4,4,2)+0.1)
dev.off()

median(metrics[[1]]$lagFast)
median(metrics[[1]]$lagSlow)

## ecdf plot for both sensor combinations
png(paste0(imgFolder, "figXX_ecdfBothSensor.png"), width=16, height=12, 
    units="cm", res=350)
layout(matrix(1:4, nrow=2, byrow=TRUE))
for(s in c("L8S2", "All")){
  metrics <- lapply(consecObs, plotMedianLags, l=0.6, sensorCombo=s, 
                    window=60, windowType, plotBox=FALSE, plotECDF=TRUE)
}
dev.off()

par(mar=c(5,4,4,2)+0.1)
layout(matrix(1:1))

# Tangentially, plot lag by size class.
## Here is consecObs=1 for an example (not the best plot, not sure if necessary)
foldPath <- paste0("data/trainingPars/")
f <- list.files(foldPath)
f <- f[grepl("metrics|obsDates", f)]
f <- f[grepl("conObs1", f)]
plotLagSizeClass(f, l, sensorCombo = "L8S2")

##----------------------------------------------------------------------------##
# Get mean or median time (days) over observations
consecObs <- 1
foldPath <- paste0("trainingPars/")
f <- list.files(foldPath)
type <- "mean"
datesFiles <- f[grepl("obsDates", f) & !grepl("S1", f) & 
                  grepl(paste0("conObs", consecObs), f)]
daysByObs <- avgDaysPerObs(datesFiles, backfill=FALSE, window, type)

# obs corresponding with window
daysByObs[time >= window, rep, by=.(sensor)][, min(rep), by=.(sensor)]

library(ggdark)

cols <- c("#0000FF", "#33CCCC", "#00CC66")
ggplot(daysByObs) +
  aes(x = rep, y = time, color = sensor) +
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=cols) +
  scale_x_continuous(breaks = function(x){
    unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))}) +
  labs(title=paste0(type, " days per sensor observation"),
       x="Observations since disturbance",
       y="Days since disturbance") +
  geom_hline(yintercept=60, linetype="dashed") +
  dark_theme_classic()
  # theme_dark()