##########################################################
## Purpose: Make metadata plot for paper1
## Run medium:
##  - Mac or PC
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.2.2, Jan 2023
##########################################################
# TO RUN THIS SCRIPT AND MAKE PLOT,
## need to first run train7_plotAccuracy with plotF1=TRUE and plotMeta=TRUE,
## such that these are in the global environment

library(groundhog)
groundhog.library(data.table, "2023-07-31")
groundhog.library(plotrix, "2023-07-31")
dt <- fread("data/metadata/methodsComparePlot.csv")

prepMeta <- function(dt, constrain){
  dtNRT <- dt[grepl("NRT", type)][-c(1,2), ]
  col <- viridis::viridis(nrow(dtNRT)+2)
  dtNRT[, `:=` (st = paste0("scenario", 1:.N), col = col[1:(length(col)-2)])]
  
  prep <- dtNRT[, .(timeMin, timeMean, timeMedian, f1, overallAcc, users, st, 
                    col)]
  
  prep[, `:=` (typeTime = ifelse(is.na(timeMean) & 
                                   !is.na(timeMedian), "median", "mean"),
               typeAcc = ifelse(!is.na(f1), "f1",
                                ifelse(!is.na(overallAcc), "overall", "users")))
  ][, `:=` (time = as.numeric(gsub("NA", "", paste0(timeMean, timeMedian))),
            acc = as.numeric(gsub("NA", "", paste0(f1, overallAcc, users))))]
  
  lookup <- data.table(typeTime=c("mean", "median", "median", "mean", "median"),
                       typeAcc=c("overall", "f1", "overall", "f1", "users"),
                       pch=c(1, 17, 2, 16, 25))
  prep[, pch := lookup$pch[match(paste0(typeTime, typeAcc), 
                                 paste0(lookup$typeTime, lookup$typeAcc))]]
  
  if(constrain){
    # prep[, timePlot := ifelse(time==60, 41, 
    #                           ifelse(time==72, 42,
    #                                  ifelse(time==126, 45, time)))]
    prep[, timePlot := ifelse(time >= 60 | is.na(time), 43, time)]
  } else {
    prep[, timePlot := time]
  }
  return(prep)
}
prepResults <- function(results){
  g <- results$lagsFull
  medLag <- g[, median(value, na.rm=TRUE), by=.(lambda, sensor)]
  colnames(medLag)[3] <- "lagMed"
  
  f <- results$tabF1PR
  medF1 <- f[, median(f1), by=.(lambda, sensor)]
  colnames(medF1)[3] <- "f1Med"
  
  medF1[, lagMed := medLag$lagMed[match(paste0(lambda, sensor),
                                        paste0(medLag$lambda, medLag$sensor))]]
  # cl <- colorRampPalette(c("violet", "orange", "gold"))
  medF1[, col := c(rep("#383aaa", nrow(medF1)/2),
                   rep("#CC33FF", nrow(medF1)/2))]
}

meta <- prepMeta(dt, constrain=TRUE)
medDT <- prepResults(results)

## first image: take the best F1 score per lambda and plot with associated
## median lag
tiff(paste0(imgFolder, "figXX_metadata1.tiff"), width=15, height=12, 
    units="cm", res=350, compression="lzw")
par(mar=c(5,4,2,2))
plot(NA, ylim=range(meta$acc), xlim=c(0, max(meta$timePlot)+2),
     ylab="Accuracy (%)", xlab="Detection lag (days)", xaxt="n")
axis(1, at=c(seq(0, 50, 10)), labels=c("0", "10", "20", "30", "40", "50"))
# axis.break(1,36,style="slash") 

points(meta$timePlot, meta$acc, pch=16)
sapply(1:16, function(s){
  label <- s
  pos <- 4
  if(s %in% c(10,14)){
    pos <- 1
  } else if(s==5){
    label <- paste0(s, "*")
  } else if(s %in% c(6,8)){
    pos <- 2
  } else if(s %in% c(9,16)){
    label <- paste0(s, "^")
  }
  text(meta[s,timePlot], meta[s,acc], label=label, pos=pos)
})

sapply(c("L8S2", "All"), function(s){
  pd <- medDT[sensor==s, ][, lambda := as.character(lambda)]
  focalLam <- "0.6"
  
  if(s=="L8S2") pch <- 17 else pch <- 15
  
  points(pd[lambda==focalLam, lagMed], pd[lambda==focalLam, f1Med], pch=pch, 
         col=unique(pd$col), cex=1)
  
  ## range of possible latencies over all lambda (horiz line where lambda=0.6)
  segments(x0=min(pd$lagMed), x1=max(pd$lagMed), col=unique(pd$col), 
           y0=pd[lambda==focalLam, f1Med], lty=1)
  
  ## range of possible F1s over all lambda (vertical line where lambda=0.6)
  segments(y0=min(pd$f1Med), y1=max(pd$f1Med), col=unique(pd$col), 
           x0=pd[lambda==focalLam, lagMed], lty=1)
  
  # if(s=="L8S2") pch <- 2 else pch <- 5
  # points(pd$lagMed, pd$f1Med, pch=pch, 
  #        col=adjustcolor(pd$col, alpha=0.5))
})

legend("bottomright", legend=c("MO", "MM", "Previous study"), bty="n",
       pch=c(17,15,16), col=c("#383aaa", "#CC33FF", "black"))
dev.off()
#

## second image: same as above but only plot 
png(paste0(imgFolder, "figXX_metadata2.png"), width=15, height=12, 
    units="cm", res=350)
par(mar=c(5,4,2,2))
plot(meta$time, meta$acc, pch=meta$pch, col="black", xlim=c(0, max(meta$time)),
     ylab="Accuracy (%)", xlab="Detection lag (days)")
sapply(c("L8S2", "All"), function(s){
  pd <- medDT[sensor==s & lambda==0.6, ]
  points(pd$lagMed, pd$f1Med, pch=8, 
         col=adjustcolor(pd$col, alpha=0.8))
})

legend("topright", legend=c("MO", "MM", "Mean-F1", "Median-F1", 
                            "Mean-overall", "Median-overall"), bty="n",
       pch=c(8,8,16,17,1,2), col=c("orange", "blue", rep("black", 4)))
dev.off()


############################################################################### 
g <- results$lags
f <- results$tabF1PR
f[, day := rep(1:60, nrow(f)/60)]
# g[, min(value, na.rm=TRUE), by=.(lambda, sensor)] # this is all 0, so do below
test <- f[day==1, max(f1), by=.(sensor)]

medResults <- lapply(c("L8S2", "All"), function(Y){
  f <- f[sensor==Y]
  
  g1 <- g[sensor==Y, quantile(value, probs=c(0.1, 0.9), na.rm=TRUE), by=lambda]
  g1[, prob := rep(c(0.1,0.9), .N/2)]
  setnames(g1, old="V1", new="day")
  g1[, f1 := f$f1[match(paste0(lambda,"-", day), 
                        paste0(f$lambda,"-", f$day))]]
  
  if(any(is.na(g1$f1))){
    tab <- g1[is.na(f1)]
    avg <- sapply(1:nrow(tab), function(X){
      dayUp <- f[lambda==tab[X, lambda] & day==ceiling(tab[X, day]), f1]
      dayDown <- f[lambda==tab[X, lambda] & day==floor(tab[X, day]), f1]
      return(mean(c(dayUp, dayDown), na.rm=TRUE))
    })
    g1[is.na(f1), ]$f1 <- avg
  }
  return(g1[, sensor := Y])
})
medResults <- rbindlist(medResults)

low <- medResults[prob==0.25]
high <- medResults[prob==0.75]
points(low$day, low$f1, pch=8, col=c(rep(adjustcolor("black", alpha=0.4), 20),
                                     rep(adjustcolor("blue", alpha=0.4), 20)))
points(high$day, high$f1, pch=8, col=c(rep(adjustcolor("black", alpha=0.4), 20),
                                       rep(adjustcolor("blue", alpha=0.4), 20)))
# points(rep(1, nrow(test)), test$V1, col=c("black", "blue"), pch=9)

## now plot the literature metadata
points(dtFull$timeMin, dtFull$f1_overallAcc, col=dtFull$col, pch=16)
points(dtFull$timeMin, dtFull$f1_overallAcc, pch=1, 
       col=adjustcolor("black", alpha=0.7))
legend("bottomright", legend=c(unique(dtNRT$year), 
                               "This study (SAR)", "This study (no SAR)"), 
       pch = c(rep(16, 7), 8,8), col=c(unique(dtNRT$col), "blue", "black"), 
       bty="n")

