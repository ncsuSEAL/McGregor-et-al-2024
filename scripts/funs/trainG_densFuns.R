createDensFuns <- function(dat, stable, monitor, saveFile, fileStable="", 
                           fileMonitor=""){
  if(stable){
    d0 <- c(dat[dist==0, ewma], dat[grepl("a", pointid) & dist==1, ewma])
    # hist(d0, breaks=50, probability=TRUE)
    # lines(density(d0))
    densEStable <- approxfun(density(d0, adjust=1.5, from=min(d0), 
                                     to=quantile(d0, 0.5)), rule=2, yleft=0)
    if(saveFile) save(densEStable, file=fileStable)
  }
  
  if(monitor){
    d1 <- dat[!grepl("a", pointid) & dist==1, ewma]
    
    # hist(d1, prob=TRUE, breaks=50)
    # lines(density(d1))
    
    densEDist <- approxfun(density(d1), rule=2, yright=0)
    if(saveFile) save(densEDist, file=fileMonitor)
  }
  return(list(densEStable=densEStable, densEDist=densEDist))
}
