##########################################################
## Purpose: Functions to mask forest in order to look at strata across training
##          data
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.3, July 2022
## Last modified: August 2022
##########################################################
maskNonForest <- function(rastPath, outPath, plotRast){
  land <- rast(rastPath)
  
  #these are the strata we want keep
  # lowEver <- land == 1 
  # upEver <- land == 2
  # decMix <- land == 4
  # decDry <- land == 6
  
  ## Note that masking like in this way just means we're creating two categories, 
  ## either 1 (the classes we want) or 0 (everything else). This doesn't retain
  ## the individual classes, however.
  # maskRast <- land %in% c(1,2,4,6)
  # out <- mask(land, maskRast)
  
  # Crop to only focus on northern Myanmar
  ## roughly equivalent to 20 degrees N
  landExt <- ext(land)
  landExt[3] <- 2250000 #change ymin (this is correct for 32647)
  
  print("Cropping to only northern Myanmar")
  out <- crop(land, landExt)
  
  # Do masking
  print("Masking non-forested areas")
  m <- c(0, NA,
         3, NA,
         5, NA,
         7, NA,
         8, NA)
  rclmat <- matrix(m, ncol=2, byrow=TRUE)
  out <- classify(out, rclmat)
  
  # Write to file
  print("Reprojecting to 46N (EPSG: 32646). This will take a second.")
  utm46 <- project(out, "EPSG: 32646", method="near")
  writeRaster(utm46, outPath, datatype="INT1U", overwrite=TRUE)
  
  print("Finished")
  
  if(plotRast) plot(utm46)
  return(print("Done"))
}
