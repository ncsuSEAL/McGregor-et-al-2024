##########################################################
## Purpose: Functions to help identify other landscape regions
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, July 2022
## Last modified: Aug 2022
##########################################################

## -----------------------------------------------------------------------------
## getBasePoint =  identify centroid of new region and convert to extent vector
## buffExtent = take basic extent and buffer by 1km, then project to lat/lon.
##              This is so that we can directly use the coordinates in GEE
## plotRegions = plot alternate region locations over northern myanmar
## createBorderShp = create vector shapefile of landscape region extent (rect)
## -----------------------------------------------------------------------------
getBasePoint <- function(X, nCol, nRow){
  tile <- crop(s2, s2[s2$Name==X])
  plot(tile)
  plot(b, add=TRUE, col=rev(viridis::viridis(50)))
  plot(tile, add=TRUE)
  out <- click()
  
  centroid <- as.numeric(out)
  
  #xmin, xmax, ymin, ymax
  fullEx <- c(centroid[1]-((nCol*10)/2), centroid[1]+((nCol*10)/2),
              centroid[2]-((nRow*10)/2), centroid[2]+((nRow*10)/2))
  
  return(fullEx)
}
buffExtent <- function(X){
  Xext <- X$ext
  xmin <- Xext[1] - 1000
  xmax <- Xext[2] + 1000
  ymin <- Xext[3] - 1000
  ymax <- Xext[4] + 1000
  
  pts <- data.frame(lon=c(xmin, xmin, xmax, xmax), 
                    lat=c(ymin, ymax, ymax, ymin))
  shp <- vect(pts, crs="EPSG: 32646")
  shp <- project(shp, "EPSG: 4326")
  return(ext(shp))
}
plotRegions <- function(regionN, basePlot, regionNames, plotCrop){
  path <- "data/dataMyanmar/spatial/"
  rastImg <- rast(paste0(path, "forestMapMasked.tif"))
  
  if(basePlot=="forest"){
    plot(rastImg, col=rev(viridis::viridis(50)), type="classes",
         levels=c("Lowland evergreen", "Upland evergreen", 
                  "Mixed deciduous", "Dry deciduous"))
  } else if(basePlot=="elev"){
    rastEl <- rast(paste0(path, "covariates/dem/demMMR_500.tif"))
    bor <- vect(paste0(path, "/covariates/MMR_adm0.shp"))
    rastEl <- mask(rastEl$elevation, bor)
    rastEl <- project(rastEl, "EPSG:32646")
    plot(rastEl)
    sbar(100000, xy="bottomleft", label=c("", "100 km", ""))
  }
  
  pol <- lapply(regionN, function(X){
    extent <- appLocations[[X]]$ext
    test <- t(rbind(c(extent[1], extent[1], extent[2],extent[2]),
                    c(extent[3], extent[4], extent[4],extent[3])))
    colnames(test) <- c("x", "y")
    polySq <- terra::vect(test, type="polygons",
                          crs=appLocations[[X]]$crs)
    polys(polySq)
    text(polySq, labels=regionNames[X], halo=TRUE, pos=3, cex=1.5, offset=1)
    return(polySq)
  })
  
  if(plotCrop){
    sapply(1:length(pol), function(X){
      if(X==1){
        pax <- list(side=c(1:2))
        xy <- c(920000, (2750000-625))
      } else if(X==2){
        pax <- list(side=c(1:2), 
                    yat=seq(2605000, 2620000, by=5000),
                    ylabs=seq(2605000, 2620000, by=5000))
        xy <- "bottomleft"
      }
      
      png(paste0("writings/paper1/figures/figXX_forestType", X, ".png"), 
          res=350, width=14.7, height=10, units="cm")
      plot(crop(rastImg, pol[[X]]), 
           col=rev(c("#FFFFCC", "#CCCC66", "#006633", "#33CC00")), 
           type="classes", colNA="#333333", pax=pax,
           plg=list(cex=0.7), 
           levels=c("Lowland evergreen", "Upland evergreen", "Mixed deciduous", 
                    "Dry deciduous"))
      sbar(2500, xy=xy, label=c("", "2.5 km", ""), 
           col="white", cex=0.8)
      dev.off()
    })
  }
  return("Done yay")
}
plotForestType <- function(base){
  path <- "data/dataMyanmar/spatial/"
  rastImg <- rast(paste0(path, "forestMapMasked.tif"))
  
  bound <- vect(paste0(path, "covariates/MMR_adm0.shp"))
  pts <- vect(base, geom=c("coordX", "coordY"), crs="EPSG: 4326")
  pts <- project(pts, rastImg)
  bound <- project(bound, rastImg)
  
  rastImg <- crop(rastImg, bound)
  
  png(paste0("writings/paper1/figures/figXX_myanmarBorder.png"), 
      res=350, width=6, height=10, units="cm")
  plot(bound)
  dev.off()
  
  png(paste0("writings/paper1/figures/figXX_forestTypeTraining.png"), 
      res=350, width=16, height=12, units="cm")
  plot(rastImg, 
       col=rev(c("#FFFFCC", "#CCCC66", "#006633", "#33CC00")), 
       type="classes", 
       plg=list(cex=0.7), 
       levels=c("Lowland evergreen", "Upland evergreen", "Mixed deciduous", 
                "Dry deciduous"))
  points(pts, col="orange", pch=16)
  points(pts, col="black", pch=1)
  dev.off()
}
createBorderShp <- function(region, crs="EPSG:4326"){
  load("dissertation/data/myanmar/appLocations.Rdata")
  num <- appLocations[[region]]$ext
  
  m <- matrix(c(num[1], num[3],
                num[1], num[4],
                num[2], num[3],
                num[2], num[4]),
              ncol=2, byrow=TRUE)
  v <- vect(m, crs="EPSG:32646")
  e <- ext(v)
  p <- as.polygons(e, crs="EPSG:32646")
  p <- project(p, crs=crs)
  plot(p)
  writeVector(p, 
              paste0("dissertation/data/myanmar/appLowEver/region", region,
                     "Border.shp"), overwrite=TRUE)
  return(print("Done. yay"))
}