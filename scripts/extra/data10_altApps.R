##########################################################
## Purpose: Identify alternate regions (outside Chatthin) to run the 
##          landscape code over, ideally covering different strata
## Run medium:
##  - Mac or PC
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Aug 2022
## Last modified: Aug 2022
##########################################################
library(terra)
library(data.table)
source("scripts/funs/dataC_idRegions.R")
crsProj <- "EPSG: 32646"

# step 1
# First, identify the S2 tiles that contain training points along with solid
## clusters (as good as possible) of the main forest strata. I do this to avoid
## having to download new S2 tiles and process with sen2cor.

b <- rast("data/dataMyanmar/spatial/forestMapMasked.tif")
plot(b, col=rev(viridis::viridis(50)))

s2 <- vect("data/spatial/s2Tiles/sentinel_2_index_shapefile.shp")
s2 <- project(s2, crsProj)
s2 <- crop(s2, b)
plot(s2, add=TRUE)
text(s2)

pts <- fread("data/dataMyanmar/trainingDataPoints.csv")
pts <- vect(pts, geom=c("coordX", "coordY"), crs="EPSG: 4326")
pts <- project(pts, crsProj)
plot(pts, add=TRUE)

#123 = lowland ever
#80 = upland ever
#43 = dry dec and mix dec
#53 = dry dec and mix dec
tiles <- s2$Name
targetTiles <- c(tiles[123], tiles[80], tiles[43], tiles[53])

appExtents <- lapply(targetTiles, getBasePoint, nCol=2773, nRow=2135)
# save(appExtents, file="data/dataMyanmar/appExtents.Rdata")
################################################################################
# Step 2
# Now, bring in the chatthin numbers and append to the appExtents list.
# absolute ext numbers for chatthin come from below
## chatthin #46QGM
chatthin <- terra::vect(shapePath)
chatthin <- terra::project(chatthin, plotCRS)
extent <- as.vector(terra::ext(terra::buffer(chatthin, width=100)))

extent <- c(744424.3, 772144.6, 2599900.5, 2621249.6) #xmin xmax ymin ymax

appExtents <- append(appExtents, list(extent))
names(appExtents) <- c("lowEver", "upEver", "mixDry1", "mixDry2", "chatthin")

appLocations <- list(
  lowEver=list(ext=appExtents$lowEver, s2tile="47RLH", crs="EPSG: 32646"),
  upEver=list(ext=appExtents$upEver, s2tile="46RGP", crs="EPSG: 32646"),
  mixDry1=list(ext=appExtents$mixDry1, s2tile="46QFL", crs="EPSG: 32646"),
  mixDry2=list(ext=appExtents$mixDry2, s2tile="46QHL", crs="EPSG: 32646"),
  chatthin=list(ext=appExtents$chatthin, s2tile="46QGM", crs="EPSG: 32646"))
save(appLocations, file="dissertation/data/dataMyanmar/appLocations.Rdata")

print("Please copy the file to SEAL/Ian/dissertation/data/myanmar/")
################################################################################
# Step 3

# While we could avoid downloading new S2 tiles, we will need to download L8
## imagery. However, we can do this in a smarter / faster way by cropping the L8
## images directly in GEE and thus download smaller files. To do so, we need the
## extents of the newly identified regions.

# get buffered extents for downloading L8 data from GEE
library(terra)
load("data/dataMyanmar/appLocations.Rdata")
sapply(appLocations, buffExtent)

## space used (matrix of 2773 col x 2135 rows at 10m res) = 59203.55 ha
2773*2135 # n pixels
27730*21350 # sq m
(27730*21350) * 0.000001 # sq km
(27730*21350) * 0.0001 # ha

################################################################################
library(data.table)
library(terra)
# Step 4: Plot northern myanmar with forest type and training data
## Fig 1 in paper1
base <- fread("data/dataMyanmar/trainingPars/trainingDataPoints.csv")
plotForestType(base, figType="pub")

# step 5: Plot northern myanmar with forest type and app regions outlined
## Fig 3 in paper1 and fig 1 in paper2
load("data/dataMyanmar/appLocations.Rdata")
# nm <- paste0("Region", 1:5)
nm <- c("Region 1", NA, NA, NA, "Region 2")

basePlot <- "elev" # forest or elev
plotRegions(regionN=c(1,5), basePlot, regionNames=nm, plotCrop=TRUE)

