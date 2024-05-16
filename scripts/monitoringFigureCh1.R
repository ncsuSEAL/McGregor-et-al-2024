##########################################################
## Purpose: Create Figure S10
## Run medium: Mac
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.3.1 April 2024
##########################################################
library(groundhog)
groundhog.library(terra, "2023-07-31")
groundhog.library(data.table, "2023-07-31")
path <- "data/spatial/"
demPath <- "covariates/dem/demMMR_merge.tif"

## NOTE: you will need to download the Hansen data from the GLAD repository
##      matching the file names below in order to run this script.

prepGLAD <- function(path, demPath){
    loss <- rast(paste0(path, "glad/Hansen_GFC-2021-v1.9_lossyear_30N_090E.tif"))
    cover <- rast(paste0(path, 
                        "glad/Hansen_GFC-2021-v1.9_treecover2000_30N_090E.tif"))
    bound <- vect(paste0(path, "covariates/MMR_adm0.shp")) #myanmar border
    dem <- rast(paste0(path, demPath)) #for cropping

    lossCrop <- crop(loss, bound)
    coverCrop <- crop(cover, bound)

    # mask to only keep pixels with >=30% tree cover
    coverCrop[coverCrop < 30] <- NA

    ## `masked` is a tif that only includes pixels that were at least 
    ## 30% tree cover in 2000
    masked <- mask(lossCrop, coverCrop, maskvalues=NA)

    ## we then mask to Myanmar's borders
    masked <- mask(masked, bound)
    masked <- crop(masked, dem)

    ## Retain disturbances only from monitoring year
    m <- c(0, 18, 0,              # mask undisturbed cells
        19, 19, 1,   # keep any disturbances over previous years
        20, 22, 0)   # mask this year + future years
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    dist19 <- classify(masked, rclmat, right=NA)

    return(dist19)
}

## get glad data, mask to myanmar border and >= 30% canopy cover and 2019 loss year
glad19 <- prepGLAD(path, demPath)

##

## show only dLoc
floc <- fread("data/trainingPars/trainingDataPoints.csv")
dloc <- floc[!grepl("a", pointid)]
uloc <- floc[grepl("a", pointid)]
v <- vect(dloc, geom=c("coordX", "coordY"), crs="EPSG:4326")
u <- vect(uloc, geom=c("coordX", "coordY"), crs="EPSG:4326")


## Make the plot
e <- ext(c(97.9079944169027, 98.0843053062245, 25.4895581045826, 25.7394596330535))
cols <- c("grey95", "violet")
r1 <- glad19
r1[r1==0] <- NA
v <- as.points(r1)

png("figures/figXX_monitorMap.png", res=350, 
    width=20, height=10, units="cm")
layout(matrix(1:2, nrow=1))
plot(glad19, legend=FALSE, main="GLAD Loss Year 2019", col=cols)
points(r1, col="violet", cex=0.05)
lines(as.polygons(e), lwd=0.5)
plot(crop(glad19, e), col=cols, 
     plg=list(legend=c("Unchanged", "Disturbed")))
add_legend(98.09, 25.72, legend=c("dLoc", "uLoc"), pch=16, 
           col=c("black", "blue"), bty="n", xpd=TRUE, cex=0.8)
points(v, cex=0.7)
points(u, col="blue", cex=0.7)
dev.off()


## Get numbers for caption

# change all 0s to NA
m <- c(0, NA)
rclmat <- matrix(m, ncol=2)
loss19 <- classify(dist19, rclmat)

# how many disturbances that we found are actually found by GLAD in 2019?
gv <- data.table(extract(loss19, v))
nrow(gv[!is.na(Layer_1)]) / nrow(gv)

# how many uLocs that we identified were labeled as disturbances by GLAD in 2019?
gu <- data.table(extract(loss19, u))
nrow(gu[!is.na(Layer_1)]) / nrow(gu)

