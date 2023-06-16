##########################################################
## Purpose: Download images from PLANET via creating geojsons of wanted points,
##          submitting orders to PLANET itself, then downloading the URLs once
##          the orders have successfully completed.
## Run medium:
##  - PC VScode terminal (faster than Mac vscode bc already connected to SEAL)
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, May 2022
## Last modified: July 2022
##########################################################

# Before running the code, we need to first activate conda, and we will open
## radian. Note that later on in the script we will need to get out of radian,
## but this is explained below.

## If you need assistance with this, please see the github page:
## https://github.com/tyson-swetnam/porder

# conda activate planet_orders
# radian

library(data.table)
library(terra)
library(parallel)

################################################################################
# Step 1a: Define variables
################################################################################
#Set the file paths and names of files/folders to be created 
# base_dir <- "/Volumes/SEAL/IanMcGregor/planetImagesNew"
base_dir <- "//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor/planetImagesNew"
csv_name <- "trainingDataPoints.csv"
shp_dir <- file.path(base_dir, "shapefiles")
geojson_dir <- file.path(base_dir, "geojsons") 
idlist_dir <- file.path(base_dir, "idlists")
orders_txt <- file.path(base_dir, "successful_orders.txt")
out_dir <- file.path(base_dir, "planetscope")

#Params for functions
bufferSize <- 100 #buffer (meters) of coordinate to make polygon to clip to
inputCRS <- "EPSG:4326" #crs of input lon/lat coords 
utmCRS <- "EPSG:32646" #desired output utm crs

#Params for planet API arguments
date_ordered <- "2022-07-13" # Date before you ordered
date_today <- "2022-07-15" # Date after you ordered the imagery
min_cloud_cover <- 0
max_cloud_cover <- 1 # number between 0.1-1 for % cloud cover
area_overlap <- 99 # ensure images overlap 100% of each buffered coordinate
days_to_consider <- 1 # N days before and after dist to look for images
num_to_select <- 1 # N images pre and post dist to order (e.g. N=2 =4 img total)

################################################################################
# Step 1b: Fixed variables
################################################################################
#Fixed params for planet api
max_imgs <- 500 # Max number of images to get per field
max_order <- 500 # Max number of assets per order for the department's license

#4-band imagery is required to get the pre-calculated NDVI band
##You can use "PSScene3Band" if you only want RGB bands 
item_type <- "PSScene4Band"

#NOTE: don't include UDM2 here if you also want images prior to UDM2 release
asset_type <- "analytic_sr,udm" 
bundle <- "analytic_sr_udm2,analytic_sr" # get SR and UDMs

# this extracts NDVI of the clipped images in a single zip file of all assets
# ops <- "clip zipall ndvi"

# as of Oct 2021 we have to download ndvi and the SR bands separately.
ops <- "clip zipall"

################################################################################
# Step 2: Define functions
################################################################################
# Convert lon/lat coordinates to shapefile using terra
coord2Shp <- function(pointNum, dt, bufferSize, shp_dir, inputCRS, utmCRS){
  
  # Create a shapefile (spatVector)
  subDT <- dt[pointNum, ]
  coords <- t(matrix(c(subDT$coordX, subDT$coordY)))
  centroid <- vect(coords, type="point", crs=inputCRS)
  
  # Convert to UTM so we can build the buffer in meters
  centroid_utm <- project(centroid, utmCRS)
  buff <- buffer(centroid_utm, width=bufferSize)
  
  # Convert to polygons
  buffPoly <- as.polygons(buff)
  
  # Write out buffered points as shapefile
  writeVector(buff, filename=paste0(shp_dir, "/", subDT$pointid, ".shp"),
              filetype="ESRI Shapefile", overwrite=TRUE)
  return(print(paste0("Finished with ", subDT$pointid)))
}

# Identify the PLANET images to order, and then actually order them
## NOTE: If you had more than 500 images per geojson, you would have to add a 
## command to split the idlist into separate .csvs of max length 500
getPSid <- function(pt, imgs_to_get, all_geojsons, idlist_dir, 
                    days_to_consider, item_type, asset_type, 
                    min_cloud_cover, max_cloud_cover, area_overlap,
                    num_to_select, bundle, ops){
  # Get the data for this coordinate
  coord_data <- imgs_to_get[pointid==pt, ]
  coord_file <- all_geojsons[grepl(paste0(pt, ".geojson"), all_geojsons)]
  
  # Make file paths for the csvs below
  pre_csv <- paste0(idlist_dir, "/", pt, "_pre.csv") # pre-dist IDs
  post_csv <- paste0(idlist_dir, "/", pt, "_post.csv") # post-dist IDs
  to_order_csv <- paste0(idlist_dir, "/", pt, "_toOrder.csv") # IDS to order
  
  # NOTE while the start date in `porder` is exclusive,
  ## the end date in `porder` is inclusive (!!!)
  # Identify the N (days_to_consider) images before the disturbance, and record
  ## the IDS in a csv
  pre_disturb_start <- coord_data$datePre - days_to_consider
  # pre_disturb_end <- as.Date("2019-09-17")
  system(paste0("porder idlist --input ", coord_file, " --start ", 
                pre_disturb_start, " --end ", pre_disturb_end, " --item ", 
                item_type, " --asset ", asset_type, " --outfile ", pre_csv, 
                " --cmin ", min_cloud_cover, " --cmax ", max_cloud_cover, 
                " --overlap ", area_overlap))
  system(paste0("porder idlist --input ", coord_file, " --start ", 
                pre_disturb_start, " --end ", (coord_data$datePre), " --item ", 
                item_type, " --asset analytic_sr", " --outfile ", pre_csv))
  # system(paste0(porder --asset analytic_sr --outfile <bla>))
  
  # Identify the N (days_to_consider) images after the disturbance, and record
  ## the IDS in a csv
  post_disturb_end <- coord_data$dateDist + days_to_consider -1
  system(paste0("porder idlist --input ", coord_file, " --start ", 
                (coord_data$dateDist - 1), " --end ", post_disturb_end, " --item ", 
                item_type, " --asset ", asset_type, " --outfile ", post_csv, 
                " --cmin ", min_cloud_cover, " --cmax ", max_cloud_cover, 
                " --overlap ", area_overlap))
  
  # Read in the csv files from above
  if(file.info(pre_csv)$size > 0){
    pre_dist_IDs <- fread(pre_csv, header=FALSE)
  } else {
    return(print(paste0("No Pre images detected for ", pt)))
  }

  if(file.info(post_csv)$size > 0){
    post_dist_IDs <- fread(post_csv, header=FALSE)
  } else {
    return(print(paste0("No Post images detected for ", pt)))
  }
  
  # Identify only the IDs of the most recent N (num_to_select) images both
  ## before and after the disturbance. E.g. if N = 2, then retain 4 IDS in total
  selected_IDs <- c(head(pre_dist_IDs[[1]], n = num_to_select),
                    tail(post_dist_IDs[[1]], n = num_to_select))
  
  # Write selected IDs to a new .csv file 
  fwrite(as.data.table(selected_IDs), file=to_order_csv, col.names=FALSE)
  
  # Order them from PLANET
  system(paste0("porder order --name ", pt, " --idlist ", to_order_csv, 
                " --item ", item_type, " --bundle ", bundle, " --boundary ", 
                coord_file, " --op ", ops))
}

# Download the orders to zip folders
downloadOrders <- function(date_ordered, date_today, orders_txt, out_dir){
  # # Get the successful orders and save their summaries to a .txt file
  # system(paste0("porder ostate --state success --start ", date_ordered, 
  #               " --end ", date_today, " >> ", orders_txt))
  
  # Then, extract the URLs of successful orders from the text file and use the 
  ## API to download the orders from the URLs 
  summary_txt <- read.delim(orders_txt, sep="|", header=FALSE)
  
  # Use Planet's Orders API to download the orders with the URLs
  urls <- lapply(summary_txt, function(orderData){
    return(sub(".*(https:\\S+).*", "\\1", orderData))
  })
  urls_to_order <- urls[["V4"]]
  urls_to_order <- urls_to_order[4:(length(urls_to_order)-1)]
  
  # Order the PS images for all of the data! 
  # Because we're running this on a login node, I will keep this serial rather than parallel
  sapply(urls_to_order, function(url){
    system(paste0("porder multipart --url ", url, " --local ", out_dir))
  })
}

# Wrapper for unzipPlanet and formats the ids it's going over
extractImgs <- function(out_dir){
  zips <- list.files(out_dir, pattern="*.zip")
  
  # Extract point ID to create non-zip folder
  id <- unique(as.vector(sapply(zips, function(W){
    strsplit(W, "_[[:alnum:]]{4}")[[1]][1]
  })))
  
  sapply(id, unzipPlanet, zips, out_dir)
}

# Unzip the downloaded folders and put the constituent files in folders labeled
## by the point id
unzipPlanet <- function(Q, zips, out_dir){
  files <- zips[grepl(paste0(Q, "_"), zips)]
  
  if(length(files)>1){
    test1 <- file.info(paste0(out_dir, "/", files[1]))
    test2 <- file.info(paste0(out_dir, "/", files[2]))
    
    if(test1$ctime < test2$ctime){
      fold1 <- "ndvi"; fold2 <- "SR"
    } else {
      fold1 <- "SR"; fold2 <- "ndvi"
    }
    
    dirOut1 <- paste0(out_dir, "/", Q, "_", fold1) #create separate folder
    dirOut2 <- paste0(out_dir, "/", Q, "_", fold2) #create separate folder
    
    unzip(paste0(out_dir, "/", files[1]), exdir=dirOut1, overwrite=TRUE)
    unzip(paste0(out_dir, "/", files[2]), exdir=dirOut2, overwrite=TRUE)
  } else {
    dirOut <- paste0(out_dir, "/", Q) #create separate folder
    unzip(paste0(out_dir, "/", files), exdir=dirOut, overwrite=TRUE)
  }
}

################################################################################
# Step 3: Run the code
################################################################################
# Import data and create folders
## Import the .csv of data
imgs_to_get <- fread(file.path(base_dir, csv_name))

## Ensure that the dates are in date format
imgs_to_get$dateDist <- as.Date(imgs_to_get$dateDist, format="%d/%m/%Y")
imgs_to_get$datePre <- as.Date(imgs_to_get$datePre, format="%d/%m/%Y")

## Filter out undisturbed points
imgs_to_get <- imgs_to_get[!grepl("a", pointid), ]
pointids <- unique(imgs_to_get$pointid)

# If the output folders don't not yet exist, create them
if(!dir.exists(geojson_dir)) dir.create(geojson_dir)
if(!dir.exists(shp_dir)) dir.create(shp_dir)
if(!dir.exists(idlist_dir)) dir.create(idlist_dir)
if(!dir.exists(out_dir)) dir.create(out_dir)
if(!file.exists(orders_txt)) file.create(orders_txt)

# Create shapefiles of polygons to be used for planet ordering
## NB: If you changed the buffer size, MUST re-run this and geojson code!!!
cl <- makeCluster(detectCores()-1)
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(terra))
clusterExport(cl, c("bufferSize", "imgs_to_get", "shp_dir", "coord2Shp"))
parallel::parSapply(1:nrow(imgs_to_get), coord2Shp, dt=imgs_to_get, bufferSize, 
                    shp_dir, inputCRS, utmCRS, cl=cl)
stopCluster(cl)

# Use PLANET API to identify images and order them.
## Make sure you're already in your conda environment with porder installed
### If you haven't already, please go read the notes on first-time use of this 
### script here: 
### https://github.com/ncsuSEAL/sealR/blob/master/downloadingPlanet/firstTime.md

system("planet init")

# Convert all .shp files to .geojson
system(paste0("porder convert --source ", shp_dir, " --destination ", geojson_dir))

# Now list all the geojsons and order them
all_geojsons <- list.files(path=geojson_dir, pattern="*.geojson", 
                           full.names=TRUE)
psIDs <- lapply(pointids, getPSid, imgs_to_get, all_geojsons, idlist_dir,
                    days_to_consider, item_type, asset_type, 
                    min_cloud_cover, max_cloud_cover, area_overlap,
                    num_to_select, bundle, ops)

################################################################################
# Step 4: Post-ordering
################################################################################
# To check on the progress, uncomment the line with the state of interest 
## (state = queued, running, success, failed, partial) 
# system(paste0("porder ostate --state state --start ", date_ordered, 
# " --end ", date_today))

# Create a text file of the URLs to download
## NOTE that this does not work in `radian` as of July 2022 due to the ">>"
## In order to run this, need to open a new terminal in VScode, activate conda,
## initialize planet, then copy paste the result of this.
### If the issue in radian is fixed, uncomment these lines within the 
### downloadOrders function.
paste0("porder ostate --state success --start ", date_ordered, 
                " --end ", date_today, " >> ", orders_txt)

# Download the orders to zip files
## NOTE this does not work in radian as of July 2022. Need to open new terminal,
## activate conda again, initialize planet again, then run this
downloadOrders(date_ordered, date_today, orders_txt, out_dir)

# Unzip the files and put into their own folders
extractImgs(out_dir)


bleh <- rbindlist(lapply(completed, function(X){
  out <- ifelse(file.info(X)$size > 0, "good", "missing")
  pt <- gsub("//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor/planetImagesNew/idlists/", 
          "", X)
  return(data.table(pt=pt, status=out))
  }))
bleh <- bleh[status == "missing", ][order(pt, )]
fwrite(bleh, "missingPlanet.csv")
