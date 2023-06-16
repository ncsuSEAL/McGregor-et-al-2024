##########################################################
## Purpose: Extract the results from a .txt file created 
##          in the execution of S2_runSen2Cor_batch.R 
##          (input file: sen2cor_out_<date>_<time>.txt)
## Run medium
##  - PC (could do Mac, but PC already connected to SEAL so is faster)
## Creators: Izzi Hinks, irhinks@ncsu.edu
## Editors: Ian McGregor
## System: R Version 4.0.2, July 2021
##########################################################
library(stringr)

# Set base path
# setwd("Z:/IanMcGregor/")

###################################################################
## full missing images
### NOTE if sen2cor has passed over an img, then a folder IS created, but the
### actual processing may not have occurred. Need to see error messages.
newTiles <- c("46QEK", "47QLD", "47QKE", "46RFN", "47QKG", "47RKK", "46RGR", 
              "47QKD", "47RKL", "46QHK", "47RLK", "46QHM")
miss <- lapply(newTiles, function(X){
  l1c <- list.files(file.path("s2Data/L1C", newTiles[1]))
  l1c <- gsub("_N.*", "", l1c)
  
  l2a <- list.files(file.path("s2Data/L2A", newTiles[1]))
  l2a <- gsub("_N.*", "", l1c)
  
  miss <- setdiff(l1c, l2a)
  return(miss)
})
names(missImg) <- newTiles

###################################################################
txt_file <- read.delim("s2Data/txt_files/sen2cor_out.txt", header=FALSE)

# Gather data about the batch conversions 
# Num of finished conversions
finished_conversions <- sum(str_count(txt_file$V1, pattern = "finished"))

finished_strings <- txt_file[,1][grepl("finished", txt_file[,1])]
error_strings <- txt_file[,1][grepl("Error", txt_file[,1])]

# Date and time of first conversion
first_finish <- as.POSIXct(strsplit(strsplit(finished_strings[1], split = "processing at ")[[1]][2], split = " with command")[[1]][1], format="%Y-%m-%d %H:%M:%S")  
# Date and time of last conversion
last_finish <- as.POSIXct(strsplit(strsplit(finished_strings[finished_conversions], split = "processing at ")[[1]][2], split = " with command")[[1]][1], format="%Y-%m-%d %H:%M:%S")  
# Time taken to process the images (thus far, if the job is still running)
time_diff <- round(as.numeric(difftime(last_finish, first_finish, units="hours")),2)

# Calculate average processing times 
imgs_per_hour <- round(finished_conversions / time_diff, 2)
mins_per_img <- round(time_diff / finished_conversions * 60, 2)

# Print results
sprintf("%d images have been converted in %s hours; the first image finished at %s and the last finished at %s", finished_conversions, time_diff, first_finish, last_finish)
sprintf("That's an average of %f images processed per hour, with each image taking approximately %f minutes!", imgs_per_hour, mins_per_img)

# Edits to come on July 26: Izzi will add a piece about the number of L2A images that have empty folders (e.g., AUX_DATA, IMG_DATA, etc.) or other issues

### make sample images to test HPC warning
base <- "//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor/s2Data" #for windows

imgs <- readRDS(paste0(base, "/all_sen2cor_commands5.rds"))
imgs <- imgs[1:30]
imgs <- gsub("/L2A/", "/L2AtestHPC/", imgs)
saveRDS(imgs, file=paste0(base, "/testImgsHPC.rds"))


##########################################################################################
# Run above code in loop and extract every image that had an error.
## Then get the commands of those only and write a new Rds file to run them through
## sen2cor again (after adding an AUX_DATA folder)

base <- "rsstu/users/j/jmgray2/SEAL/IanMcGregor/s2Data" # if in SEAL dir #for mac
base <- "//oitrspprd.hpc.ncsu.edu/rsstu/users/j/jmgray2/SEAL/IanMcGregor/s2Data" #for windows

files <- list.files(file.path(base, "txt_files"))
files <- files[!grepl("prev", files)]

err <- lapply(1:length(files), function(X){
  txt_file <- read.delim(file.path(base, "txt_files", files[X]), header=FALSE)
  error_strings <- txt_file[,1][grepl("Error", txt_file[,1])]
  
  imgs <- as.vector(sapply(error_strings, function(Y){
    return(strsplit(strsplit(Y, "with ")[[1]][2], " due")[[1]][1])
  }))
  
  return(imgs)
})

## now add in missing AUX_DATA folders
allErr <- as.vector(unlist(err))

## now get commands for these images and make new rds file
commFiles <- list.files(base, pattern="commands")
allErr <- paste0(allErr, collapse="|")
redoComm <- sapply(1:length(commFiles), function(X){
  commFile <- readRDS(file.path(base, commFiles[X]))
  b <- commFile[grepl(allErr, commFile)]
  return(b)
})

out <- as.vector(unlist(redoComm))
saveRDS(out, file=paste0(base, "/sen2cor_commands_redo1.rds"))

