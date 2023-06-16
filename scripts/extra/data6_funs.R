##########################################################
## Purpose: Functions to run sen2cor software
##          Companion to data6_runSen2Cor...R scripts
## Creators: Ian McGregor, imcgreg@ncsu.edu
## Editors: Izzi Hinks, Xiaojie Gao
## System: R Version 3.6.3, July 2021, updated Nov 22
##########################################################

calcNDVI <- function(im, bands, specFolder){
  dirImg <- file.path(specFolder, im, "GRANULE")
  imgName <- list.files(dirImg)
  dirBands <- file.path(dirImg, imgName, "IMG_DATA/R20m")
  files <- list.files(dirBands, pattern=paste0(bands, collapse="|"),
                      full.names = TRUE)
  r1 <- rast(files[grepl("B04", files)])
  r2 <- rast(files[grepl("B8A", files)])
  r <- c(r1, r2)
    
  # calculate ndvi
  ## NIR - red / nir + red
  ndvi <- (r[[2]] - r[[1]]) / (r[[2]] + r[[1]])
  scl <- rast(files[3])
    
  ndviRast <- c(ndvi, scl)
  names(ndviRast)[1] <- gsub("B8A", "ndvi", names(ndviRast)[1])
    
  rastName <- gsub("SAFE", "tif", im)
  writeRaster(ndviRast, file.path(specFolder, rastName), overwrite=TRUE)

  ## fully delete L2A folder
  unlink(file.path(specFolder, im), recursive=TRUE)
}

# main function to run the Sen2Cor software
runSen2Cor <- function(command, txt_out_file, rmpi){
       tile_name <- strsplit(strsplit(command, "L2A/")[[1]][2], " --")[[1]][1]
       img_name <- strsplit(strsplit(command, paste0(tile_name, "/"))[[1]][2], "/")[[1]][1]

       if(rmpi){
        sprintf("Processing command %s for tile %s on rank %d on node %s with %d cores", command, tile_name, mpi.comm.rank(), Sys.info()[c("nodename")], parallel::detectCores())
       } else {
        sprintf("Processing command %s for tile %s on node %s with %d cores", command, tile_name, Sys.info()[c("nodename")], parallel::detectCores())
       }
       
       write(paste0(img_name, " started processing at ", Sys.time(), " with command: ", command), file=txt_out_file, append = TRUE)
       
       # RUN SEN2COR
       system(command) # uncomment to run the commands
       
       # Give output status (either "finished" or "error")
       specFolder <- strsplit(strsplit(command, "output_dir ")[[1]][2], " --resolution")[[1]][1]
       allFolders <- list.files(specFolder)
       imgl2a <- gsub("L1C", "L2A", img_name)
       imgNew <- gsub(paste0("_N.*"), "", imgl2a)
       
       latestImg <- data.table::last(allFolders[grepl(imgNew, allFolders)])
       fullPath <- paste0(specFolder, "/", latestImg, "/GRANULE")
       granFolder <- list.files(fullPath)
       nFolders <- length(list.files(paste0(fullPath, "/", granFolder, "/IMG_DATA")))
       nImgs10 <- length(list.files(paste0(fullPath, "/", granFolder, "/IMG_DATA/R10m")))
       nImgs20 <- length(list.files(paste0(fullPath, "/", granFolder, "/IMG_DATA/R20m")))
       
       if(nFolders<2 | nImgs10==0 | nImgs20==0){
         write(paste0("Error with ", img_name, ": no output images. Time of failure = ", Sys.time()), file=txt_out_file, append = TRUE)
       } else {
        ## calculate ndvi and save only NDVI + SCL bands
        calcNDVI(im=latestImg, bands=c("B04", "B8A", "SCL"), specFolder)
        write(paste0(img_name, " finished processing at ", Sys.time(), " with command: ", command), file=txt_out_file, append = TRUE)
       }
}

# function to set up parallelization of runSen2Cor
convertAllS2 <- function(nRun, rmpi){
  # Get vector of sen2cor commands (one command per tile)
  all_coms <- readRDS(paste0("s2Data/all_sen2cor_commands",
                             as.numeric(nRun), ".rds")) 
  txt_dir <- "s2Data/txt_files"
  txt_out <- file.path(txt_dir, paste0("sen2cor_out.txt"))
  
  # If txt_files dir doesn't exist, create it
  if (!dir.exists(txt_dir)) {
    dir.create(txt_dir)
  }
  
  # Create file for this specific execution of Sen2Cor
  if(!file.exists(txt_out)) file.create(txt_out)
     ## if want to do a couple test images (choose random images)
     # all_coms <- sample(all_coms, 10)
     
     ## if you want, run the function without actually converting the images to ensure
     ## multiple nodes are being used as directed (see `stdout` file)
     
     if(rmpi){
      cl <- makeCluster((mpi.universe.size()-1), type='MPI')
     } else {
      cl <- parallel::makeCluster(detectCores())
     }
     clusterEvalQ(cl, library(data.table))
     clusterEvalQ(cl, library(terra))
     clusterExport(cl, c("calcNDVI"), envir=environment())
     clusterOut <- clusterApply(cl=cl, x=all_coms, runSen2Cor, txt_out_file=txt_out, rmpi=rmpi)
     
     # Delete xml log file (path only has 1 GB memory)
     ## NOTE this is not strictly necessary now as we have changed the path of the log output from the L2A_process.py script itself.
     ## The new path of the log files is in the Sen2Cor home directory, which has as much storage as we need. The code here is being
     ## left for posterity.
     
     # homePath <- "/home/imcgreg/sen2cor/2.9/log"
     # if(length(list.files(homePath))>0) {
     #   deleteXML <- paste0("rm -r ", paste0(homePath, "/*")) #delete all files in the folder
     #   system(deleteXML)
     # }
     
     cat(unlist(clusterOut), sep='\n')
     snow::stopCluster(cl)

     if(rmpi){
      Rmpi::mpi.exit()
     }
}