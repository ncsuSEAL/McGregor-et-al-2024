################################################################################
## Purpose: Create the necessary commands for batch-running Sen2Cor (v2.9) in 
##          HPC (companion script = "S2_runSen2Cor_batch.R")
## Run medium: 
##  - HPC submitted job with accompanying .csh file: same name, in sen2cor_rmpi/
##  - PC if want, just add a `-1` to the `detectCores()` command below
## Creators: Ian McGregor, imcgreg@ncsu.edu
## Editors: Izzi Hinks, Xiaojie Gao
## System: R Version 3.6.3, July 2021
################################################################################
library(data.table)
library(parallel)

## -----------------------------------------------------------------------------
## Function = purpose
##
## containsImg = makes sure that all the necessary L1C directories / band images 
##              for L2A conversion are present
## lookS2 = create the explicit commands to run sen2cor for each tile.
## -----------------------------------------------------------------------------
containsImg <- function(f, dirName){
  sub <- list.files(file.path(dirName, f, "GRANULE"), full.names = TRUE)
  if (length(sub) == 0) return(FALSE)
  
  if(all(dir.exists(file.path(sub, "IMG_DATA")))){
    if(length(list.files(file.path(sub, "IMG_DATA"))) >= 7){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}
lookS2 <- function(w, outputDir=L2A, L1C_dir=L1C){
  tileDirHPC <- file.path(L1C_dir, w)
  L2A_outDirHPC <- file.path(outputDir, w)
  sen2cor_comPath <- file.path("Sen2Cor-02.09.00-Linux64/bin/L2A_Process")
  
  # Get all files in this tile's directory
  fileListOrig <- list.files(tileDirHPC)
  if (!is.null(regExp)) fileListOrig <- grep(regExp, fileListOrig, value = TRUE)
  
  # Add the tile directory within the L2A directory (if it doesn't already exist) 
  if (!dir.exists(L2A_outDirHPC)) {
    dir.create(L2A_outDirHPC)
  }
  
  ## when duplicated img name, only keep those that have larger file size 
  ### (indicates full data)
  nm <- sapply(fileListOrig, function(N){
    return(strsplit(N, "_2.{7}T.{6}.SAFE")[[1]][1])
  })
  nm <- unique(nm)
  
  fileList <- lapply(nm, function(G){
    files <- fileListOrig[grepl(G, fileListOrig)]
    
    if(length(files > 1)){
      sizeDups <- sapply(files, function(f){
        path <- list.files(file.path(tileDirHPC, f), all.files=TRUE, 
                           recursive=TRUE, full.names=TRUE)
        return(sum(file.info(path)$size))
      })
      return(files[which.max(sizeDups)])
    } else {
      return(files)
    }
  })
  
  # Create the commands to run Sen2Cor on each image in this tile's directory 
  eachImg <- lapply(fileList, function(imgPath){
    
    # If the necessary data for conversion are present
    if(containsImg(imgPath, tileDirHPC)) {
      datastrip_dir <- file.path(tileDirHPC, imgPath, "DATASTRIP")
      
      # Create the AUX_DATA dir (in 2 locations) for the image 
      ## (if it doesn't already exist)
      if (!dir.exists(file.path(tileDirHPC, imgPath, "AUX_DATA"))) {
        dir.create(file.path(tileDirHPC, imgPath, "AUX_DATA"))
      }
      
      subDir <- list.files(file.path(tileDirHPC, imgPath, "GRANULE"))
      if(!dir.exists(file.path(tileDirHPC, imgPath, "GRANULE", 
                               subDir, "AUX_DATA"))){
        dir.create(file.path(tileDirHPC, imgPath, "GRANULE", 
                             subDir, "AUX_DATA"))
      }
      
      if(dir.exists(datastrip_dir)){
        subfolder <- list.files(datastrip_dir)[1]
        # If the QI_DATA folder doesn't exist, this is an older L1C image; 
        # If get random errors in future sen2cor iteration, then 
        # possibly need to process it with Sen2Cor v2.5.5
        if (!dir.exists(file.path(datastrip_dir, subfolder, "QI_DATA"))) {
          sprintf("Image %s lacks a QI_DATA file", imgPath) 
          
          # See if adding the QI_DATA folder fixes the problem 
          if (!dir.exists(file.path(datastrip_dir, subfolder, "QI_DATA"))) {
            dir.create(file.path(datastrip_dir, subfolder, "QI_DATA"))
          }
        }
        
        # Create command to process the image with Sen2Cor v2.9
        ## note that from the user guide, if we do not explicitly include 
        ### `--resolution`, the default is to process the 20m bands first, 
        ### followed by the 10m band. This is the same workflow as when the 
        ### resolution is specified as 10m. Including `--resolution 10` is thus 
        ### a user preference.
        
        ## Also note that specifying the GIP paths is technically not necessary, 
        ### as earlier tests (admittedly not fully thorough) indicated the files
        ### were still used even when not specified (i.e. resultant images were 
        ### differenced with a result of 0), but they can be included for 
        ### peace of mind.
        img <- file.path(tileDirHPC, imgPath)
        addInsRes <- "--resolution 10"
        addInsGIP <- paste0("--GIP_L2A ", "Sen2Cor-02.09.00-Linux64/lib/python2.7/site-packages/sen2cor/cfg/L2A_GIPP.xml")
        addInsPB <- paste0("--GIP_L2A_PB ", "Sen2Cor-02.09.00-Linux64/lib/python2.7/site-packages/sen2cor/cfg/L2A_PB_GIPP.xml")
        addinDebug <- "--debug" #from the software guide, I think adding in debug means not running toolbox mode.
        com <- paste(sen2cor_comPath, "--output_dir", L2A_outDirHPC,
                     addInsRes, addInsGIP, addInsPB, img)
        
        imgAcq <- substr(imgPath, 12, 26) #image acquisition date
        prodID <- substr(imgPath, nchar(imgPath)-19, nchar(imgPath)-5) #prod ID
        return(data.table(imgAcq=imgAcq, prodID=prodID, com=com))
      } else{
        # If DATASTRIP dir doesn't exist for this image
        return(NA)
      }
    } else {
      return(NA) # No image
    }
  })
  
  noNA <- eachImg[!is.na(eachImg)]
  allImg <- as.data.table(rbindlist(noNA))
  
  ## format the data table
  allImg <- allImg[,`:=` (dateImg = as.POSIXct(imgAcq, format="%Y%m%dT%H%M%S"),
                          dateProd = as.POSIXct(prodID, format="%Y%m%dT%H%M%S"))
  ][order(dateProd), ]
  
  return(allImg[,com])
}

################################################################################
# Define necessary variables
L1C <- "s2Data/L1C"
L2A <- "s2Data/L2A"
tiles <- list.files(L1C)
## new tiles from Oct 2022
# tiles <- c("46QEK", "46QHK", "46QHM", "46RFN", "46RGR", "47QKD", "47QKE",
#             "47QKG", "47QLD", "47RKK", "47RKL", "47RLK")
rds_file <- "s2Data/all_sen2cor_commands"
regExp = "S2[A-Z]_MSIL1C_.*.SAFE$"

# path only needed for creating the command line commands
basePathHPC <- "/rsstu/users/j/jmgray2/SEAL/IanMcGregor"

cl <- makeCluster(detectCores())
clusterEvalQ(cl, library(data.table))
clusterExport(cl, list("cl", "tiles", "L1C", "L2A", "regExp", 
                       "containsImg", "lookS2"))

out <- parLapply(cl, X=tiles, lookS2, L2A, L1C)
stopCluster(cl)

## Data table output: list of data.tables (one for each tile), each with 
## columns: imgAcq, prodID, com, dateImg, dateProd
# f <- rbindlist(out)
# all_coms <- f[,com]

## Vector output: list of commands (one for each tile)
all_coms <- unlist(out)

## bring in current text output if available to filter
done <- read.delim("s2Data/txt_files/sen2cor_out.txt", header=FALSE)

test <- done[,1][grepl("finished", done[,1])]
vec <- lapply(test, function(X){
  b <- strsplit(X, ".SAFE")[[1]]
  bleh <- b[grepl("^S2.*", b)]
  return(bleh)
  })
vec <- unlist(vec)

vecFinished <- sapply(vec, function(X){
  out <- which(grepl(X, all_coms))
  return(out)
})
vecFinished <- as.numeric(vecFinished)

finished <- all_coms[vecFinished]
notFinished <- all_coms[!(all_coms %in% finished)]

# If you have >1000 images, then separate into chunks of equal size
## Lisa recommends chunks of <= 1000 commands per job.
## If you don't care about the amount of jobs you have, splitting up more is helpful
tt <- split(notFinished, ceiling(seq_along(notFinished)/25))

# Save each chunk of commands to an RDS file
sapply(1:length(tt), function(X){
  saveRDS(tt[[X]], file = paste0(rds_file, X, ".rds"))
})
