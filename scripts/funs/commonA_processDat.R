##########################################################
## Purpose: Process SR data from Landsat-8, Sentinel-2, and Sentinel-1
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 4.1.1, Nov 2021
## Last modified: Jan 2022
##########################################################

## -----------------------------------------------------------------------------
## Function name and purpose:
##
## parseL8SR_pixel = parse the QA_PIXEL bit from Landsat 8 SR
## parseL8SR_radsat = parse the QA_RADSAT bit from Landsat 8 SR
## parseL8SR_aerosol = parse the SR_QA_AEROSOL bit from Landsat 8 SR
## qa_remap = label sensor obs as good or bad depending on QA flags
## calcIndex = calculate NDVI for S2 and L8 **assuming the bandNames are in 
##              order (e.g., c(Band 4, Band 5)**
## filterThenIndex = screen out obs based on output from qa_remap and invalid values
## -----------------------------------------------------------------------------
parseL8SR_pixel <- function(x){
  # QA_PIXEL: We only want to keep the best quality images, so only those that are
  ## clear and have low confidences
  
  # Binary
  ## Bit 0 - if pixel is fill, then true
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  ## Bit 1 - if dilated cloud, then true
  dilatedCloud <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  ## Bit 2 - if cirrus, then true
  cirrus <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  ## Bit 3 - if cloud, then true
  cloud <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  ## Bit 4 - if cloud shadow, then true
  cloudShadow <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  ## Bit 5 - if snow, then true
  snow <- ifelse(bitwAnd(bitwShiftR(x,5),1), TRUE, FALSE)
  ## Bit 6 - if clear, then true
  clear <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  ## Bit 7 - if water, then true
  water <- ifelse(bitwAnd(bitwShiftR(x,7),1), TRUE, FALSE)
  
  # Confidences
  ## Confidences should be interpreted as the answer to the question: "What are 
  ### the chances I will see X outside?", with X being cloud, cloud shadow, etc. 
  
  ## Bits 8-9 - if cloud conf low or no level set, then false
  ### 0=no level set, 1=low, 2=medium, 3=high
  cloudConf <- ifelse(bitwAnd(bitwShiftR(x,8), 3) == 1, FALSE, TRUE)
  ## Bits 10-11 - if cloud shadow confidence low, then false
  ### 0=no level set, 1=low, 2=reserved, 3=high
  cloudShadowConf <- ifelse(bitwAnd(bitwShiftR(x,10), 3) == 1, FALSE, TRUE)
  ## Bits 12-13 - if snow/ice confidence low, then false
  ### 0=no level set, 1=low, 2=reserved, 3=high
  snowConf <- ifelse(bitwAnd(bitwShiftR(x, 12), 3) == 1, FALSE, TRUE)
  ## Bits 14-15 - if low cirrus confidence, then false; 
  ### 0=no level set, 1=low, 2=reserved, 3=high
  cirrusConf <- ifelse(bitwAnd(bitwShiftR(x,14), 3) == 1, FALSE, TRUE)
  
  return(list(fill=fill, dilatedCloud=dilatedCloud, cirrus=cirrus, cloud=cloud, 
              cloudShadow=cloudShadow, snow=snow, clear=clear, water=water,
              cloudConf=cloudConf, cloudShadowConf=cloudShadowConf, 
              snowConf=snowConf, cirrusConf=cirrusConf)
  )
}
parseL8SR_radsat <- function(x){
  # QA_RADSAT: We only want best images, so no saturation and no occlusion
  
  #is saturated?
  b1 <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  b2 <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  b3 <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  b4 <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  b5 <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  b6 <- ifelse(bitwAnd(bitwShiftR(x,5), 1), TRUE, FALSE)
  b7 <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  #band 8 is not used
  b9 <- ifelse(bitwAnd(bitwShiftR(x,8), 1), TRUE, FALSE)
  terrainOcclusion <- ifelse(bitwAnd(bitwShiftR(x,11), 1), TRUE, FALSE)
  
  return(list(
    b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, b6=b6, b7=b7, b9=b9,
    terrainOcclusion=terrainOcclusion
  ))
}
parseL8SR_aerosol <- function(x){
  # SR_QA_AEROSOL: We want best images, so no fill, no water, and low aerosol
  ## difference with climatology if correction applied (see user guide link in
  ## L8 section below)
  
  # Bit 0: if fill, then true
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  # Bit 2: if water, then true
  water <- ifelse(bitwAnd(bitwShiftR(x, 2), 1), TRUE, FALSE)
  #aerosol level; 0=climatology (no correction), 1=low, 2=med, 3=high)
  aerosolLow <- ifelse(bitwAnd(bitwShiftR(x, 6), 3) < 2, TRUE, FALSE)
  return(list(fill=fill, water=water, aerosolLow=aerosolLow))
}
qa_remap <- function(G, qa=""){
  if(qa=="scl"){
    # https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-2-msi/level-2a/algorithm
    return(ifelse(G %in% c(4,5), "good", "bad"))
  }
  
  # https://www.usgs.gov/media/files/landsat-8-9-collection-2-level-2-science-product-guide
  if(qa=="pixel"){
    bit_output <- parseL8SR_pixel(G)
    cond_false <- c("fill", "dilatedCloud", "cirrus", "cloud", "cloudShadow", 
                    "snow", "water", "cloudConf", "cloudShadowConf", 
                    "snowConf", "cirrusConf")
    cond_true <- c("clear")
    
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0 &
                    sum(bit_output[[cond_true]]) == 1, "good", "bad"))
  }
  
  if(qa=="radsat"){
    bit_output <- parseL8SR_radsat(G)
    cond_false <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7",
                    "b9", "terrainOcclusion")
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0, "good", "bad"))
  }
  
  if(qa=="aerosol"){
    ### pixel = 0 && radsat = 0. We are not filtering out based on aerosol due to low numbers of obs
    bit_output <- parseL8SR_aerosol(G)
    cond_false <- c("fill", "water")
    cond_true <- c("aerosolLow")
    #in other words, the data is valid (TRUE), if aerosol is low.
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0 &
                    sum(bit_output[[cond_true]]) == 1, "good", "bad"))
  }
}
calcIndex <- function(input, s=sat, b=bandNames, vegIndex){
  
  # no scaling for S2 bc already done in L2A conversion
  bandA <- input[[b[1]]]
  bandB <- input[[b[2]]]
  
  # scaling comes from Table 6-1 https://www.usgs.gov/media/files/landsat-8-9-collection-2-level-2-science-product-guide
  if(grepl("landsat", s)){
    bandA <- bandA*0.0000275 + (-0.2)
    bandB <- bandB*0.0000275 + (-0.2)
  }
  
  if(vegIndex=="ndvi") index <- (bandB - bandA) / (bandB + bandA)
  if(vegIndex=="evi2") index <- 2.5*((bandB - bandA) / (bandB + 2.4*bandA + 1))
  return(index)
}
filterThenIndex <- function(input, sat, bandNames, vegIndex){
  
  if(grepl("landsat", sat)){
    qaBands <- bandNames[!grepl("B", bandNames)]
    indexBands <- bandNames[grepl("B", bandNames)]
    
    qaBands <- qaBands[!grepl("AEROSOL", qaBands)]
    qaNames <- gsub("QA_", "", qaBands)
    
    #filter out values that are not valid (outside of valid range based on 
    ##user guide - see above in calcIndex)
    for(j in indexBands){
      input[[j]][input[[j]] < 7273 | input[[j]] > 43636] <- NA
    }
    
    ## QA PIXEL valid range
    input$QA_PIXEL[input$QA_PIXEL < 21824 | input$QA_PIXEL > 65534] <- NA
    
    ## QA RADSAT valid range
    input$QA_RADSAT[input$QA_RADSAT < 0 | input$QA_RADSAT > 3829] <- NA
  } else {
    indexBands <- bandNames[1]
    qaBands <- bandNames[2]
    qaNames <- bandNames[2]
  }
  
  for(i in 1:length(qaBands)){
    look <- data.table(v=sort(unique(as.numeric(input[[qaBands[i]]]))))
    look[, qa := sapply(v, qa_remap, qa=tolower(qaNames[i]))] # apply qa_map function
    
    for(j in indexBands){
      #screen out qa values that are invalid
      input[[j]][input[[qaBands[i]]] %in% look[qa == "bad", v]] <- NA 
    }
  }
  
  if(grepl("sentinel2", sat)){
    return(input$NDVI)
  } else {
    return(calcIndex(input, s=sat, b=bandNames, vegIndex))
  }
}
