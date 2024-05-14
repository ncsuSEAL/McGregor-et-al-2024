# How to navigate this repo
This repository contains the scripts and workflow used for McGregor et al. 2024.

## Analyses
The step-by-step workflow for how the scripts should be run is in `workflow.md`. Below are some notes about each part of the analysis.

### Focal locations
The entire analysis can be run as-is using the directory structure of this repository. The data for the focal locations were originally downloaded from Google Earth Engine and are included here as csv files. If you want to inspect how the data was obtained in the first place, please see the "Sample points" and "Sentinel-2 imagery" sections below.

### Landscape application (monitoring grid)
Due to the size of the monitoring grids and the memory needed for processing, all code for this analysis was run on the high performance computing (HPC) cluster at NC State. We are providing the code so it can be parsed if needed. 
- We are providing the HPC submission parameters used in hopes that they are useful (`scripts/hpcPars/`). The file names match the associated scripts. 
- Note the file paths in these `.csh` files will need to be changed if you want to directly use them.

For the landscape code to run as-is, virtual raster tiles (VRTs) need to be created from the sensor data (`scripts/app0B_prepareVRT.R`).
- The landscape application involves all images from Landsat-8 and Sentinel-2 for the timeframe of this study (June 2015 - January 2020). The script to download Landsat-8 images is included in the `scripts/extra/` folder (Sentinel-2 was the same manual process as for the training data).

## Obtaining data
### Sample points (Step 1 of workflow)
Examples of deforestation were identified from [PlanetScope](https://www.planet.com) imagery, verified in Google Earth Engine (GEE) and we recorded metadata (see `data/trainingPars/trainingDataPoints.csv`). We then used GEE to download sensor data for Landsat-8, Sentinel-2, and Sentinel-1.

### Sentinel-2 imagery
This research began in 2019, and training data was obtained for that year as a monitoring year. However, at the time of this paper's review and publication, GEE did not have Sentinel-2 data for nothern Myanmar prior to December 2018 because it was not been released / processed by European Space Agency. Because of this, we had to manually download the Sen2Cor application and process the imagery to L2A ourselves. While we have provided the workflow we used, we recommend using some other repository of L2A data if possible (i.e. save yourself the trouble if you can)!