# How to navigate this repo
This repository contains the scripts used for McGregor et al. 2023.

## Analyses
The step-by-step workflow for how the scripts should be run is in `workflow.md`. Below are some notes about each part of the analysis.

### Training data
As formatted, the training data analysis can be run as-is with no additional input because the sensor data is housed as separate csv files. If you want to inspect how the data was obtained in the first place, please see the "Sample points" and "Sentinel-2 imagery"" below.

#### Extra
Explanatory scripts for the initial identification and processing of the training data (as well as processing Sentinel-2 images) have been added here (e.g. `data1` - `data10`). The training data analysis in this repository can be run without these (note that Figure 10 in the paper is from `scripts/extra/data10...R`).

After the training data has been processed and analyzed, further visualization can be done by making diagnostic plots (Step 4 in workflow). Note this was not directly included in the paper.

### Landscape application (monitoring grid)
Due to the size of the monitoring grids and the memory needed for processing, all code for this analysis was run on the high performance computing (HPC) cluster at NC State. We are providing the code so it can be parsed if needed. 
- Technically, the code can be run as it is on a normal machine, but it will take a much longer time to finish. We are providing the HPC submission parameters used in hopes that they are useful (`scripts/hpcPars/`). The file names match the associated scripts. Note the file paths in these `.csh` files will need to be changed if you want to directly use them.

NOTE: for the landscape code to run as-is, virtual raster tiles (VRTs) need to be created from the sensor data (`scripts/app0B_prepareVRT.R`).
- The landscape application involves all images from Landsat-8 and Sentinel-2 for the timeframe of this study (June 2015 - January 2020). We have provided the cropped images already, but the script to download Landsat-8 images is included in the `scripts/extra/` folder (Sentinel-2 was the same manual process as for the training data).

## Obtaining data
### Sample points (Step 1 of workflow)
Training data (examples of deforestation) were identified from [PlanetScope](https://www.planet.com) imagery, verified in Google Earth Engine (GEE) and we recorded metadata (see `data/trainingPars/trainingDataPoints.csv`). We then used GEE to download sensor data for Landsat-8, Sentinel-2, and Sentinel-1.

### Sentinel-2 imagery
This research began in 2019, and training data was obtained for that year as a monitoring year. However, (as of Summer 2023) GEE does not have Sentinel-2 data for nothern Myanmar prior to December 2018 because it has not been released / processed by European Space Agency. Because of this, we were forced to manually download the Sen2Cor application and process imagery ourselves.
- It is  recommended that future users of the method in this repository avoid this path if avoidable, and try to use available L2A data on GEE or another platform. If unavoidable, we hope that providing the workflow here can help.