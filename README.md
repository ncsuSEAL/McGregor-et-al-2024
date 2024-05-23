# A multi-source change detection algorithm supporting user customization and near real-time deforestation detections
This repository contains the scripts and workflow used for [McGregor et al. 2024](https://www.sciencedirect.com/science/article/pii/S003442572400213X?via%3Dihub). I have verified that the initial analysis runs without issue, and the landscape analysis *should* run similarly, but that has not been tested (or altered) since ~spring 2023. Please feel free to open an issue if something pops up.

## Analysis
The step-by-step workflow for how the scripts should be run is in `workflow.md`. There are 2 steps to the analysis:
- The initial analysis with sample points (focal locations) covers the bulk of the paper and can be run via the script `train0_prepRunCode.R`.
- The landscape analysis is the extension of the sample points and refers to Section 2.3 in the paper. Please see notes about it below.

### Running the code
The scripts were originally developed using RStudio, but this repository was developed using R in a conda environment (see `deforMonitorEnv.yml`). The necessary packages with the correct versions are listed in `train0_prepRunCode.R` via the `groundhog` package.

## Notes
### Terminology
The entire analysis can be run as-is using the directory structure of this repository. The data for the focal locations (labeled as "trainingDataPoints") were originally downloaded from Google Earth Engine and are included here as csv files. If you want to inspect how the data was obtained in the first place, please see the "Sample points" and "Sentinel-2 imagery" sections below.
- Note that in the code there are several references to "training data". This is not training data in the sense of a machine learning application; rather, this was my way of differentiating variables at the start of the research, and the name stuck. While the paper originally used this terminology, we realized from reviews that it was confusing because of the discrepancy between machine learning applications and what we actually did in this study. The verbiage in the paper was changed to "focal locations", but the scripts still retain the original phrasing.

### Landscape application (monitoring grid)
Landscape regions were initially labeled as 1-5, but ultimately 2-4 were dropped from the final analysis. In the paper, Region 2 is the original Region 5; this is why the landscape code has the different numbers.

Due to the size of the monitoring grids and the memory needed for processing, all code for this analysis was run on the high performance computing (HPC) cluster at NC State. We are providing the code so it can be parsed if needed. 
- We are providing the HPC submission parameters used in hopes that they are useful (`scripts/hpcPars/`). The file names match the associated scripts. 
- Note the file paths in these `.csh` files will need to be changed if you want to directly use them.

For the landscape code to run as-is, virtual raster tiles (VRTs) need to be created from the sensor data (`scripts/app0B_prepareVRT.R`).
- The landscape application involves all images from Landsat-8 and Sentinel-2 for the timeframe of this study (June 2015 - January 2020). The script to download Landsat-8 images is included in the `scripts/extra/` folder (Sentinel-2 was the same manual process as for the training data).

### "Extra" code
The folder `scripts/extra` contains ancillary code to the primary analysis. This includes things like data identification, formatting, and necessary processing prior to running the main analysis. For the sake of this repository, these are considered extra because the sensor and focal location data is supplied here a priori.

## Obtaining data
### Sample points (Step 1 of workflow)
Examples of deforestation were identified from [PlanetScope](https://www.planet.com) imagery, verified in Google Earth Engine (GEE) and we recorded metadata (see `data/trainingPars/trainingDataPoints.csv`). We then used GEE to download sensor data for Landsat-8, Sentinel-2, and Sentinel-1.

### Sentinel-2 imagery
This research began in 2019, and training data was obtained for that year as a monitoring year. However, at the time of this paper's review and publication, GEE did not have Sentinel-2 data for northern Myanmar prior to December 2018 because it was not been released / processed by European Space Agency. Because of this, we had to manually download the Sen2Cor application and process the imagery to L2A ourselves. While we have provided the workflow we used, we recommend using some other repository of L2A data if possible (i.e. save yourself the trouble if you can)!
