# Workflow steps with corresponding script names and locations

## Description
The steps are listed in the order that the scripts should be run. 

Medium = Mac, PC, Either, or HPC
- "Either" refers to either Mac or PC

Figure # | Script
---|---|
1 | 
2 | draw.io
3 | `data10_altApps.R`
4 | `train6_accuracy.R`
5 | `train7_plotAccuracy.R`
6 | `train7_plotAccuracy.R`
7 | `train7_plotAccuracy.R`
8 | `SEAL/app3_dailyProbMaps.R`
9 | `paper1_metadata.R`
10 | `train9_plotDiagnostics.R`
S1 | `train9_plotDiagnostics.R`
S2 | equations
S3 | `train2_sdInflation.R`
S4 | `train9_plotDiagnostics.R`
S5 | `train7_plotAccuracy.R`
S6 | `train7_plotAccuracy.R`
S7 | `train7_plotAccuracy.R`
S8 | `train7_plotAccuracy.R`
S9 | `train8_crossVal.R`

## Training data analysis
Step #|Task|Medium|Script|Output
---|---|---|---|---|
||**Step 3: Training Exploratory: Compare different sensor combinations for training data**|
1|Pre-process training data, compute ts models, calculate residuals, aggregate to 1 ts | Either | `scripts/train1_getResiduals.R` | Rdata files
2|Calculate seasonality adjustment based on raw z-scores | Either | `scripts/train2_sdInflation.R` | vector in Rdata file
3|Re-run step 1a now using the seasonal adjustment. | Either | `scripts/train1_getResiduals.R` | Rdata files
4|Calculate monitoring discount factor to use in 1e and landscape. | Either | `scripts/train2_sdInflation.R` | number; update in `args.R`
5|Define probability functions for each lambda per sensor combination | Either | `scripts/train3_createProbFuns.R` | .Rdata files
6|Calculate ewma, logMod, and probs for a vector of possible lambdas | Either | `scripts/train4_runLambdas.R` | csv files (separate for each lambda)
7|Prep accuracy metrics for using a sensitivity analysis over all thresholds| Either | `scripts/train5_metricsThresholds.R` | csv file
8|Calculate overall accuracy metrics and plot comparison panels| Either | `scripts/train6_accuracy.R` | csv files
9|Choose a threshold based on output from Step 8 and re-run Step 7 using that singular threshold. | | |
10|Plot results from second run of Step 7 to visualize F1, PR, and lag differences btwn lambdas and sensors | Either | `scripts/train7_plotAccuracy.R` | plots
11|Choose the best lambda and sensor combination for landscape application based on output of Step 10|
12|Run k-fold cross-validation for chosen sensor combo and chosen lambda. | Either | `scripts/train8_crossVal.R` | plots


## Landscape application
NOTE - Step 5-1 and 5-2 don't need to be run; the images are already in this repository.

Step #|Task|Medium|Script|Output
---|---|---|---|---|
||**Step 5: Application: Obtain data**|
1|If not already done from training data, download S2 images for Chatthin and convert using sen2cor (see above)|
3|Download L8 images to drive folder | GEE | `downloadL8App.js` |L8 images
5|Create text files that are needed for gdal to build the VRTs | Either | `scripts/app0B_prepareVRT.R` | text files
6|Create VRTs for each sensor's images | PC | `scripts/app0B_prepareVRT.R` | VRTs
| | | | |
||**Step 6: Application: Process data**|
0|Make sure all necessary data has been copied over to `SEAL/dissertation/myanmar/trainingPars/`, especially train1, train2, and train3 outputs.
1|Read in VRT data and create matrices for all of Chatthin | HPC job| `scripts/app1_createNDVImat.R`; `SEAL/Ian/app1_createNDVImat.csh` |matrix binary files
2|Process the landscape data and get prob ts using same functions as training data | HPC | `scripts/app2_processLandProbs.R`; `SEAL/Ian/app2_processLandProbs.csh` | prob ts files
3|Create daily landscape maps as separate pngs | HPC interactive session (fastest), or Mac/PC. Submitted job doesn't work for some reason | `scripts/app3_dailyProbMaps.R`; `SEAL/Ian/app3_dailyProbMaps.csh` with `bayes=FALSE; binaryDist=FALSE` | png maps
4|Convert pngs into gif | Either, but not HPC| `scripts/app3_dailyProbMaps.R`; `SEAL/Ian/app3_dailyProbMaps.csh` with `bayes=FALSE; binaryDist=FALSE` | gif
5|Calculate ratio of disturbed pixels per day | HPC job | `scripts/app3_dailyProbMaps.R`; `SEAL/Ian/hpcApp3_dailyProbMaps.csh` with `bayes=FALSE; binaryDist=FALSE` | .Rdata file (vector)
| | | | |
||**Step 7: Application: Analyze results**|
1|Analyze results of landscape application by spotchecking pixels | Either | `scripts/app4_landscapeSpotCheck.R`, `funs/commonE_plotSummaries.R`| none
2|Plot the daily dist ratio of all regions together | Either | `scripts/app4_landscapeSpotCheck.R` | plots
3|Create table of validation pixels and dates to look over | PC | `scripts/app5_validation.R` | csv file
4|Record validation metrics from looking through planet imagery | Manual | no script | updated csv file
5|Calculate validation accuracy metrics | Either | `scripts/app5_validation.R` | metrics


## Extras
These scripts are not necessary to run the analysis using the data already included in the repository

### Obtain data for training data analysis
Step #|Task|Medium|Script|Output
---|---|---|---|---|
||**Step 1: Identify and download training data**|
1|Identify training data via PLANET and GEE for verification |GEE| `GEE/data0_verifyTrainingData.js` | data
2|Format training data from google forms| Either | `scripts/data1_formatTrainingData.R` | data
2a|Run Step 2-2 *now* to see if you are filtering out any points. This way you don't have to waste time getting more data from GEE. | | |
3|Download satellite data by uploading the csv to GEE |GEE| `GEE/data2_getImagery.js` | csv files with band data per training point for the sensors
| | | | |
||**Step 2: Process S2 using Sen2Cor**|
1|Don't need to do this if GEE has full S2 ts. Please read `github/hpcHelp/Sen2Cor_HPC/howToGuide.md` and `hpc/IanMcGregor/S2_0workflow.R` for context|
2|Determine which S2 tiles to use | Either | `scripts/data3_s2TilesL1C.R` | csv of tiles by point (`trainingDataS2Tiles.csv`)
3|Download S2 images | PC VScode | `scripts/data3_s2TilesL1C.R` | S2 images for training data
4|Merge S2 DEMs | PC VScode | `scripts/data4_mergeS2DEMs.ipynb` | large DEM over northern Myanmar
5|Create the terminal commands for running sen2cor | PC VScode | `SEAL/Ian/sen2cor_rmpi/data5_createS2Orders.R`; `data5_createS2Orders.csh` | text file with terminal commands
6|Batch run sen2cor, convert images from L1C to L2A, calculate NDVI, and delete L1C images to save memory. There is both a manual and auto version depending on how want to submit HPC jobs| HPC job | `SEAL/Ian/sen2cor_rmpi/data6_runSen2Cor_batch*.R`; `data6_runSen2Cor_batch*.csh`|L2A images
7|Check progress of sen2cor |PC| `SEAL/Ian/sen2cor_rmpi/data7_trackSen2CorProgress.R` | analysis
7b|Calc NDVI and delete L1C if this didn't happen within the sen2cor script | HPC | `SEAL/Ian/scripts/ch1_ndviL2A.R`; `SEAL/Ian/ch1_ndviL2A.csh` | L2A tifs
8|Extract data from L2A images and put in same format as downloaded GEE data | HPC job | `SEAL/Ian/data8_extractS2Bands.csh`; `scripts/data8_extractS2Bands.R`; `scripts/funs/dataA_extractS2.R` |csv
9|Identify forest strata covered by training data | PC | `scripts/data9_forestType.R` + `funs/dataB_maskForest.R` | cropped tif of strata
10|Identify other regions for landscape application | Either | `scripts/data10_altApps.R` + `funs/dataC_idRegions.R` | Rdata file of extents

### Further visualization of training data
Step #|Task|Medium|Script|Output
---|---|---|---|---|
||**Step 4: Training: Create summary plots with chosen lambda / sensor combo**|
0|Use best lambda and sensor combination to make logistic model for landscape application (non-backfilled data) | Either | `scripts/train6_createLogMod.R` | .Rdata file
0|^Note as of March 21 2023 we have abandoned log mod method |||
2|Download planet images for diagnostic plots | PC (faster) on VScode | `scripts/train9_getPlanetImagery.R` | tif files
3|Make diagnostic plots (8) for training data, including PLANET before/after images | Either | `train10_plotDiagnostics.R` | plots

### OBtain data for landscape application analysis
See Step 5 above in landscape application
