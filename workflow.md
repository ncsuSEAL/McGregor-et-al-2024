# Workflow steps with corresponding script names and locations

## Description
This is the same workflow as in the flowchart I've created but the associated scripts are included along with their general location. Each section's steps are in order.

Medium = Mac, PC, Either, or HPC
- "Either" refers to either Mac or PC

Figure # | Script
---|---|
1 | `data10_altApps.R`
2 | draw.io
3 | `data10_altApps.R`
4 | `train6_accuracy.R`
5 | `train7_plotOthers.R`
6 | `train7_plotAccuracy.R`
7 | `paper1_metadata.R`
8 | `SEAL/app3_dailyProbMaps.R`
9 | `train7_plotAccuracy.R`
S1 | `train10_plotDiagnostics.R`
S2 | `train10_plotDiagnostics.R`
S3 | `train2_sdInflation.R`
S4 | `train10_plotDiagnostics.R`
S5 | `train7_plotAccuracy.R`
S6 | `train7_plotAccuracy.R`
S7 | `train7_plotOthers.R`
S8 | `train8_crossVal.R`
S9 | `train10_plotDiagnostics.R`
S10| `monitoringFigureCh1.R`

Step #|Task|Medium|Script|Output
---|---|---|---|---|
||**Step 1: Identify and download training data**|
1|Identify training data via PLANET and GEE for verification |GEE| `GEE/data0_verifyTrainingData.js` | data
2|Format training data from google forms| Either | `scripts/data1_formatTrainingData.R` | data
2a|Run Step 2-2 *now* to see if you are filtering out any points. This way you don't have to waste time getting more data from GEE. | | |
3|Download satellite data by uploading the csv to GEE |GEE| `GEE/data2_getImagery.js` | csv files with band data per training point for the sensors
| | | | |
||**Step 2: Process S2 using Sen2Cor**|
1|Don't need to do this if GEE has full S2 ts|
2|Determine which S2 tiles to use | Either | `scripts/data3_s2TilesL1C.R` | csv of tiles by point (`trainingDataS2Tiles.csv`)
3|Download S2 images | PC VScode | `scripts/data3_s2TilesL1C.R` | S2 images for training data
4|Merge S2 DEMs | PC VScode | `scripts/data4_mergeS2DEMs.ipynb` | large DEM over northern Myanmar
5|Create the terminal commands for running sen2cor | PC VScode | `SEAL/Ian/sen2cor_rmpi/data5_createS2Orders.R`; `data5_createS2Orders.csh` | text file with terminal commands
6|Batch run sen2cor, convert images from L1C to L2A, calculate NDVI, and delete L1C images to save memory. There is both a manual and auto version depending on how want to submit HPC jobs| HPC job | `SEAL/Ian/sen2cor_rmpi/data6_runSen2Cor_batch*.R`; `data6_runSen2Cor_batch*.csh`|L2A images
7|Check progress of sen2cor |PC| `SEAL/Ian/sen2cor_rmpi/data7_trackSen2CorProgress.R` | analysis
7b|Calc NDVI and delete L1C if this didn't happen within the sen2cor script | HPC | `SEAL/Ian/scripts/ch1_ndviL2A.R`; `SEAL/Ian/ch1_ndviL2A.csh` | L2A tifs
8|Extract data from L2A images and put in same format as downloaded GEE data | HPC job | `SEAL/Ian/data8_extractS2Bands.csh`; `scripts/data8_extractS2Bands.R`; `scripts/funs/dataA_extractS2.R` |csv
9|Identify forest strata covered by focal locations | PC | `scripts/data9_forestType.R` + `funs/dataB_maskForest.R` | cropped tif of strata
10|Identify other regions for landscape application | Either | `scripts/data10_altApps.R` + `funs/dataC_idRegions.R` | Rdata file of extents
| | | | |
||**Step 3: Training Exploratory: Compare different sensor combinations for training data**|
0| Run through 1 once for L8S2 + All, then run 1b for 30,60,90 days for span of lambda from 0-1 (0, 0.025, 0.05, 0.1, 0.2, 0.5, 1). To do so, set the window in `args.R` before each run of 1b.
1|Pre-process training data, compute ts models, calculate residuals, aggregate to 1 ts | Either | `scripts/train1_getResiduals.R` | Rdata files
2|Calculate seasonality adjustment based on raw z-scores | Either | `scripts/train2_sdInflation.R` | vector in Rdata file
3|Re-run step 1 now using the seasonal adjustment. | Either | `scripts/train1_getResiduals.R` | Rdata files
4|Calculate monitoring discount factor to use in Step 6-2 and landscape. | Either | `scripts/train2_sdInflation.R` | number; update in `args.R`
5|Define probability functions for each lambda per sensor combination | Either | `scripts/train3_createProbFuns.R` | .Rdata files
6|Calculate ewma, logMod, and probs for a vector of possible lambdas | Either | `scripts/train4_runLambdas.R` | csv files (separate for each lambda)
7|Prep accuracy metrics for using a sensitivity analysis over all thresholds| Either | `scripts/train5_metricsThresholds.R` | csv file
8|Calculate overall accuracy metrics and plot comparison panels| Either | `scripts/train6_accuracy.R` | csv files
9|Choose a threshold based on output from Step 8 and re-run Step 7 using that singular threshold. | |`scripts/train5_metricsThresholds.R` |
10|Plot results from Step 9 to visualize F1, PR, and lag differences btwn lambdas and sensors | Either | `scripts/train7_plotAccuracy.R` | plots
11|Choose the best lambda and sensor combination for landscape application based on output of Step 10|
12|Run k-fold cross-validation for chosen sensor combo and chosen lambda. | Either | `scripts/train8_crossVal.R` | plots
| | | | |
||**Step 4: Training: Apply chosen lambda / sensor combo, and create summary plots**|
0|Use best lambda and sensor combination to make logistic model for landscape application (non-backfilled data) | Either | `scripts/train6_createLogMod.R` | .Rdata file
0|^Note as of March 21 2023 we have abandoned log mod method |||
2|Download planet images for diagnostic plots | PC (faster) on VScode | `scripts/train9_getPlanetImagery.R` | tif files
3|Make diagnostic plots (8) for training data, including PLANET before/after images | Either | `train10_plotDiagnostics.R` | plots
| | | | |
||**Step 5: Application: Obtain data**|
1|If not already done from training data, download S2 images for Chatthin and convert using sen2cor (see above)|
2|Batch export S1 ARD images to Google drive folder | GEE | `app0A_exports1ARD.js` | S1 images
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
