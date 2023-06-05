# Workflow steps with corresponding script names and locations

## Description
This is the same workflow as in the flowchart I've created but the associated scripts are included along with their general location. Each section's steps are in order.

Medium = Mac, PC, Either, or HPC
- "Either" refers to either Mac or PC

Step #|Task|Medium|Script|Output
---|---|---|---|---|
||**Step 1: Identify and download training data**|
1| NA since already have data
| | | | |
||**Step 2: Process S2 using Sen2Cor**|
1| NA since already have data
| | | | |
||**Step 3: Training Exploratory: Compare different sensor combinations for training data**|
1| Run through 1a once for L8S2 + All, then run 1b for 60 days for span of lambda from 0-1 (0, 0.025, 0.05, 0.1, 0.2, 0.5, 1). To do so, set the window in `args.R` before each run of 1b.
1a|Pre-process training data, compute ts models, calculate residuals, aggregate to 1 ts | Mac is fine, could do HPC job | `scripts/train1_getResiduals.R`; `SEAL/Ian/diss/train1_getResiduals.csh` | Rdata files
1b|Calculate seasonality adjustment based on raw z-scores | Either | `scripts/train2_sdInflation.R` | vector in Rdata file
1c|Re-run step 1a now using the seasonal adjustment. | Either | `scripts/train1_getResiduals.R` | Rdata files
1d|Calculate monitoring discount factor to use in 1e and landscape. | Either | `scripts/train2_sdInflation.R` | number; update in `args.R`
1e|Define probability functions for each lambda per sensor combination | Either | `scripts/train3_createProbFuns.R` | .Rdata files
1f|Calculate ewma, logMod, and probs for a vector of possible lambdas | Either | `scripts/train4_runLambdas.R` | csv files (separate for each lambda)
1g|Calculate accuracy metrics for all lambdas and sensors | Either | `scripts/train5_metricsThresholds.R` | csv files
1h|Plot results to visualize F1, PR, and lag differences btwn lambdas and sensors | Either | `scripts/train6_plotAccuracy.R` | plots
2|Run Step 1c again using new lambda vector and best window based on output plots from 1h.
3|Using output of step 2, compare F1 scores and other metrics via plots | Either | `scripts/train6_plotAccuracy` | plots
4|Choose the best lambda and sensor combination for landscape application|
5|Using output of step 3, run k-fold cross-validation for sensor combo and chosen lambda. | Either | `scripts/train7_crossVal.R` | plots
| | | | |
||**Step 4: Training: Apply chosen lambda / sensor combo, and create summary plots**|
0|Use best lambda and sensor combination to make logistic model for landscape application (non-backfilled data) | Either | `scripts/train6_createLogMod.R` | .Rdata file
0|^Note as of March 21 2023 we have abandoned log mod method |||
2|Download planet images for diagnostic plots | PC (faster) on VScode | `scripts/train8_getPlanetImagery.R` | tif files
3|Make diagnostic plots (8) for training data, including PLANET before/after images | Either | `train9_plotDiagnostics.R` | plots
| | | | |
||**Step 5: Application: Obtain data**|
1| NA
| | | | |
||**Step 6: Application: Process data**|
0| NA
| | | | |
||**Step 7: Application: Analyze results**|
1| NA
