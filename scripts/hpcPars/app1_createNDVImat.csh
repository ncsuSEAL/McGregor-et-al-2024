#!/bin/csh
#
# Script:  R_loops.csh
# Usage: For submitting multiple batch R jobs to LSF,
#     an example of looping over years and models
# Author:  NCSU HPC
#
# 
## To run this (and submit 18 identical jobs), type the following in terminal after cd to my root folder in HPC
## Note here we are doing 18 straight tasks, with the actual nrows and start row hard-coded in the script itself.
# cd to the right folder first: cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
#
# NOTE this script MUST be in root folder (above) because of the source file paths used by the VRTs
#
#     ./app1_createNDVImat.csh 36 #for 60 rows at a time. this should be changed to 18 if doing 120 rows at a time
#              Note that script must have execute permissions, i.e.,
#     chmod u+x R_loops.csh

conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3

# Specify number of jobs to submit
set nIter = $1

# Initialize iteration loop counter 
set iter = 1

while ($iter <= $nIter)
  echo "Submitting job iter = $iter"

  bsub -n 16 -W 60 -R "rusage[mem=16000]" -x -J land1_iter=$iter -oo jobs/app1_createNDVImat/outLand1_iter=$iter -eo jobs/app1_createNDVImat/errLand1_iter=$iter "Rscript ./dissertation/scripts/app1_createNDVImat.R $iter"

  @ iter++
  
end

conda deactivate