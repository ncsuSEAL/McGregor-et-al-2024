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
# NOTE because app1 script MUST be in my root folder (above), I'm keeping this csh file there as well for
# consistency. Technically this could be in the dissertation folder.
#
#     ./app2_processLandProbs.csh 36 #this is for 60 rows at a time. Change to 18 if doing 120 rows at a time
#  Note that script must have execute permissions, i.e.,
#     chmod u+x R_loops.csh

setenv OMP_NUM_THREADS 1
conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3
 
# Specify number of jobs to submit
set nIter = $1

# Initialize iteration loop counter 
set iter = 1

# Mem of 65GB is for doing allProbs. 
# If doing bayesUpdate or normalTest alone, specify 80GB for 5k pts or 96GB for 33k
## MAKE SURE -n matches the nCores argument in the script!!!
while ($iter <= $nIter)
    echo "Submitting job iter = $iter"

    bsub -n 20 -W 60 -q short -R "rusage[mem=100GB]" -R "span[hosts=1]" -J land2_iter=$iter -oo jobs/app2_processLandProbs/outLand2_iter=$iter -eo jobs/app2_processLandProbs/errLand2_iter=$iter "Rscript ./dissertation/scripts/app2_processLandProbs.R $iter"
 
    @ iter++
    
end

conda deactivate