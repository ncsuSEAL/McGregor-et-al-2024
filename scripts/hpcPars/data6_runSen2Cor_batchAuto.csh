#!/bin/csh
#
# Script:  R_loops.csh
# Usage: For submitting multiple batch R jobs to LSF,
#     an example of looping over years and models
# Author:  NCSU HPC
#
# 
## To run this, first cd to the folder here and change the "6" to however many
## command RDS files we have. Then copy paste the csh line into the terminal.
#     cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
#     ./sen2cor_rmpi/data6_runSen2Cor_batchAuto.csh 10
#  Note that script must have execute permissions, i.e.,
#     chmod u+x R_loops.csh

## IF USING RMPI
### do not include -x in the submitted job. Otherwise, if not using rmpi, then it is necessary
# 

## Time needed based on how many commands are in each rds file (from data5_createOrders.R)
### 1000 images = 5760 (4 days)
### 500 images = 2880 (48 hours)
### 250 images = 1440 (24 hours)
### 125 images = 720 (12 hours)
### 63 images = 360 (6 hours)
### 32 images = 180 (3 hours)

cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor

# module load openmpi-gcc/openmpi1.8.4-gcc4.8.2
conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3

# Specify number of jobs to submit from the first argument above
set nIter = $1

# Initialize iteration loop counter 
set iter = 1

# Full loop
while ($iter <= $nIter)
  echo "Submitting job iter = $iter"

  bsub -n 16 -W 180 -R "rusage[mem=64000]" -R "select[model==E52640v2 || model==E52650v2 || model==E52650v3 || model==E52650v4 || model==E52680v2 || model==E52690v2]" -x -J sen2corBatch$iter -oo sen2cor_rmpi/jobRunSen2Cor/sen2corOutBatch=$iter -eo sen2cor_rmpi/jobRunSen2Cor/sen2corErrBatch=$iter "Rscript ./sen2cor_rmpi/data6_runSen2Cor_batch.R $iter"

  ## when running Rmpi, there should NOT be -x
  # bsub -n 16 -W 2880 -R "rusage[mem=64000]" -R "select[model==E52640v2 || model==E52650v2 || model==E52650v3 || model==E52650v4 || model==E52680v2 || model==E52690v2]" -J sen2corBatch$iter -oo sen2cor_rmpi/jobRunSen2Cor/sen2corOutBatch=$iter -eo sen2cor_rmpi/jobRunSen2Cor/sen2corErrBatch=$iter "mpirun -n 1 Rscript ./sen2cor_rmpi/data6_runSen2Cor_batchRmpi.R $iter"

  @ iter++
  
end

conda deactivate