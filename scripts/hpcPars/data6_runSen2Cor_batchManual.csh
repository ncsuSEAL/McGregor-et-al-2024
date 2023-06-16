#!/bin/tcsh
#BSUB -n 16                     # specify number of cores per node (either 16-16-64000 or 32-32-130000)
#BSUB -W 1440                   # request 4 days for running if 1000 imgs per rds file. Otherwise, do less.
#BSUB -R "rusage[mem=64000]"    # memory allocation of 4GB (approximate memory needed for processing 1 L1C image) 
                                ## multiplied by the number of cores per node requested (in MB) = 64000
#BSUB -R "select[model==E52640v2 || model==E52650v2 || model==E52650v3 || model==E52650v4 || model==E52680v2 || model==E52690v2]"
#BSUB -x                        # exclusive use of node
#BSUB -J sen2cor_batchNew_6        # specify job NAME
#BSUB -o sen2cor_rmpi/jobRunSen2Cor/sen2corOutBatch=6            # output file (%J will be replaced by job name [number])
#BSUB -e sen2cor_rmpi/jobRunSen2Cor/sen2corErrBatch=6            # error file (%J will be replaced by job name [number])
# setenv TMP /scratch
# setenv TMPDIR /scratch
# setenv TEMP /scratch

## to run this script, first cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
## then do bsub < sen2cor_rmpi/data6_runSen2Cor_batchManual.csh
conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3
Rscript sen2cor_rmpi/data6_runSen2Cor_batch.R 
conda deactivate

## bsub -n 16 -W 5760 -R "rusage[mem=64000]" -x "select[model==E52640v2 || model==E52650v2 || model==E52650v3 || model==E52650v4 || model==E52680v2 || model==E52690v2]" -J sen2corBatch$iter -oo sen2cor_rmpi/jobRunSen2Cor/out -eo sen2cor_rmpi/jobRunSen2Cor/err mpirun -n 1 "Rscript sen2cor_rmpi/data6_runSen2Cor_batchManual.R" 