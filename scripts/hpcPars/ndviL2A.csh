#!/bin/tcsh
#BSUB -n 16                     # specify number of cores per node (either 16-16-64000 or 32-32-130000)
#BSUB -W 240                   # request 4 days for running if 1000 imgs per rds file. Otherwise, do less.
#BSUB -R "rusage[mem=70000]"    # max mem on previous jobs was 69GB and 124.9 GB
                                ## multiplied by the number of cores per node requested (in MB) = 64000
#BSUB -R "select[model==E52640v2 || model==E52650v2 || model==E52650v3 || model==E52650v4 || model==E52680v2 || model==E52690v2]"
#BSUB -x                        # exclusive use of node
#BSUB -J ndviCon        # specify job NAME
#BSUB -o jobs/ch1_ndviConOut            # output file (%J will be replaced by job name [number])
#BSUB -e jobs/ch1_ndviConErr           # error file (%J will be replaced by job name [number])

## to run this script, first cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
## then do bsub < ndviL2A.csh
conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3
Rscript ./dissertation/scripts/ndviL2A.R 
conda deactivate