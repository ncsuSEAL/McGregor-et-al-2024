#!/bin/tcsh
#BSUB -n 16                     # specify number of cores per node
#BSUB -W 30                   # request 4 days (5760) for running (max for general HPC)
#BSUB -R "rusage[mem=16GB]"    # memory allocation based on previous job submissions
#BSUB -R "span[hosts=1]"
#BSUB -q short
##BSUB -x                        # exclusive use of node
#BSUB -J dailyDistRatio1      # specify job NAME
#BSUB -o jobs/app3out              # output file (%J will be replaced by job name [number])
#BSUB -e jobs/app3err              # error file (%J will be replaced by job name [number])

# cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
# bsub < app3_dailyProbMaps.csh
conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3
Rscript ./dissertation/scripts/app3_dailyProbMaps.R
conda deactivate