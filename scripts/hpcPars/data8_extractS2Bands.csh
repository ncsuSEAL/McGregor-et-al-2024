#!/bin/tcsh
#BSUB -n 16                    # specify number of cores per node
                                # TO RUN THIS FILE, cd to IanMcGregor/ and type
                                # bsub < data8_extractS2Bands.csh in the terminal
#BSUB -R "rusage[mem=8000]"    # this doesn't use a ton of memory since just doing NDVI calculation
#BSUB -R select[avx2]           # request one of the newer CPUs (otherwise will take >24 hours)
#BSUB -x                        # exclusive use of node
#BSUB -W 60                    # when using avx2, this only needs >2 <3 hours, otherwise will be >24
                                # can be 4 days (max for general HPC) or 10 days (14400 for CNR queue)
#BSUB -J extractS2              # specify job ID
#BSUB -o jobs/data8_extractS2Bands/out_extractS2_%J  # output file (%J will be replaced by job name [number])
#BSUB -e jobs/data8_extractS2Bands/err_extractS2_%J             # error file (%J will be replaced by job name [number])

## to run this, do cd /rsstu/users/j/jmgray2/SEAL/IanMcGregor
## then do bsub < data8_extractS2Bands.csh
conda activate /usr/local/usrapps/jmgray2/imcgreg/env_dissPareto3
Rscript ./dissertation/scripts/data8_extractS2Bands.R
conda deactivate
