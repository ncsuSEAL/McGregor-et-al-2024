#!/bin/tcsh
#BSUB -n 11              	# request 11 cores
#BSUB -W 60            	# request 60 minutes for running
#BSUB -q shared_memory  	# ensure all cores are in same node
#BSUB -J get_s2_orders         	# specify job ID
#BSUB -o jobS2Orders/stdout.%J      	# output file (%J will be replaced by job name)
#BSUB -e jobS2Orders/stderr.%J      	# error file (%J will be replaced by job name)

# conda activate /Your/Path # if you'd like, but it's not necessary
module load R
R CMD BATCH --vanilla data5_createS2Orders.R