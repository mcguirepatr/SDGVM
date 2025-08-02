#!/bin/bash


#===============================================================================

#SBATCH -A #ACCOUNT#                      ### account to run job under 
#SBATCH -J #SIM#_pp_SDGVM                 ### jobname 
#SBATCH -p #PARTITION#                    ### partition to run job in (sinfo)
#SBATCH --mem=#MEMORY#                    ### <number>[unit], when 0, all available memory requested
#SBATCH -N #NODES#                        ### nodes requested for job
#SBATCH -n #NP#                           ### tasks job will run per node
#SBATCH -c 1                              ### cpu cores each task requires
#SBATCH -t #WALLTIME#                     ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#EMAIL#
#SBATCH -o #RUNDIR#/logs/#PROJECT#_#SIM#_pp.o%j

#===============================================================================


cd #RUNDIR#

echo "`date` -- PP EXECUTION HAS STARTED" 
./run_postprocessing.bs
echo "`date` -- PP EXECUTION HAS FINISHED" 


### END ###
