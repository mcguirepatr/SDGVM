#!/bin/bash


#===============================================================================

#SBATCH -A #ACCOUNT#                      ### account to run job under
#SBATCH -J #SIM#_#PROJECT#_SDGVM          ### jobname 
#SBATCH -p #PARTITION#                    ### partition to run job in (sinfo)
#SBATCH --mem=0                           ### <number>[unit], when 0, all available memory requested
#SBATCH -N #NODES#                        ### nodes requested for job
#SBATCH -n #NP#                           ### tasks job will run per node
#SBATCH -c 1                              ### cpu cores each task requires
#SBATCH -t 02:00:00                       ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#EMAIL#
#SBATCH -o #RUNDIR#/logs/#PROJECT#_#SIM#.o%j

#===============================================================================


cd #RUNDIR#

echo "`date` -- PP EXECUTION HAS STARTED" 
./run_postprocessing.bs
echo "`date` -- PP EXECUTION HAS FINISHED" 


### END ###
