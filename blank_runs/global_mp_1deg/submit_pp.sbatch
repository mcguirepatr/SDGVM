#!/bin/bash


#===============================================================================

#SBATCH --account=cades-ccsi  ### appears redundant with below
#SBATCH -A ccsi               ### account to run job under
#SBATCH -p testing            ### partition to run job in (sinfo)
#SBATCH --mem=0               ### <number>[unit], when 0, all available memory requested
#SBATCH -N 1                  ### nodes requested for job
#SBATCH -n 32                 ### tasks job will run per node
#SBATCH -c 1                  ### cpu cores each task requires
#SBATCH -t 00:10:00           ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walkerap@ornl.gov
#SBATCH -o #RUNDIR#/logs/${SLURM_JOB_NAME}.o%j

#===============================================================================
#setenv OMP_NUM_THREADS 1



cd #RUNDIR#

./run_postprocessing.bs
echo "`date` -- PP EXECUTION HAS FINISHED" 


### END ###
