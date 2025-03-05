#!/bin/bash


#===============================================================================

#SBATCH -A                       ### account to run job under
#SBATCH -J S3_TRENDYv12_SDGVM          ### jobname 
#SBATCH -p par-multi                    ### partition to run job in (sinfo)
#SBATCH --mem=64000                    ### <number>[unit], when 0, all available memory requested
#SBATCH -N 4                        ### nodes requested for job
#SBATCH -n 8                           ### tasks job will run per node
#SBATCH -c 1                              ### cpu cores each task requires
#SBATCH -t 48:00:00                       ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcguirepatr@gmail.com
#SBATCH -o /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv12/S3//logs/TRENDYv12_S3.o%j

#===============================================================================


cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv12/S3/

echo "`date` -- SDGVM EXECUTION HAS STARTED" 
./run_SDGVM.bs
echo "`date` -- SDGVM EXECUTION HAS FINISHED" 

./run_postprocessing.bs
echo "`date` -- PP EXECUTION HAS FINISHED" 



### END ###
