#!/bin/bash


#===============================================================================

##SBATCH -A                       ### account to run job under
#SBATCH -J FS5_pp10CGVM          ### jobname 
#SBATCH -p par-single                    ### partition to run job in (sinfo)
#SBATCH --mem=64000                    ### <number>[unit], when 0, all available memory requested
##SBATCH --dependency=afterok:42259901     ### wait for FS5 run to finish
#SBATCH -N 1                        ### nodes requested for job
#SBATCH -n 8                           ### tasks job will run per node
#SBATCH -c 1                              ### cpu cores each task requires
#SBATCH -t 08:00:00                     ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcguirepatr@gmail.com
#SBATCH -o /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v2/FS5/logs/TRENDYv12_FS5_v10c3.o%j
#SBATCH -e /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v2/FS5/logs/TRENDYv12_FS5_v10c3.e%j

#===============================================================================

module load jasr
cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v2/FS5/

echo "`date` -- PP EXECUTION HAS STARTED" 
./run_postprocessing_v10c3.bs
echo "`date` -- PP EXECUTION HAS FINISHED" 


### END ###
