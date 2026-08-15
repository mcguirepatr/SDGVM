#!/bin/bash


#===============================================================================

##SBATCH -A                       ### account to run job under
#SBATCH -J FS4_pp10CGVM          ### jobname 
#SBATCH -p standard                    ### partition to run job in (sinfo)
#SBATCH --qos=high                 ### partition to run job in (sinfo)
#SBATCH --account=landsurf_rdg
#SBATCH --mem=64000                    ### <number>[unit], when 0, all available memory requested
#SBATCH -N 1                        ### nodes requested for job
#SBATCH -n 8                           ### tasks job will run per node
#SBATCH -c 1                              ### cpu cores each task requires
#SBATCH -t 24:00:00                     ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcguirepatr@gmail.com
#SBATCH -o /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9e/FS4/logs/TRENDYv12_FS4_v10c.o%j
#SBATCH -e /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9e/FS4/logs/TRENDYv12_FS4_v10c.e%j

#===============================================================================

module load jasr
cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9e/FS4/

echo "`date` -- PP EXECUTION HAS STARTED" 
./run_postprocessing_v10c.bs
echo "`date` -- PP EXECUTION HAS FINISHED" 


### END ###
