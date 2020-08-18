#!/bin/bash


#===============================================================================

#SBATCH -A ccsi                           ### account to run job under
#SBATCH -J SDGVM_#PROJECT#_archive        ### jobname 
#SBATCH -p batch                          ### partition to run job in (sinfo)
#SBATCH --mem=0                           ### <number>[unit], when 0, all available memory requested
#SBATCH -N 1                              ### nodes requested for job
#SBATCH -n 32                             ### tasks job will run per node
#SBATCH -c 1                              ### cpu cores each task requires
#SBATCH -t 05:00:00                       ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=walkerap@ornl.gov
#SBATCH -o #PROJDIR#/tools/logs/%x.o%j

#===============================================================================


cd #PROJDIR#

echo "`date` -- SDGVM nc ZIP HAS STARTED" 
./gzip_ncdf.bs
echo "`date` -- SDGVM nc ZIP HAS FINISHED" 

echo "`date` -- SDGVM ARCHIVE HAS STARTED" 
./archive.bs
echo "`date` -- SDGVM ARCHIVE HAS FINISHED" 



### END ###
