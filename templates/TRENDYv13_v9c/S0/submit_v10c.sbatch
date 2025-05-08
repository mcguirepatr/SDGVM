#!/bin/bash


#===============================================================================

##SBATCH -A                       ### account to run job under
#SBATCH --job-name=S010CNDY2024_SDGVM          ### jobname 
#SBATCH -p standard                    ### partition to run job in (sinfo)
#SBATCH --qos=standard                  ### partition to run job in (sinfo)
#SBATCH --account=landsurf_rdg
#SBATCH --array=1-100                      ### tasks job will run 
#SBATCH --ntasks-per-node=1
#SBATCH -t 24:00:00                       ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcguirepatr@gmail.com
#SBATCH -o /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9c/S0/logs/TRENDY2024_S0_v10c.o%A-%a
#SBATCH -e /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9c/S0/logs/TRENDY2024_S0_v10c.e%A-%a

#===============================================================================


cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9c/S0 

module load jaspy/3.7/r20210320

echo "`date` -- SDGVM EXECUTION HAS STARTED" 
#./run_SDGVM.bs

# script to run SDGVM according to the sites, resolution, and cores specified in sites.bs 
# designed to be run from within the directory of a specific, single run

MODEL=/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/sdgvm_v8t_debug/
INPUT=input_global_v5a10g10cS0.dat

SPINDIR='/work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9c/spin_accel/output_v10c'
OUTPUT='/work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v9c/S0/output_v10c'


# start script
# clear output directories
#rm -r ./$OUTPUT
#rm -r ./src


# copy source code and compiled model run file
#mkdir ./$OUTPUT
#mkdir ./src
#cp $MODEL/*.f    ./src
#cp $MODEL/sdgvm0 ./src


# sites variables
#SDGVM_SITES=/gws/nopw/j04/nexcs/pmcguire/sdgvmD/data//sites/global_1deg_CRUJRA2018.dat
#SDGVM_SITES=/gws/nopw/j04/nexcs/pmcguire/sdgvmD/data//sites/global_1deg_CRUNCEPv7.dat
#SDGVM_N_SITES=( 1980 625 620 610 615 515 515 515 500 510 575 575 560 430 500 500 520 525 550 575 530 630 540 540 455 445 490 510 400 400 550 511 )
#SDGVM_SITES=/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/CRUJRA2021/SDGVMdata/monthly/0.5deg_v2a5/lat_lon_unformatted_1deg.dat #this file name is too long
SDGVM_SITES=/work/scratch-pw2/pmcguire/TRENDYv13/db/sites/global_1deg_CRUJRA2024.dat
#SDGVM_SITES=/gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/sites/global_1deg_CRUJRA2024.dat
#SDGVM_N_SITES=( 1980 625 620 610 615 515 515 515 500 510 575 575 560 430 500 500 520 525 550 575 530 630 540 540 455 445 490 510 400 400 550 511 )
SDGVM_N_SITES=( \
199 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 \
183 183 183 183 183 183 183 183 183 183 ) 

#
# run loop 
s4=0
((e=${#SDGVM_N_SITES[*]}))
s=${SLURM_ARRAY_TASK_ID}
for t in `seq 1 $s` ;do
 ((s1=$t - 1))
 ((s3=1 + $s4))
 ((s4=$s4 + ${SDGVM_N_SITES[$s1]}))
done

mkdir $OUTPUT/grid$s
echo "./src/sdgvm0 $INPUT $SPINDIR/grid$s $OUTPUT/grid$s  $SDGVM_SITES $s3 $s4 > $OUTPUT/g${s}out.txt 2> $OUTPUT/g${s}err.txt"
      ./src/sdgvm0 $INPUT $SPINDIR/grid$s $OUTPUT/grid$s  $SDGVM_SITES $s3 $s4 > $OUTPUT/g${s}out.txt 2> $OUTPUT/g${s}err.txt


echo "`date` -- SDGVM EXECUTION HAS FINISHED" 

#./run_postprocessing.bs
#echo "`date` -- PP EXECUTION HAS FINISHED" 
### END ###
