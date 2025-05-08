Notes:  
  The version numbers below refer to:
      v13  = TRENDYv13 (2024)
      v10c = code version number for the current version of input/output files 
      v8y  = code version number for the current version of SDGVM FORTRAN code
      v9c  = template installation version number from TRENDYv12 to copy from 
      v11c = installation version number to copy from 
      v11d = installation version number to copy to 

#Steps required to install and run the code on JASMIN:
 #Copy or checkout the repo https://bitbucket.org/mcguirepatrUReading/sdgvm_transitions2/src/bugfix_negsoilc/
 #into the following directory:
  /gws/nopw/j04/nexcs/pmcguire/TRENDYv10/sdgvm_v8y_debug/

 #open a screen session & rsync the db:
   screen
   rsync -auh --info=progress2 /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/ /work/scratch-pw2/pmcguire/TRENDYv13/db/
   #quit the screen session with CTRL-A-D, as the rsync might take a while
   #we can check the completion of the rsync later with screen -r

 #Prior to this rsync, you might want to update some of the files in the db. For example, for TRENDYv12 & TRENDYv13, we have updated the landcover database
 # so that it uses the NETCDF file directly, but binning by two is required beforehand.
 #These files have already been downloaded and processed for 2024, and are available in /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/
      cd /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/landcover/hyde/LUH2/v2h_2024/
      #copy your states.nc and transitions.nc from the website to the directory above
 #The half-resolution version is computed with:
      module load jaspy # on JASMIN
      cdo gridboxmean,2,2 states.nc states2b.nc
      cdo gridboxmean,2,2 transitions.nc transitions2b.nc

 #Also, prior to this rsync, make sure that you have downloaded and done the R preprocessing of the JRU/CRA data.
 #  The 2024 copy of the preprocessed data is found here:
    /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/clim/global/crujra_2024/30min
 #The preprocessing scripts are found here:
  cd /gws/nopw/j04/nexcs/pmcguire/TRENDYv10/sdgvm_v8y_debug/reformat/
 # The current preprocessing script for this has a few most-recent versions, i.e.:
  process_CRUJRAmet_4SDGVM-2_v2b4.R 
  process_CRUJRAmet_4SDGVM-2_v2c.R 
 #And they can be run with:
  screen
  module load jasr
  Rscript process_CRUJRAmet_4SDGVM-2_v2c.R  > process_CRUJRAmet_4SDGVM-2_v2c.log 2> process_CRUJRAmet_4SDGVM-2_v2c.err &
  #quit the screen session with CTRL-A-D, as the preprocessing  might take a while
  #we can check the completion of the preprocessing later with screen -r


 #compile the SDGVM FORTRAN code:
  cd /gws/nopw/j04/nexcs/pmcguire/TRENDYv10/sdgvm_v8y_debug/
  module load jaspy
  make sdgvm0
  module unload jaspy

 #Create & modify the runs in the landsurf_rdg group workspace 
 cd /gws/nopw/j04/landsurf_rdg/pmcguire/sdgvmR2
 #this contains the template installation from TRENDYv12
 cp -pr /gws/nopw/j04/nexcs/pmcguire/TRENDYv10/sdgvm_v8y_debug/TRENDYv13_v9c TRENDYv13_v9c
 #this contains the installation for TRENDYv13
 cp -pr /gws/nopw/j04/nexcs/pmcguire/TRENDYv10/sdgvm_v8y_debug/TRENDYv13_v11c TRENDYv13_v11d
 cd TRENDYv13_v11d 
 ./create_runs.scr #create the TRENDYv12 runs in the (spin_accel,spin_accelFS,S0,S1,S2,S3,S5,FS4,FS5) dirs
 ./modify_runs.scr #update the runs for TRENDYv13

 #Copy the runs to the much faster scratch-pw2 HDD & run them:
 cp -pr ../TRENDYv13_v11d/  /work/scratch-pw2/pmcguire/sdgvmR/.
 cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v11d/
 ./run_runs.scr

 #We can monitor the runs:
 squeue -u YOUR_USER_NAME 

 #After the runs are finished (this could take a day or more), check the output files & look at the NETCDF files in some of the run directories:
 cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v11d/
 grep End S3/output_v10c/g*out* |wc #This should return 100; anything less than that means some of the 100 array jobs failed somewhere
 ncview S3/output_v10c/SDGVM_S3_cSoil.nc
 #repeat for the other runs (S0,S1,S2,FS4,FS5,spin_accel,spin_accelFS)

 #run the results plotting script to compare runs
 cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v11d/results
 screen 
 ./plot_TRENDY_v13c.bs > logs/plot_TRENDY_v13c.log 2> logs/plot_TRENDY_v13c.err & 
 #quit the screen session with CTRL-A-D, as the plotting might take a while
 #we can check the completion of the plotting script later with screen -r
 #The output plots are in the 250306a directory (as directed by the .bs script), and can be viewed (or downloaded).
 #Viewing is possible with the display command.

