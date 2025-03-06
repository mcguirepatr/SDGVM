Notes:  
  The version numbers below refer to:
      v13  = TRENDYv13 (2024)
      v10c = code version number for the current version of input/output files 
      v8y  = code version number for the current version of SDGVM FORTRAN code
      v9c  = template installation version number from TRENDYv12 to copy from 
      v11c = installation version number to copy from 
      v11d = installation version number to copy to 

#Steps required to install and run the code on JASMIN:
 #open a screen session & rsync the db:
   screen
   rsync -auh --info=progress2 /gws/nopw/j04/nexcs/pmcguire/TRENDYv13/db/ /work/scratch-pw2/pmcguire/TRENDYv13/db/
   #quit the screen session with CTRL-A-D, as the rsync might take a while
   #we can check the completion of the rsync later with screen -r

 #Create & modify the runs in the landsurf_rdg group workspace 
 cd /gws/nopw/j04/landsurf_rdg/pmcguire/sdgvmR2
 cp -pr TRENDYv13_v11c TRENDYv13_v11d
 cd TRENDYv13_v11d 
 ./create_runs.scr
 ./modify_runs.scr

 #Copy the runs to the much faster scratch-pw2 HDD & run them:
 cp -pr ../TRENDYv13_v11d/  /work/scratch-pw2/pmcguire/sdgvmR/.
 cd /work/scratch-pw2/pmcguire/sdgvmR/TRENDYv13_v11d/
 ./run_runs.scr
