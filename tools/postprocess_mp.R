######################
#
# Read parallel SDGVM output and stich together
#
# AWalker
# Jun 2015
#
######################

rm(list=ls())



# Set default argument values
########################

# output timesteps to process
annual  <- T
monthly <- T
daily   <- T

# stich output files from all grids into a single file
stich   <- T

# delete sub-grid files once processesed
delete  <- T 

# write CMOR netcdf output
netcdf  <- F

# simulation res 1x1 degree, F - 3.75 x 2.5
deg1    <- T

# main directory
dir  <- '~/models/SDGVM/'

# source code directory
fd   <- paste(dir,'src/sdgvm/tools/',sep='/')

# model project directory
wd   <- paste(dir,'run/',sep='/')

# simulation directory
sim  <- c('blank/')

# index array to index the above 'sim' vector
pia  <- 1

# start year and number of years of data in daily and monthly output
sty  <- 1901
ny   <- 112

# number of parallel grid directories 
grids  <- 30

# number of cores to run the analysis over 
cores  <- 32

#file names etc
fend   <- '.dat'
ncf    <- 'SDGVM'
ncfend <- '.nc'

# netcdf files to create
ncdf_avars <- c('cVegpft')
ncdf_mvars <- c('pr','evapotranspft')
#ncdf_avars <- c('cVeg','cLitter','cSoil','cVegpft','fFire','fLuc','cLeaf','cRoot','burntArea')
#ncdf_mvars <- c('tas','pr','rsds','mrro','evapotrans','gpp','ra','npp','rh','nbp','lai',
#                'evapotranspft','transpft','swepft','gpppft','npppft','nbppft','tran')
#'mrso','landCoverFrac'

# nc file parameters
mis_val  <- -99999
nsites   <- 1548
lon      <- 3.75
lat      <- 2.5
pftnames <- c('BARE','CITY','C3','C3crop','C4','C4crop','Dc_Bl','Dc_Nl','Ev_Bl','Ev_Nl')

if(deg1){
  nsites <- 15417
#  nsites <- 18316
  lon    <- 1
  lat    <- 1  
}



### Parse command line arguments   
##########################
# any one of the above objects can be specified as a command line argument using the syntax:
# Rscript <nameofthisscript> "<object1><-<value1>" "<object2><-<value2>"
# e.g. Rscript postprocess_mp.R "sty<-2000" "dir<-'/home/alp/models/SDGVM'"
if(length(commandArgs(T))>=1) {
  for( ca in 1:length(commandArgs(T)) ) {
    eval(parse(text=commandArgs(T)[ca]))
  }
}



### Read in functions 
##########################
setwd(fd)
source('read_SDGVM_output.R')



### Start Program
##########################
wd_list <- paste(wd,sim,'output/',sep='/')
print(wd_list)

# stich mp data
if(stich) {
  lapply(wd_list[pia],stich_sdgvm_mp_apply,
         grids=grids,mc=T,
         annual=annual,monthly=monthly,daily=daily,
         mc.cores=cores,styear=sty,nyears=ny)

  # remove grid output now that combined files have been created
  if(delete) {
    for( wdc in wd_list[pia]) {
      setwd(wdc)
      system("for i in grid*; do mv $i/simulation.dat $i/simulation.txt; done")
      system("for i in grid*; do mv $i/site_info.dat  $i/site_info.txt; done")
      system("for i in grid*; do mv $i/diag.dat       $i/diag.txt; done")
      system("rm ./grid*/*.dat")
  }} 
}

# convert data to CMOR netcdf output
if(netcdf) source('functions_netcdf.R')
if(netcdf) lapply(wd_list[pia],write_sdgvm_netcdf,
                  afiles=ncdf_avars,
                  mfiles=ncdf_mvars,
                  nsites=nsites,nyears=ny,styr=sty,lon=lon,lat=lat)
  
