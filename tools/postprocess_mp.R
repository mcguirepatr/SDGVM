######################
#
# Read parallel SDGVM output and stich together
#
# AWalker
# Jun 2015
#
######################

rm(list=ls())



# Set default object values
########################

# output timesteps to process
annual  <- T
monthly <- T
daily   <- T

# stich output files from all grids into a single file
stich   <- T

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

# nc file parameters
mis_val <- -99999
nsites  <- 1548
lon     <- 3.75
lat     <- 2.5

if(deg1){
  nsites <- 15417
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
#source('functions_netcdf.R')



### Start Program
##########################
wd_list <- paste(wd,sim,'output/',sep='/')

# stich mp data
if(stich) lapply(wd_list[pia],stich_sdgvm_mp_apply,
                 grids=grids,mc=T,
                 annual=annual,monthly=monthly,daily=daily,
                 mc.cores=cores,styear=sty,nyears=ny)

# convert data to CMOR netcdf output
if(netcdf) lapply(wd_list[pia],write_sdgvm_netcdf,
                  afiles=afiles,
                  mfiles=mfiles,
                  nsites=nsites,nyears=ny,lon=lon,lat=lat)

  
  
  
  
