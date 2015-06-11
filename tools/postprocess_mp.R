######################
#
# Read parallel SDGVM output and stich together
#
# AWalker
# Oct 2013
#
######################

rm(list=ls())

# Initialisation
########################

library(parallel)

dir <- '~/models/SDGVM/'
fd  <- paste(dir,'src/sdgvm/tools/',sep='')
setwd(fd)
source('read_SDGVM_output.R')
source('functions_netcdf.R')

wd <- '/mnt/disk2/Research_Projects/leaf_trait_SDGVM/simulations/TERRABITES/'

#number of parallel grid directories 
grids  <- 15

#type of data to process
stich   <- T
netcdf  <- T
deg1    <- T

#file names etc
fend   <- '.dat'
ncf    <- 'SDGVM'
ncfend <- '.nc'
sim    <- c('T1_131001_GLC2000','T1_131001_GLOB2009','T1_140129','T1_noC4','T1','T2','T3')
pia    <- 7
# pia    <- 5

# nc file parameters
# naming convention SDGVM_Tx_v.nc
# one file per variable
mis_val <- -99999


### Start Program
##########################
wd_list <- paste(wd,sim,'/output/',sep='')
if(deg1) wd_list <- paste(wd,sim,'/output_1deg/',sep='')


# stich mp data
if(stich) lapply(wd_list[pia],stich_sdgvm_mp_apply,
                 grids=grids,annual=F,monthly=T,daily=F,procs=4)

# convert data to CMOR netcdf output
nyears <- 110
nsites <- 1548
lon    <- 3.75
lat    <- 2.5

if(deg1){
  nsites <- 15417
  lon    <- 1
  lat    <- 1  
}

if(netcdf) lapply(wd_list[pia],write_sdgvm_netcdf,
                  afiles=afiles,
                  mfiles=mfiles,
                  nsites=nsites,nyears=nyears,lon=lon,lat=lat)

  
  
  
  
  
  
  
  
  
  
  
