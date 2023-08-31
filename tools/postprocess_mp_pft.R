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
daily   <- F

# stich output files from all grids into a single file, use multiple processors
stich   <- T
mc      <- T 

# delete sub-grid files once processesed
delete  <- F 

# write CMOR(ish) netcdf output
netcdf  <- F

# simulation res 1x1 degree, F - 0.5 x 0.5 
deg1    <- T

# main directory
#dir  <- '~/models/SDGVM/'
dir <- '/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/'

# source code tools directory
#fd   <- paste(dir,'src/sdgvm/tools/',sep='/')
fd   <- paste(dir,'sdgvm/tools/',sep='/')

# directory in which simulation directory lives
wd   <- paste(dir,'run/',sep='/')

# simulation directory
sim  <- c('blank/')

# index array to index the above 'sim' vector
pia  <- 1

# start year and number of years of data in SDGVM output files
sty  <- 1700 
ny   <- 321

# years of output requested for netcdf
outsyear <- NULL 
outeyear <- NULL 

# number of parallel grid directories 
grids  <- 32 

# number of cores to run the analysis over 
cores  <- 8 

#file names etc
fend   <- '.dat'
ncf    <- 'SDGVM'
ncfend <- '.nc'

root_dir <- '/work/scratch-pw/pmcguire/TRENDY2021_v1/'
output_dir <- 'output_6_7_500/'


# netcdf files to create
ncdf_avars <- c('cVeg','cLitter','cSoil','cRoot','burntArea','landCoverFrac')
#ncdf_avars <- c('cVeg','cLitter','cSoil','fFire','fLuc','cLeaf','cRoot','burntArea')
#ncdf_avars <- c('fFire','fLuc')
#ncdf_avars <- 'cVegpft' 
#ncdf_avars <- 'fLeach' 
#ncdf_avars <- 'pot_evapotrans' 
#ncdf_avars <- NULL 

#ncdf_mvars <- c('tas','pr')
#ncdf_mvars <- 'snow_depthpft' 
#ncdf_mvars <- c('gpppft','tran','npppft')
ncdf_mvars <- c('tas','pr','rsds','mrro','mrso','evapotrans','gpp','ra','npp','rh','nbp','lai',
                'evapotranspft','transpft','snow_depthpft','gpppft','npppft','tran','laipft')
#ncdf_mvars <- c('tas','pr','rsds','mrro','mrso','evapotrans')
#ncdf_mvars <- c('gpp','ra','npp','rh','nbp','lai',
#                'evapotranspft','transpft','snow_depthpft','gpppft','npppft','tran','landCoverFrac')
#ncdf_mvars <- 'nbp' 
#ncdf_mvars <- 'laipft' 

# nc file parameters
mis_val   <- -99999
nsites    <- 62220 
lon       <- 0.5 
lat       <- 0.5
pftnames  <- c('BARE','CITY','C3','C3crop','C4','C4crop','Dc_Bl','Dc_Nl','Ev_Bl','Ev_Nl')
#pftnames  <- c('BARE','CITY','C3')

# specifiy a variable to process, this should be the filename not including the extension
# - used to test whether the outputting is working correctly 
var       <- NULL

# set variable for global attributes in netcdf files
# person, email, institution
#person      <- 'Anthony P. Walker'
#email       <- 'walkerap@ornl.gov'
#institution <- 'Oak Ridge National Laboratory'
person      <- 'Patrick C. McGuire'
email       <- 'mcguirepatr@gmail.com'
institution <- 'University of Reading'

# project name
project   <- 'TRENDYv10, 2021'


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

if(!annual)  ncdf_avars <- NULL
if(!monthly) ncdf_mvars <- NULL

if(deg1){
  #nsites <- 15417
  #nsites <- 15729
  #nsites <- 15789
  lon    <- 1
  lat    <- 1  
}

if(is.null(outsyear)) outsyear <- sty
if(is.null(outeyear)) outeyear <- sty + ny - 1
print('')
print('')
print(sty)
print(ny)

### Read in functions 
##########################
setwd(fd)
source('read_SDGVM_output_pft.R')
if(netcdf) source('functions_netcdf.R')



### Start Program
##########################
#wd_list <- paste(wd,sim,'output_6_7/',sep='/')
#wd_list <- paste(wd,sim,'output_6_7_500/',sep='/')
wd_list <- paste(root_dir,sim,output_dir,sep='/')
print(wd_list)

# if single variable specified process that and nothing else
if(!is.null(var)) {
  monthly <- daily <- annual <- F
  monthly <- grepl('monthly',var) 
  daily   <- grepl('daily',var)
  annual  <- !(daily|monthly)
  pft     <- grepl('[A-Z]',var)
  wd      <- wd_list[1]  
  ifile   <- paste(var,'.dat',sep='')

  print(c(annual,monthly,daily,pft))
  if(annual)       lapply(ifile, stitch_annual,    grids, wd )
  if(monthly&!pft) lapply(ifile, stitch_subannual, grids, wd, atr=12,  ad=2, styear=sty, nyears=ny )
  if(daily&!pft)   lapply(ifile, stitch_subannual, grids, wd, atr=360, ad=1, styear=sty, nyears=ny )
  if(monthly&pft)  lapply(ifile, stitch_subannual, grids, wd, atr=12,  ad=1, styear=sty, nyears=ny )
  if(daily&pft)    lapply(ifile, stitch_subannual, grids, wd, atr=360, ad=1, styear=sty, nyears=ny )
 
  stich <- F
}


# stich mp data
if(stich) {
  lapply(wd_list[pia], stich_sdgvm_mp_apply,
         grids=grids, mc=mc,
         annual=annual, monthly=monthly, daily=daily,
         mc.cores=cores, styear=sty, nyears=ny )

}
print("finished stich_sdgvm_mp_apply")


# remove grid output now that combined files have been created
if(delete) {
  for( wdc in wd_list[pia]) {
    setwd(wdc)
    system("for i in grid*; do mv $i/simulation.dat $i/simulation.txt; done")
    system("for i in grid*; do mv $i/site_info.dat  $i/site_info.txt; done")
    system("for i in grid*; do mv $i/diag.dat       $i/diag.txt; done")
    system('for i in grid*; do for f in init*.dat; do mv $i/$f $i/"`basename "$f" .dat`.txt"; done; done')
    system("rm ./grid*/*.dat")
    system('for i in grid*; do for f in init*.txt; do mv $i/$f $i/"`basename "$f" .txt`.dat"; done; done')
}} 


# convert data to CMOR netcdf output
outnyears <- outeyear - outsyear + 1
if(netcdf) lapply(wd_list[pia], write_sdgvm_netcdf,
                  afiles=ncdf_avars,
                  mfiles=ncdf_mvars,
                  mc=mc, procs=cores,
                  nsites=nsites, nyears=outnyears, styr=sty, lon=lon, lat=lat,
                  osyr=outsyear, person=person, email=email, institution=institution )



### END ###  
