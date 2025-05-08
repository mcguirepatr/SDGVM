##############################
#
# Open CRUJRA netcdf files & process for SDGVM 
#
# AWalker
# July 2022
#
##############################

### CRU-JRA files from ...
#
#   see protocol for link
#
#   6 hourly data 1901 - 2021
#   0.5 x 0.5 degree resolution
#   no leap years
#
#   CRU units
#   tmp    - K
#   Q      - g/g
#   rain   - mm/6hrs
#   swdown - J/m2/s or W/m2

# FYI - the new radiation files have changed to units from J/m/step to J/m/s
# 
# Cheers,
# Luke

### SDGVM output units and format
#
#   tmp - monthly mean - oC * 10; 4 column integers, no spaces
#   hum - monthly mean - %  * 10; 4 column integers, no spaces
#   prc - monthly sum  - mm     ; 4 column integers, no spaces
#   swr - monthly mean of daily 24 hr mean  - W/m2?? ; 6 column real nos 6.2, 1 space

# Outputs varying site faster than year, then runs a script which reads output and reorganises by year varying fastest for SDGVM 


rm(list=ls())

### Functions
#######################################

wd <- '/mnt/disk2/script_library/R/map_tools'
setwd(wd)
source('regrid.R')



mat_SDGVM_write <- function(v, ncol=12, w=4, ofile, r=0, s=' ' ) { 
  # designed for use with apply
  # appends a vector to the output file unless vector is all NAs
  # written in near SDGVM format, only difference is columns are separated by a blank space
  if(sum(is.na(v))!=length(v))  write(format(round(v,r),width=w), ofile, append=T, ncolumns=ncol, sep=s )
}

SDGVM_write <- function(df, ofile, w=4, r=0, s='' ) {
  # writes a matrix in SDGVM input format
  write.table(format(round(df,r),width=w), ofile, row.names=F, col.names=F, sep=s, quote=F )  
}

conv_sh_to_rh <- function(q, t, p=101325 ) {
  # convert specific humidity to relative
  # use Buck 1981 for sat vp (as does SDGVM)
  
  # q - water:dry air mass ratio, g/g
  # t - temperature,              K
  # p - pressure,                 Pa
  
  paw  <- p / (0.62198/q + 1)
  paws <- 6.1121 * exp(17.502*(t-273.15) / (240.97 + t - 273.15)) * 100
  #print(c(paw,paws))
  paw/paws * 1000      # units for SDGVM 0.1 % or parts per mil
}

# conversion functions for CRU-NCEP to SDGVM
conv_K_to_oC <- function(t,dummy)    (t - 273.15) * 10       # units 0.1 oC
conv_prc     <- function(prc,dummy)  prc                     # mm
# conv_sw      <- function(swr,d_in_m) swr / (d_in_m*24*3600)  # W/m2
conv_sw      <- function(swr, dummy ) swr                    # W/m2
conv_fd      <- function(fd, dummy )  fd                     # unitless, proportion



### Initialise
#######################################

library(lattice)
library(grid)
library(parallel)
library(ncdf4)
library(gdata)

# paths
src  <- '/mnt/disk2/script_library/R/map_tools'
date <- Sys.Date()
dir  <- '/mnt/disk3/'
wd   <- paste(dir,'/Databases/CRU/CRU-JRA/2022/',sep='')
wdd  <- paste(wd,'/trendy/',sep='')
wdo  <- paste(wd,'/SDGVMdata/',sep='')
fdir <- wd 
# setwd(wdd)

# filename variables
fname   <- 'crujra.v2.3.5d'
vars    <- c('tswrf','tmp','spfh','pre','fd')
svars   <- c('swr','tmp','hum','prc','fdf')
ncvars  <- vars
time    <- 1901:2021
fsuf    <- '365d.noc'
fend    <- 'nc'
zend    <- 'gz'
send    <- '.dat'
lon_res <- 720
lat_res <- 360
  
# time variables
year      <- 1460
nyears    <- 1
days_in_m <- c(31,28,31,30,31,30,31,31,30,31,30,31)

write_mask <- F


### Rough working
#######################################
v <- 4
t <- 1

setwd(paste0(wd,'/',vars[v]))
ifile    <- paste(fname,vars[v],time[t],sep='.')
system(paste('gunzip -c',paste(ifile,fend,zend,sep='.'),' > ',paste(ifile,fend,sep='')))
mynetcdf <- nc_open(paste(ifile,fend,sep=''))

# xm is a 3 dimensional array - by lon, lat and time (in 6 hrly intervals)
x3d  <- ncvar_get(mynetcdf, ncvars[v], start=c(1,1,500), count=c(lon_res,lat_res,4) )
x3d[which(x3d< -100)] <- NA

# x is a 2-d matrix by lon and lat averaged across the time dimension of xm
x2d  <- apply(x3d,c(1,2),sum,na.rm=F)
summary(as.vector(x2d))

levelplot(x2d,
          col.regions=rev(heat.colors(50)),contour=F,
          main=paste(vars[v],time[t]),ylab='',xlab='',
          panel=function(...){
            panel.fill(col='black')
            panel.levelplot(...)
          })

# close file
nc_close(mynetcdf)

# remove unzipped file
system(paste('rm',paste(ifile,fend,sep='.')))



# create SDGVM maskmap & sites list 
################################

# load regridding functions
setwd(src)
source('regrid.R')

# write maskmap and lat lon file
if(write_mask) {
  
  # open netcdf file
  v <- 1
  t <- 1  
  setwd(paste0(wd,'/',vars[v]))
  ifile    <- paste(fname,vars[v],time[t],sep='.')
  system(paste('gunzip -c',paste(ifile,fend,zend,sep=''),' > ',paste(ifile,fend,sep='')))
  mynetcdf <- nc_open(paste(ifile,fend,sep=''))
  
  # remove unzipped file
  system(paste('rm',paste(ifile,fend,sep='')))
  
  # extract mask
  # mask <- ncvar_get(mynetcdf,'mask',start=c(1,1),count=c(lon_res,lat_res))
  x3d  <- ncvar_get(mynetcdf, ncvars[v], start=c(1,1,1), count=c(lon_res,lat_res,1460) )

  # x is a 2-d matrix by lon and lat averaged across the time dimension of xm
  mask  <- apply(x3d,c(1,2),mean,na.rm=T)
  summary(as.vector(mask))
  dim(mask)
  mask[which(!is.na(mask))] <- 1:( prod(dim(mask)) - length(which(is.na(mask))) )
  
  # mask <- mask + 1
  mask[,69]

  # plot mask native orientation  
  levelplot(mask,ylim=c(lat_res,1),col.regions=rev(heat.colors(50)),contour=F,
            panel=function(...){
              panel.fill(col='black')
              panel.levelplot(...)
            })

  # plot mask SDGVM orientation  
  mask              <- mask[,360:1]
  mask[is.na(mask)] <- 0  
  levelplot(mask,ylim=c(lat_res,1),col.regions=rev(heat.colors(50)),contour=F,
            panel=function(...){
              panel.fill(col='black')
              panel.levelplot(...)
            })
  
  # write this years data to output
  setwd(paste(wdo,'0.5deg',sep=''))
  write.table(format(t(mask),width=6),'maskmap.dat',sep='',quote=F,row.names=F,col.names=F)  
  
  # close file
  nc_close(mynetcdf)
  

  # generate sites list
  dim(mask)
  mask[mask==0] <- NA  
  
  setwd(paste(wdo,'0.5deg',sep=''))
  lat <- seq(-89.75,   89.75, 0.5 )
  lon <- seq(-179.75, 179.75, 0.5 )
  length(lat)
  length(lon)
  
  for(ns in 1:360){
    for(we in 1:720){
      if(!is.na(mask[we,ns])) write(format(round(c(lat[ns],lon[we]),3),width=8,nsmall=3), 'lat_lon_unformatted.dat', append=T, ncolumns=2, sep=' ' )
    }
  }
  
  df1 <- read.table('lat_lon_unformatted.dat')
  head(df1)
  
  library(gdata)
  write.fwf(df1, width=c(7,9), nsmall=3, 'lat_lon_formatted.dat', colnames=F, sep='' )
  
  
  # regrid to 1 degree
  mask_1deg <- reduce_res(mask)
  dim(mask_1deg)
  which(is.na(mask_1deg))
  
  lat <- seq(-89.5,   89.5, 1 )
  lon <- seq(-179.5, 179.5, 1 )
  length(lat)
  length(lon)
  
  for(ns in 1:180){
    for(we in 1:360){
      if(!is.na(x2d_1deg[we,ns])) write(format(round(c(lat[ns],lon[we]),3),width=8,nsmall=3), 'lat_lon_unformatted_1deg.dat', append=T, ncolumns=2, sep=' ' )
    }
  }
  
  df1 <- read.table('lat_lon_unformatted_1deg.dat')
  head(df1)
  
  library(gdata)
  write.fwf(df1, width=c(7,9), nsmall=3, 'lat_lon_formatted_1deg.dat', colnames=F, sep='' )
  
}



### Program proper
####################################
# read and process data for monthly output

# miss out tair & qair
via <- c(1,4,5)
via <- c(4,5)
via <- 1
vars[via]
tia <- 1:length(time)

# variable loop
v <- 1
y <- 1
for(v in via) {
  if(vars[v]=='fd')    { func <- mean ; conv <- conv_sw       ; r <- 2 ; w <- 4 }
  if(vars[v]=='tswrf') { func <- mean ; conv <- conv_sw       ; r <- 2 ; w <- 7 }
  if(vars[v]=='tmp')   { func <- mean ; conv <- conv_K_to_oC  ; r <- 0 ; w <- 4 }
  if(vars[v]=='spfh')  { func <- mean ; conv <- conv_sh_to_rh ; r <- 0 ; w <- 4 }
  if(vars[v]=='pre')   { func <- sum  ; conv <- conv_prc      ; r <- 0 ; w <- 4 }

  # dummy arrays
  month <- array(NA, dim=c(lon_res,lat_res,12) )
  x3d   <- array(NA, dim=c(lon_res,lat_res,max(days_in_m)*4) )
  x2d   <- x3d[,,1]  
    
  # working dirs
  setwd(paste0(wd,'/',vars[v]))
  write(paste(vars[v],'started',date()),'out.txt',append=T)
  
  # annual loop
  for(y in tia){
    
    # open netcdf file
    setwd(paste0(wd,'/',vars[v]))
    ifile <- 
      if((vars[v]=='tswrf') | (vars[v]=='fd')) {
        paste(vars[v],'v11',time[y],sep='_')
      } else {
        paste(fname,vars[v],time[y],fsuf,sep='.')
      }
    system(paste('gunzip -c',paste(ifile,fend,zend,sep='.'),' > ',paste(ifile,fend,sep='.')))
    mynetcdf <- nc_open(paste(ifile,fend,sep='.'))
    
    # month loop
    hoy <- 1
    #m <- 1
    for(m in 1:12) {
      
      # extract data by month
      x3d[,,1:(days_in_m[m]*4)] <- ncvar_get(mynetcdf, ncvars[v], start=c(1,1,hoy), count=c(lon_res[1],lat_res[1],days_in_m[m]*4) )
      
      # convert NaN
      if(vars[v]=='tswrf' | vars[v]=='fd') x3d[which(is.nan(x3d[,,1:(days_in_m[m]*4)]))] <- 0.0     
      
      # take monthly mean/total
      x2d[] <- apply(x3d[,,1:(days_in_m[m]*4)], c(1,2), func, na.rm=F )

      # add to monthly mean dataframe
      month[,,m] <- x2d
      
      # increment start day
      hoy   <- hoy + days_in_m[m]*4      
      
      # NA array
      x3d[] <- NA
    }

    # convert units to SDGVM units
    month[] <- conv(month, days_in_m[m] )

    # write this years data to output
    setwd(paste0(wdo,'monthly/0.5deg/'))
    apply(month, c(1,2), mat_SDGVM_write, ofile=paste(svars[v],'_byyear',send,sep=''), r=r, w=w )
    
    # close file
    nc_close(mynetcdf)
    
    # remove unzipped file
    setwd(paste0(wd,'/',vars[v]))
    system(paste('rm',paste(ifile,fend,sep='.')))
    write(paste(time[y],'completed'),'out.txt',append=T)
  }
  
  write(paste(vars[v],'completed',date()),'out.txt',append=T)
}



# Make plots
##########################

# for(m in 1:12){
#   print(
#   levelplot(month[,,m],ylim=c(lat_res[1],1),col.regions=rev(heat.colors(50)),contour=F,
#             panel=function(...){
#               panel.fill(col='black')
#               panel.levelplot(...)
#             })  
#   )
# }
# 
# levelplot(x3drh[,,1:4],ylim=c(lat_res[1],1),col.regions=rev(heat.colors(50)),contour=F,
#           panel=function(...){
#             panel.fill(col='black')
#             panel.levelplot(...)
#           })
# 
# levelplot(x2d,ylim=c(lat_res[1],1),col.regions=rev(heat.colors(50)),contour=F,
#           panel=function(...){
#             panel.fill(col='black')
#             panel.levelplot(...)
#           })
# 
# 
# levelplot(monthc[,,1],ylim=c(lat_res[1],1),col.regions=rev(heat.colors(50)),contour=F,
#                 panel=function(...){
#                   panel.fill(col='black')
#                   panel.levelplot(...)
#                 })
 
 

# Temperature and RH
##########################################

# setup arrays
montht <- month <- array(NA,dim=(c(lon_res,lat_res,12)))
x3drh <- x3dt <- x3d <- array(NA,dim=c(lon_res,lat_res,max(days_in_m)*4))
x2dt <- x2d <- x3d[,,1]  


# temperature and RH loop
tia <- 1:length(time)
vt  <- 2      # temperature
for(v in 3) { # specific humidity
  if(vars[v]=='spfh')   { func <- mean    ; conv <- conv_sh_to_rh ; r <- 0 ; w <- 4 }
  
  setwd(paste0(wd,'/',vars[v]))
  write(paste(vars[2:3],'started',date()),'out.txt',append=T)

  # annual loop
  for(y in tia){

    # open netcdf file
    ifile  <- paste(fname, vars[v],  time[y], sep='.' )
    ifilet <- paste(fname, vars[vt], time[y], sep='.' )

    # unzip files
    # open netcdf files
    setwd(paste0(wd,'/',vars[v]))
    system(paste('gunzip -c', paste(ifile,fsuf,fend,zend,sep='.'),  ' > ', paste(ifile,fsuf,fend,sep='.'))  )
    mynetcdf  <- nc_open(paste(ifile,fsuf,fend,sep='.'))
    
    setwd(paste0(wd,'/',vars[vt]))
    system(paste('gunzip -c', paste(ifilet,fsuf,fend,zend,sep='.'), ' > ', paste(ifilet,fsuf,fend,sep='.')) )
    mynetcdft <- nc_open(paste(ifilet,fsuf,fend,sep='.'))

    # month loop
    hoy <- 1
    m   <- 1
    for(m in 1:12){

      # extract data by month and convert sh to rh at each time point i.e. pre-averaging
      setwd(paste0(wd,'/',vars[v]))
      x3d[,,1:(days_in_m[m]*4)]  <- ncvar_get(mynetcdf, ncvars[v], start=c(1,1,hoy), count=c(lon_res[1],lat_res[1],days_in_m[m]*4) )
      setwd(paste0(wd,'/',vars[vt]))
      x3dt[,,1:(days_in_m[m]*4)] <- ncvar_get(mynetcdft, ncvars[vt], start=c(1,1,hoy), count=c(lon_res[1],lat_res[1],days_in_m[m]*4) )

      # # convert missing values
      # x3d[which(x3d< -1000)]   <- NA
      # x3dt[which(x3dt< -1000)] <- NA

      # convert to RH
      x3drh[] <- conv_sh_to_rh(x3d,x3dt)
      x3drh[x3drh>1000] <- 1000

      # take monthly mean
      x2d[]   <- apply(x3drh[,,1:(days_in_m[m]*4)], c(1,2), func, na.rm=F )
      x2dt[]  <- apply(x3dt[,,1:(days_in_m[m]*4)],  c(1,2), func, na.rm=F )

      # add to monthly mean dataframe
      month[,,m]  <- x2d
      montht[,,m] <- x2dt
      # increment start day
      # print(hoy)
      hoy <- hoy + days_in_m[m]*4
    }

    # convert units to SDGVM units
    montht[] <- conv_K_to_oC(montht)
    
    # write this years data to output
    setwd(paste0(wdo,'monthly/0.5deg/'))
    apply(month,c(1,2),  mat_SDGVM_write, ofile=paste(svars[v], '_byyear',send,sep=''), r=r, w=w )
    apply(montht,c(1,2), mat_SDGVM_write, ofile=paste(svars[vt],'_byyear',send,sep=''), r=r, w=w )

    # close files
    nc_close(mynetcdf)
    nc_close(mynetcdft)

    # remove unzipped files
    setwd(paste0(wd,'/',vars[vt]))
    system(paste('rm',paste(ifilet,fsuf,fend,sep='.')))
    setwd(paste0(wd,'/',vars[v]))
    system(paste('rm',paste(ifile,fsuf,fend,sep='.')))
    write(paste(time[y],'completed'),'out.txt',append=T)
  }
  write(paste(vars[2:3],'completed',date()),'out.txt',append=T)
}



# check written output
##################################


# setwd(paste0(wdo,'monthly/0.5deg/'))
# testdf <- as.matrix(read.table(paste(svars[v],'_byyear',send,sep='')))
# 
# class(testdf)
# dim(testdf)
# dim(testdf)[1]/length(time)

# m1 <- matrix(mprc[,1],nrow=720)
# levelplot(m1,ylim=c(lat_res,1),col.regions=rev(heat.colors(50)),contour=F,
#           panel=function(...){
#             panel.fill(col='black')
#             panel.levelplot(...)
#           })
# matrix(1:20,nrow=4)



# Read files organised by year and write files organised by site
###############################
years <- length(time)
setwd(paste0(wdo,'monthly/0.5deg/'))

via <- 1:length(vars)
for( v in via) {
  # read by year file
  ifile      <- paste(svars[v],'_byyear',send,sep='')
  databyyear <- as.matrix(read.table(ifile))
  # set up subscripts 
  l          <- length(databyyear[,1])
  sites      <- l/years
  ss         <- as.vector(apply(as.matrix(1:sites),1,seq,to=l,by=sites))
  # re-arrange and write data by site
  if(vars[v]=='tswrf') { r <- 2 ; w <- 7 }
  if(vars[v]=='fd')    { r <- 2 ; w <- 4 }
  if(vars[v]=='tmp')   { r <- 0 ; w <- 4 }
  if(vars[v]=='spfh')  { r <- 0 ; w <- 4 }
  if(vars[v]=='pre')   { r <- 0 ; w <- 4
    # cap precip at 9999 mm mon-1
    if(any(databyyear>9999)) {
      nms  <- length(which(databyyear>9999))
      over <- sum(databyyear[which(databyyear>9999)] - 9999)
      write(paste('prc sites and months over 9999: n,',nms,'; total over (mm),',over),'out.txt',append=T)
      databyyear[which(databyyear>9999)] <- 9999
    }
  }
  ofile      <- paste(svars[v],send,sep='') 
  SDGVM_write(databyyear[ss,],ofile,r=r,w=w)
}



### END ###