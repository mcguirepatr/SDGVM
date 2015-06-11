#########################
#
# Convert SDGVM ASCII to netcdf
#
# AWalker 
# Feb 2015
#
#########################

library(ncdf4)
library(lattice)
source('params_netcdf.R')
  
### Creat a netcdf from scratch
# wd     <- wd_list[pia]
# mfiles <- mfiles
# mc     <- F
write_sdgvm_netcdf <- function(wd,afiles=NULL,mfiles=NULL,dfiles=NULL,
                               mc=F,procs=4,...){
  
  # Function to concatenate SDGVM output when run across multiple cores
  ########################
  
  # annual files just open and stich together rowise
  # monthly open, convert structure and stich rowise
  # daily open, convert structure and stich rowise
  
  
  # Initialisation
  ########################
  setwd(wd)
  if(grepl('T1',wd)) fref <- 'T1'
  if(grepl('T2',wd)) fref <- 'T2'
  if(grepl('T3',wd)) fref <- 'T3'
  
  #read output file names
  annual <- monthly <- daily <- F
  if(!is.null(dfiles)) daily   <- T
  if(!is.null(mfiles)) monthly <- T
  if(!is.null(afiles)) annual  <- T
           
  
  # data processing
  ########################
  if(mc){
    if(annual)  mclapply(afiles,make_netcdf_TERRABITES,mc.cores=procs,annual=annual,fref=fref,...)
    if(monthly) mclapply(mfiles,make_netcdf_TERRABITES,mc.cores=procs,monthly=monthly,fref=fref,...)
    if(daily)   mclapply(dfiles,make_netcdf_TERRABITES,mc.cores=procs,daily=daily,fref=fref,...)    
  } else {
    if(annual)  lapply(afiles,make_netcdf_TERRABITES,annual=annual,fref=fref,...)
    if(monthly) lapply(mfiles,make_netcdf_TERRABITES,monthly=monthly,fref=fref,...)
    if(daily)   lapply(dfiles,make_netcdf_TERRABITES,daily=daily,fref=fref,...)        
  }
}


make_netcdf_TERRABITES <- function(var,annual=F,monthly=F,daily=F,...){
  # opens data and calls the make netcdf function
  
  # open data
  fname <- paste(var$file,'.dat',sep='')
  if(monthly) fname <- paste('monthly_',var$file,'.dat',sep='')
  if(daily)   fname <- paste('daily_',var$file,'.dat',sep='')
  cs <- 4
  if(daily) cs <- 3
  
  if(length(var$file)==1){
    df <- read.table(fname)    
  } else {  
    # if data need to be combined to create a composit variable
    for(file in fname){
      df1 <- read.table(file)
      if(file==fname[1]) {
        df <- df1
        rm(df1)
        # the below line sums the data in the multiple files - could make more flexible
      } else {
        df[,cs:length(df)] <- df[,cs:length(df)] + df1[,cs:length(df1)]        
        rm(df1)
      } 
    }
    # need to write a check in here
  }
  
  # scale variable to correct output units
  df[,cs:length(df)] <- df[,cs:length(df)] * var$scale
  
  # call function to make netcdf file for dataset
  make_netcdf(df,var=var$name,fname=var$name,lname=var$lname,units=var$units,
              annual=annual,monthly=monthly,daily=daily,cs=cs,...)
  rm(df)
}


make_netcdf <- function( df,var,fname=var,lname,units,
                         nsites=1548,nyears=110,lon=3.75,lat=2.5,mv=-99999,
                         annual=F,monthly=F,daily=F,cs=3,fref='',... ) {
  # expects df in the mp processed SDGVM output format
  
  # initialise - leap years not yet considered
  if(annual){
    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
  } else if(monthly){
    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
  } else if(daily){
    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
  }
  time_seq <- seq(st,end,st)
  
  # create the nc dimensions
  nclon   <- ncdim_def( name='lon',units='degrees_east',vals=(seq(lon/2,360-lon/2,lon)) )
  nclat   <- ncdim_def( name='lat',units='degrees_north',vals=(seq(lat/2,180-lat/2,lat)) )
  nctime  <- ncdim_def( name='time',units='hours since 1901-01-01 00:00:00',vals=time_seq,unlim=T )
  
  # create the ncdf4 object of the var
  ncvar <- ncvar_def( name=var,units=units,dim=list(nclon,nclat,nctime),missval=mv,longname=lname )
  
  # create the file given the var
  newnc <- nc_create( paste('SDGVM_',fref,'_',fname,'.nc',sep=''),ncvar )
  
  # put the data into the ncfile
  df$lats  <- (df[,1] + 90)/lat  + 1
  df$lons  <- (df[,2] + 180)/lon + 1
  da       <- array(1:(360/lon*180/lat*nyears*sa),dim=c(360/lon,180/lat,nyears*sa))
  da[]     <- NA
  time_sub <- unlist(lapply(1:nyears,slice,l=sa,nsites=nsites))
  smat     <- cbind(rep(df$lons,nyears),rep(df$lats,nyears),time_sub)
  da[smat] <- as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
  rm(df)
  rm(smat)
  rm(time_sub)
  ncvar_put(newnc,var,da)
  print(levelplot(da[,,101],main=var))
  rm(da)
  
  # set global attributes
  ncatt_put(newnc,0,attname='Conventions',attval='CF-1.4 (or close)')
  ncatt_put(newnc,0,attname='Calendar',attval='no leap years, 360 day years')
  ncatt_put(newnc,0,attname='institution',attval='Oak Ridge National Laboratory')
  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by: Anthony Walker (walkerap@ornl.gov)'))
  
  # write the file
  print(newnc)
  nc_close(newnc)  
}
slice <- function(i,l,nsites){
  s <- (i-1)*l + 1
  e <- i*l
  rep(s:e,nsites)
}


start_netcdf <- function( var,fname=var,lname,units,nsites=1548,nyears=110,lon=3.75,lat=2.5,mv=-99999,
                          annual=F,monthly=F,daily=F,fref='',... ){
  # expects df in the mp processed SDGVM output format
  
  # initialise - leap years not yet considered
  if(annual){
    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
  } else if(monthly){
    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
  } else if(daily){
    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
  }
  time_seq <- seq(st,end,st)
  
  # create the nc dimensions
  nclon   <- ncdim_def( name='lon',units='degrees_east',vals=(seq(lon/2,360-lon/2,lon)) )
  nclat   <- ncdim_def( name='lat',units='degrees_north',vals=(seq(lat/2,180-lat/2,lat)) )
  nctime  <- ncdim_def( name='time',units='hours since 1901-01-01 00:00:00',vals=time_seq,unlim=T )
  
  # create the file given the var
  newnc <- nc_create( paste('SDGVM_',fref,'_',fname,'.nc',sep=''),ncvar )
  
  # create the file given the var
  newnc <- nc_create( paste(fname,'.nc',sep=''),ncvar )
  print(newnc)
    
  # set global attributes
  ncatt_put(newnc,0,attname='Conventions',attval='CF-1.4 (or close)')
  ncatt_put(newnc,0,attname='Calendar',attval='no leap years, 360 day years')
  ncatt_put(newnc,0,attname='institution',attval='Oak Ridge National Laboratory')
  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by: Anthony Walker (walkerap@ornl.gov)'))
  
  # write the file
  nc_close(newnc)  
}

# puts data into a pre-existing netcdf
# could be modified to add data bit by bit for daily data
data_netcdf <- function( df,var,fname=var,lname,units,nsites=1548,nyears=110,lon=3.75,lat=2.5,mv=-99999,
                         annual=F,monthly=F,daily=F,cs=3 ){
  # expects df in the mp processed SDGVM output format
  
  # initialise - leap years not yet considered
  if(annual){
    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
  } else if(monthly){
    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
  } else if(daily){
    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
  }
  
  # open the file given the var
  newnc <- nc_open( paste(fname,'.nc',sep=''),write=T )
  print(newnc)
  
  ### put the data into the ncfile - this needs modification to allow multiple subsets of the data to be entered into the netcdf
  # create global array
  df$lats  <- (df[,1] + 90)/lat  + 1
  df$lons  <- (df[,2] + 180)/lon + 1
  da       <- array(1:(360/lon*180/lat*nyears*sa),dim=c(360/lon,180/lat,nyears*sa))
  da[]     <- NA
  # get subscripts
  time_sub <- unlist(lapply(1:nyears,slice,l=sa,nsites=nsites))
  smat     <- cbind(rep(df$lons,nyears),rep(df$lats,nyears),time_sub)
  # put data into matrix
  da[smat] <- as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
  rm(df)
  rm(smat)
  rm(time_sub)
  # put data into netcdf
  ncvar_put(newnc,var,da)
  print(levelplot(da[,,101],main=var))
  rm(da)
  
  # write the file
  nc_close(newnc)  
}
