#########################
#
# Convert SDGVM ASCII to netcdf
#
# AWalker 
# Feb 2015
#
#########################

library(parallel)
library(ncdf4)
library(lattice)

source('params_netcdf.R')
  
# wd     <- wd_list[pia]
# mfiles <- mfiles
# mc     <- F
write_sdgvm_netcdf <- function(wd,afiles=NULL,mfiles=NULL,dfiles=NULL,
                               mc=F,procs=4,...){
  
  # Function to take concatenated SDGVM output and process to TRENDY netcdf output format
  ########################
  
  
  # Initialisation
  ########################
  setwd(wd)
 
  # scan the directory name for the simlation id, passed thru functions to generate output filename  
  if(grepl('T1',wd)) fref   <- 'T1'
  if(grepl('T2',wd)) fref   <- 'T2'
  if(grepl('T3',wd)) fref   <- 'T3'
  if(grepl('Stest',wd)) fref <- 'Stest'
  if(grepl('S0',wd)) fref   <- 'S0'
  if(grepl('S1',wd)) fref   <- 'S1'
  if(grepl('S2',wd)) fref   <- 'S2'
  if(grepl('S3',wd)) fref   <- 'S3'
  if(grepl('S4',wd)) fref   <- 'S4'
  if(grepl('tpu0',wd)) fref <- 'tpu0'
  if(grepl('tpu1',wd)) fref <- 'tpu1'
  if(grepl('tpu2',wd)) fref <- 'tpu2'
  if(grepl('tpu3',wd)) fref <- 'tpu3'

  print('',quote=F)
  print('',quote=F)
  print(paste('Processing simulation:',fref),quote=F)
  print('',quote=F)

  
  # determine which temporal and spatial scales to process  
  annual <- monthly <- daily <- F
  if(!is.null(dfiles)) daily   <- T
  if(!is.null(mfiles)) monthly <- T
  if(!is.null(afiles)) annual  <- T
           
  
  # data processing - call make_netcdf_TERRABITES function
  ########################

  # parrallel process or not
#  if(mc){
#    if(annual)  mclapply(afiles,make_netcdf_TERRABITES,mc.cores=procs,annual=annual,fref=fref,...)
#    if(monthly) mclapply(mfiles,make_netcdf_TERRABITES,mc.cores=procs,monthly=monthly,fref=fref,...)
#    if(daily)   mclapply(dfiles,make_netcdf_TERRABITES,mc.cores=procs,daily=daily,fref=fref,...)    
#  } else {
#    if(annual)  lapply(afiles,make_netcdf_TERRABITES,annual=annual,fref=fref,...)
#    if(monthly) lapply(mfiles,make_netcdf_TERRABITES,monthly=monthly,fref=fref,...)
#    if(daily)   lapply(dfiles,make_netcdf_TERRABITES,daily=daily,fref=fref,...)        
#  }

  # parrallel process or not
  if(mc){
    if(annual)  mclapply(afiles,make_netcdf_TRENDY,mc.cores=procs,annual=annual,fref=fref,...)
  } else {
    if(annual)  lapply(afiles,make_netcdf_TRENDY,annual=annual,fref=fref,...)
  }
  if(monthly) lapply(mfiles,make_netcdf_TRENDY,monthly=monthly,fref=fref,...)
  if(daily)   lapply(dfiles,make_netcdf_TRENDY,daily=daily,fref=fref,...)        

}


#make_netcdf_TERRABITES <- function(varo,annual=F,monthly=F,daily=F,pft=NULL,...){
#  # opens data and calls the make netcdf function
#  var <- get(varo)
#  
#  # open data
#  fname <- var$file
#  if(monthly) fname <- paste('monthly',var$file,sep='_')
#  if(daily)   fname <- paste('daily',var$file,sep='_')
#  cs <- 4
#  if(daily) cs <- 3
#  
#  # read data 
#  # if SDGVM output is the same variable as required
#  if(length(var$file)==1){
#  
#    # if by pft output required open each file ans stack them rowwise
#    if(is.null(pft)) df <- read.table( paste(fname,'dat',sep='.') )    
#    else { 
#      for( p in pft ) {
#        df1 <- read.table( paste(fname,'_',p,'.dat',sep='') )
#        df  <- if( p == pft[1] ) df1 else rbind(df,df1) 
#        rm(df1)
#      }
#    }
#
#  # if SDGVM output needs to be combined to create a composite variable
#  } else {  
#  
#    for(file in fname) {
#      # if by pft output required open each file ans stack them rowwise
#      if(is.null(pft)) df <- read.table( paste(fname,'dat',sep='.') )    
#      else { 
#        for( p in pft ) {
#          df1 <- read.table( paste(fname,'_',p,'.dat',sep='') )
#          df2 <- if( p == pft[1] ) df1 else rbind(df2,df1) 
#          rm(df1)
#        }
#      }
#  
#      # sum the data in the multiple files - could make more flexible
#      if(file==fname[1]) {
#        df <- df2
#        rm(df2)
#      } else {
#        df[,cs:length(df)] <- df[,cs:length(df)] + df2[,cs:length(df1)]        
#        rm(df2)
#      } 
#    }
#    # need to write a check in here
#  }
#  
#  # scale variable to correct output units
#  df[,cs:length(df)] <- df[,cs:length(df)] * var$scale
#  
#  # call function to make netcdf file for dataset
#  make_netcdf(df,var=var$name,fname=var$name,lname=var$lname,units=var$units,
#              pft=pft,annual=annual,monthly=monthly,daily=daily,cs=cs,addnotes=var$addnotes,...)
#  rm(df)
#}
#
#
#make_netcdf <- function( df,var,fname=var,lname,units,
#                         nsites=1548,nyears=110,styr=1901,lon=3.75,lat=2.5,mv=-99999,
#                         pft=NULL,annual=F,monthly=F,daily=F,cs=3,fref='',addnotes=NULL,... ) {
#  # expects df in the mp processed SDGVM output format
#  
#  # initialise - leap years not yet considered
#  if(annual){
#    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
#  } else if(monthly){
#    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
#  } else if(daily){
#    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
#  }
#  time_seq <- seq(st,end,st)
#  
#  # create the nc dimensions
#  nclon   <- ncdim_def( name='lon',units='degrees_east',vals=(seq(lon/2,360-lon/2,lon)) )
#  nclat   <- ncdim_def( name='lat',units='degrees_north',vals=(seq(lat/2,180-lat/2,lat)) )
#  nctime  <- ncdim_def( name='time',units=paste('hours since ',styr,'-01-01 00:00:00',sep=''),vals=time_seq,unlim=T )
#  if(!is.null(pft)) ncpft   <- ncdim_def( name='vegtype',units='pft id, see notes',vals=1:length(pft) )
#  
#  # create the ncdf4 object of the var
#  if(!is.null(pft)) ncvar <- ncvar_def( name=var,units=units,dim=list(nclon,nclat,ncpft,nctime),missval=mv,longname=lname )
#  else              ncvar <- ncvar_def( name=var,units=units,dim=list(nclon,nclat,nctime),missval=mv,longname=lname )
#  
#  # create the file given the var
#  newnc <- nc_create( paste('SDGVM_',fref,'_',fname,'.nc',sep='') , ncvar )
#  
#  # create an nc ready array 
#  df$lats  <- (df[,1] + 90)/lat  + 1
#  df$lons  <- (df[,2] + 180)/lon + 1
#  a_dim    <- if(!is.null(pft)) c(360/lon,180/lat,length(pft),nyears*sa) else c(360/lon,180/lat,nyears*sa)
#  da       <- array(NA,dim=a_dim)
#
#  # put the data into an nc ready array 
#  time_sub <- unlist(lapply(1:nyears,slice,l=sa,nsites=nsites))
#  smat     <- cbind(rep(df$lons,nyears),rep(df$lats,nyears),time_sub)
#  da[smat] <- as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
#  rm(df); rm(smat); rm(time_sub)
#  
#  # put the data into the ncfile
#  #if(!is.null(pft)) ncvar_put(newnc,var,da,start=(1,1,which(pft==p),1)
#  #else     
#  ncvar_put(newnc,var,da)
#  #print(levelplot(da[,,101],main=var))
#  rm(da)
#  
#  # set global attributes
#  ncatt_put(newnc,0,attname='Conventions',attval='CF-1.4 (or close)')
#  ncatt_put(newnc,0,attname='Calendar',attval='no leap years, 360 day years')
#  ncatt_put(newnc,0,attname='institution',attval='Oak Ridge National Laboratory')
#  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by: Anthony Walker (walkerap@ornl.gov)'))
#  if(!is.null(notes)) {
#    trash <- T
#  } 
#
#  # write the file
#  print(newnc)
#  nc_close(newnc)  
#}


slice <- function(i,l,nsites){
  s <- (i-1)*l + 1
  e <- i*l
  rep(s:e,nsites)
}



#make_netcdf_TERRABITES <- function(varo,annual=F,monthly=F,daily=F,pft=NULL,...){
#make_netcdf_TRENDY <- function( df,var,fname=var,lname,units,
#                         nsites=1548,nyears=110,styr=1901,lon=3.75,lat=2.5,mv=-99999,
#                         pft=NULL,annual=F,monthly=F,daily=F,cs=3,fref='',addnotes=NULL,... ) {
make_netcdf_TRENDY <- function(varo,annual=F,monthly=F,daily=F,fref='',
                               nsites=1548,nyears=110,styr=1901,lon=3.75,lat=2.5,mv=-99999,
                               ... ) {

  # creates netcdf and call the read write function
  var <- get(varo)
  pft <- if(var$pft) pftnames else NULL 
 
  print('',quote=F)
  print('',quote=F)
  print(paste('Processing variable:',var$name),quote=F)
  print('Processing PFTs:',quote=F)
  print(pft,quote=F)
  print('',quote=F)

 
  # initialise time dimension - leap years not yet considered
  # annmonth is for variables that are annual output from SDGVM but are required in monthly format for TRENDY
  if(annual&is.null(var$annmonth)){
    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
  } else if(monthly|!is.null(var$annmonth)){
    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
  } else if(daily){
    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
  }
  # provides time at mid-point of timestep 
  time_seq <- seq(st,end,st) - st/2
  

  # create the nc dimensions
  nclon   <- ncdim_def( name='lon',units='degrees_east',vals=(seq(lon/2,360-lon/2,lon)) )
  nclat   <- ncdim_def( name='lat',units='degrees_north',vals=(seq(-90+lat/2,90-lat/2,lat)) )
  nctime  <- ncdim_def( name='time',units=paste('hours since ',styr,'-01-01 00:00:00',sep=''),vals=time_seq,unlim=T )
  if(!is.null(pft)) ncpft <- ncdim_def( name='vegtype',units='pft id, see notes',vals=1:length(pft) )
  

  # create the ncdf4 object of the var
  if(!is.null(pft)) ncvar <- ncvar_def( name=var$name,units=var$units,dim=list(nclon,nclat,ncpft,nctime),missval=mv,longname=var$lname )
  else              ncvar <- ncvar_def( name=var$name,units=var$units,dim=list(nclon,nclat,nctime),missval=mv,longname=var$lname )
  

  # create the nc file given the var
  newnc <- nc_create( paste('SDGVM_',fref,'_',var$name,'.nc',sep='') , ncvar )
  

  # set attributes
  ncatt_put(newnc,'time',attname='calendar',attval='360_day')


  # set global attributes
  ncatt_put(newnc,0,attname='title',attval='SDGVM output for TRENDYv7, 2018')
  ncatt_put(newnc,0,attname='Conventions',attval='CF-1.4 (or close)')
  ncatt_put(newnc,0,attname='calendar',attval='no leap years, 360 day years')
  ncatt_put(newnc,0,attname='institution',attval='Oak Ridge National Laboratory')
  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by: Anthony Walker (walkerap@ornl.gov)'))
  if(!is.null(var$notes)) ncatt_put(newnc,0,attname='notes',attval=var$notes)
  if(!is.null(pft))       ncatt_put(newnc,0,attname='PFTs',attval=paste(paste(pft,collapse=' '),'. These PFT distributions were derived by combining the HYDE 3.2 land-use and land-cover change database with the ESA GLCP 2014 data categorised according to SDGVM PFTs.',sep='') )
   

  # data filename 
  fname <- var$file
  if(monthly) fname <- paste('monthly',var$file,sep='_')
  if(daily)   fname <- paste('daily',var$file,sep='_')
  cs <- 3
  if(monthly) cs <- 4
  #cs <- 4
  #if(daily)   cs <- 3
  

  # data reading & writing loop
  # if by pft output required 
  if(is.null(pft)) {
    fnamefull <- paste(fname,'dat',sep='.')  
    newnc     <- readSDGVM_writeNCDF(fnamefull,var,newnc,ncvar,cs,NULL,
                                     lat,lon,nyears,nsites,sa,mv) 
  } else { 
    for( p in pft ) {
      fnamefull <- paste(fname,'_',p,'.dat',sep='')
      newnc     <- readSDGVM_writeNCDF(fnamefull,var,newnc,ncvar,cs,which(pft==p), 
                                       lat,lon,nyears,nsites,sa,mv) 
    }
  }


  # write the file
  print(newnc)
  nc_close(newnc)  
}


readSDGVM_writeNCDF <- function(fname,var,newnc,ncvar,cs,pftid=NULL,
                                lat,lon,nyears,nsites,sa,mv) {

  # if SDGVM output is the same variable as required
  if(length(var$file)==1) {
    df <- read.table( fname )    
  
  # if SDGVM output needs to be combined to create a composite variable
  } else {  

    for(file in fname) {
      df1 <- read.table( file )    

      # sum the data in the multiple files - could make more flexible
      if(file==fname[1]) {
        df <- df1
        rm(df1)
      } else {
        df[,cs:length(df)] <- df[,cs:length(df)] + df1[,cs:length(df1)]        
        rm(df1)
      } 
    }
    # need to write a check in here
  }

  print('Read file(s):',quote=F)
  print(fname,quote=F)
  print(head(df),quote=F)
  print('',quote=F)

  # for monthly grid square NBP remove fire flux, lulcc flux, and leached DOC flux
  if(fname=='monthly_nep.dat') {

    #print('Made it to routine to add C fluxes to nep')

    m1 <- as.matrix(read.table('lch.dat'))
    m2 <- as.matrix(read.table('lulccc.dat'))
    m3 <- as.matrix(read.table('fcn.dat'))

    #print('')
    #print(df[1:12,cs:length(df)])  
    #print('')
    #print(apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] )  
    #print('')
    #print(apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] )  
    #print('')
    #print(apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] ) 

    #print('')
    #print(sum(df[,cs:length(df)]))  
    #print(sum(apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 )))  
    #print(sum(apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 )))  
    #print(sum(apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 ))) 

    df[,cs:length(df)] <- df[,cs:length(df)] - 
				apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 ) 
    rm(m1,m2,m3)
  } 



  # scale variable to correct output units
  if( var$name %in% c('tas') )  df[,cs:length(df)] <- df[,cs:length(df)] + var$scale
  else                          df[,cs:length(df)] <- df[,cs:length(df)] * var$scale

  # rename lats and lons
  df$lats  <- (df[,1] + 90)/lat  + 1
  df$lons  <- (df[,2] + 180)/lon + 1

  # put the data into an nc ready array
  # if(!is.null(var$annmonth)) then need to modify time_sub to index the first month of the year only   
  time_sub <- unlist(lapply(1:nyears,slice,l=sa,nsites=nsites))
  #smat     <- cbind(rep(df$lons,nyears),rep(df$lats,nyears),time_sub)
  smat     <- if(!is.null(pftid)) cbind(rep(df$lons,nyears),rep(df$lats,nyears),1,time_sub) else 
                                  cbind(rep(df$lons,nyears),rep(df$lats,nyears),time_sub)
  # if(!is.null(pftid)) counts <- c(360/lon,length(df$lats),1,nyears*sa) 
  rm(time_sub)
  
  # create an nc ready array 
  #a_dim    <- c(360/lon,180/lat,nyears*sa)
  a_dim    <- if(!is.null(pftid)) c(360/lon,180/lat,1,nyears*sa) else 
                                  c(360/lon,180/lat,nyears*sa)
  da       <- array(mv,dim=a_dim)

  # put the data into an nc ready array 
  print(c(cs,nyears))
  as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
  da[smat] 
  da[smat] <- as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
  rm(df); rm(smat)
  
  # put the data into the ncfile
  if(!is.null(pftid)) ncvar_put( newnc,ncvar,da,start=c(1,1,pftid,1),count=a_dim )
  else                ncvar_put( newnc,ncvar,da )
  #print(levelplot(da[,,101],main=var))
  rm(da)

  newnc
}


#start_netcdf <- function( var,fname=var,lname,units,nsites=1548,nyears=110,lon=3.75,lat=2.5,mv=-99999,
#                          annual=F,monthly=F,daily=F,fref='',... ){
#  # expects df in the mp processed SDGVM output format
#  
#  # initialise - leap years not yet considered
#  if(annual){
#    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
#  } else if(monthly){
#    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
#  } else if(daily){
#    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
#  }
#  time_seq <- seq(st,end,st)
#  
#  # create the nc dimensions
#  nclon   <- ncdim_def( name='lon',units='degrees_east',vals=(seq(lon/2,360-lon/2,lon)) )
#  nclat   <- ncdim_def( name='lat',units='degrees_north',vals=(seq(lat/2,180-lat/2,lat)) )
#  nctime  <- ncdim_def( name='time',units='hours since 1901-01-01 00:00:00',vals=time_seq,unlim=T )
#  
#  # create the file given the var
#  newnc <- nc_create( paste('SDGVM_',fref,'_',fname,'.nc',sep=''),ncvar )
#  
#  # create the file given the var
#  newnc <- nc_create( paste(fname,'.nc',sep=''),ncvar )
#  print(newnc)
#    
#  # set global attributes
#  ncatt_put(newnc,0,attname='Conventions',attval='CF-1.4 (or close)')
#  ncatt_put(newnc,0,attname='Calendar',attval='no leap years, 360 day years')
#  ncatt_put(newnc,0,attname='institution',attval='Oak Ridge National Laboratory')
#  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by: Anthony Walker (walkerap@ornl.gov)'))
#  
#  # write the file
#  nc_close(newnc)  
#}
#
## puts data into a pre-existing netcdf
## could be modified to add data bit by bit for daily data
#data_netcdf <- function( df,var,fname=var,lname,units,nsites=1548,nyears=110,lon=3.75,lat=2.5,mv=-99999,
#                         annual=F,monthly=F,daily=F,cs=3 ){
#  # expects df in the mp processed SDGVM output format
#  
#  # initialise - leap years not yet considered
#  if(annual){
#    st <- 24*360    ; end <- 24*360*nyears ; sa <- 1
#  } else if(monthly){
#    st <- 24*30     ; end <- 24*360*nyears ; sa <- 12
#  } else if(daily){
#    st <- 24        ; end <- 24*360*nyears ; sa <- 360   
#  }
#  
#  # open the file given the var
#  newnc <- nc_open( paste(fname,'.nc',sep=''),write=T )
#  print(newnc)
#  
#  ### put the data into the ncfile - this needs modification to allow multiple subsets of the data to be entered into the netcdf
#  # create global array
#  df$lats  <- (df[,1] + 90)/lat  + 1
#  df$lons  <- (df[,2] + 180)/lon + 1
#  da       <- array(1:(360/lon*180/lat*nyears*sa),dim=c(360/lon,180/lat,nyears*sa))
#  da[]     <- NA
#  # get subscripts
#  time_sub <- unlist(lapply(1:nyears,slice,l=sa,nsites=nsites))
#  smat     <- cbind(rep(df$lons,nyears),rep(df$lats,nyears),time_sub)
#  # put data into matrix
#  da[smat] <- as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
#  rm(df)
#  rm(smat)
#  rm(time_sub)
#  # put data into netcdf
#  ncvar_put(newnc,var,da)
#  print(levelplot(da[,,101],main=var))
#  rm(da)
#  
#  # write the file
#  nc_close(newnc)  
#}
