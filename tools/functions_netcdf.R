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
  
  
########################
# Function to take concatenated SDGVM output and process to TRENDY netcdf output format
write_sdgvm_netcdf <- function(wd, afiles=NULL, mfiles=NULL, dfiles=NULL,
                               mc=F, procs=4, ... ) {
  
  
  
  # Initialisation
  ########################
  setwd(wd)
 
  # scan the directory name for the simlation id, passed thru functions to generate output filename  
  if(grepl('T1',wd)) fref   <- 'T1'
  if(grepl('T2',wd)) fref   <- 'T2'
  if(grepl('T3',wd)) fref   <- 'T3'
  if(grepl('Stest',wd)) fref <- 'Stest'
  if(grepl('spin_short',wd)) fref   <- 'spin_short'
  if(grepl('spin_accelFS4',wd)) fref   <- 'spin_accelFS4'
  else { if(grepl('spin_accel',wd)) fref   <- 'spin_accel'}
  if(grepl('S0',wd)) fref   <- 'S0'
  if(grepl('S1',wd)) fref   <- 'S1'
  if(grepl('S2',wd)) fref   <- 'S2'
  if(grepl('S3',wd)) fref   <- 'S3'
  if(grepl('FS4',wd)) fref   <- 'FS4' 
  else {if(grepl('S4',wd)) fref   <- 'S4'}
  if(grepl('FS5',wd)) fref   <- 'FS5'
  else {if(grepl('S5',wd)) fref   <- 'S5'}
  if(grepl('S6',wd)) fref   <- 'S6'
  if(grepl('S7',wd)) fref   <- 'S7'
  if(grepl('S8',wd)) fref   <- 'S8'
  if(grepl('S9',wd)) fref   <- 'S9'
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
           
  
  # data processing - call make_netcdf_TRENDY function
  ########################

  # parallel process or not
  if(mc){
    if(annual)  mclapply(afiles,make_netcdf_TRENDY,mc.cores=procs,annual=annual,fref=fref,...)
  } else {
    if(annual)  lapply(afiles,make_netcdf_TRENDY,annual=annual,fref=fref,...)
  }
  if(mc){
    if(monthly) mclapply(mfiles,make_netcdf_TRENDY,mc.cores=procs,monthly=monthly,fref=fref,...)
  } else {
    if(monthly) lapply(mfiles,make_netcdf_TRENDY,monthly=monthly,fref=fref,...)
  }
  if(daily)   lapply(dfiles,make_netcdf_TRENDY,daily=daily,fref=fref,...)        

}


slice <- function(i,l,nsites){
  s <- (i-1)*l + 1
  e <- i*l
  rep(s:e,nsites)
}


make_netcdf_TRENDY <- function(varo, annual=F, monthly=F, daily=F, fref='',
                               nsites=1548, nyears=110, osyr=1901, lon=3.75, lat=2.5, mv=-99999,
                               person='Anthony P. Walker', email='walkerap@ornl.gov', institution='Oak Ridge National Laboratory',  ... ) {

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
  print('')
  print('')
  print(st)
  print(end)
  print(nyears)
  time_seq <- seq(st,end,st) - st/2
  

  # create the nc dimensions
  nclon   <- ncdim_def( name='longitude',units='degrees_east',vals=(seq(-180+lon/2,180-lon/2,lon)) )
  nclat   <- ncdim_def( name='latitude',units='degrees_north',vals=(seq(-90+lat/2,90-lat/2,lat)) )
  nctime  <- ncdim_def( name='time',units=paste('hours since ',osyr,'-01-01 00:00:00',sep=''),vals=time_seq,unlim=T )
  if(!is.null(pft)) ncpft <- ncdim_def( name='PFT',units='pft id, see global attributes',vals=1:length(pft) )
  

  # create the ncdf4 object of the var
  if(!is.null(pft)) ncvar <- ncvar_def( name=var$name,units=var$units,dim=list(nclon,nclat,ncpft,nctime),missval=mv,longname=var$lname )
  else              ncvar <- ncvar_def( name=var$name,units=var$units,dim=list(nclon,nclat,nctime),missval=mv,longname=var$lname )
  

  # create the nc file given the var
  newnc <- nc_create( paste('SDGVM_',fref,'_',var$name,'.nc',sep='') , ncvar )
  

  # set attributes
  ncatt_put(newnc,'time',attname='calendar',attval='360_day')

  # determine global attributes

  # set global attributes
  ncatt_put(newnc,0,attname='title',attval=paste('SDGVM output for', project ))
  ncatt_put(newnc,0,attname='Conventions',attval='CF-1.4 (or close)')
  ncatt_put(newnc,0,attname='calendar',attval='no leap years, 360 day years')
  ncatt_put(newnc,0,attname='institution',attval=institution)
  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by:',person,paste0('(',email,')')))
  if(!is.null(var$notes)) ncatt_put(newnc,0,attname='notes',attval=var$notes)
  if(!is.null(pft))       ncatt_put(newnc,0,attname='PFTs',attval=paste(paste(pft,collapse=' '),'. PFT distributions were derived by combining the LUH2v2h land-use and land-cover change database with the ESA GLCP 2014 data categorised according to SDGVM PFTs.',sep='') )
   

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
                                     lat,lon,nyears,nsites,sa,mv,osyr,...) 
  } else { 
    for( p in pft ) {
      fnamefull <- paste(fname,'_',p,'.dat',sep='')
      newnc     <- readSDGVM_writeNCDF(fnamefull,var,newnc,ncvar,cs,which(pft==p), 
                                       lat,lon,nyears,nsites,sa,mv,osyr,...) 
    }
  }


  # write the file
  print(newnc)
  nc_close(newnc)  
}


readSDGVM_writeNCDF <- function(fname,var,newnc,ncvar,cs,pftid=NULL,
                                lat,lon,nyears,nsites,sa,mv,osyr,styr,
                                ... ) {

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



  # for monthly grid square NBP remove fire flux, lulcc flux, yield flux, and leached DOC flux
  #PCM if(fname=='monthly_nep.dat') {
  if(var$name=='nbp') { #PCM

    print('Made it to routine to add C fluxes to nep')

    m1 <- as.matrix(read.table('lch.dat'))
    m2 <- as.matrix(read.table('lulccc.dat'))
    m3 <- as.matrix(read.table('fcn.dat'))
    m4 <- as.matrix(read.table('yield.dat')) #PCM

    #print('')
    #print(df[1:12,cs:length(df)])  
    #print('')
    #print(apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] )  
    #print('')
    #print(apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] )  
    #print('')
    #print(apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] ) 
    #print('')
    #print(apply(m4[,(dim(m4)[2]-length(df)+cs):dim(m4)[2]], 2, function(v) rep(v,each=12)/12 )[1:12,] ) 

    print('')
    print(sum(df[,cs:length(df)]))  
    print(sum(apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 )))  
    print(sum(apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 )))  
    print(sum(apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 ))) 
    print(sum(apply(m4[,(dim(m4)[2]-length(df)+cs):dim(m4)[2]], 2, function(v) rep(v,each=12)/12 ))) 

    df[,cs:length(df)] <- df[,cs:length(df)] - 
				apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m4[,(dim(m4)[2]-length(df)+cs):dim(m4)[2]], 2, function(v) rep(v,each=12)/12 ) 
    rm(m1,m2,m3,m4)
  } 


  # scale variable to correct output units
  if( var$name %in% c('tas') )  df[,cs:length(df)] <- df[,cs:length(df)] + var$scale
  else                          df[,cs:length(df)] <- df[,cs:length(df)] * var$scale

  # expand annual SDGVM variables to monthly TRENDY output 
  if(!is.null(var$annmonth)) {
    df1        <- as.data.frame(matrix(0,dim(df)[1]*12,dim(df)[2]+1))
    names(df1) <- c(names(df)[1:2],'month',names(df)[cs:dim(df)[2]])
    df1[,1]    <- rep(df[,1], each=12 )
    df1[,2]    <- rep(df[,2], each=12 )
    df1[,3]    <- 1:12
    df1[,4:dim(df1)[2]] <- apply(as.matrix(df[,cs:dim(df)[2]]), 2, function(v) rep(v,each=12)/12 )
    df <- df1
    cs <- 4
    rm(df1)
  } 

  # rename lats and lons
  df$lats  <- (df[,1] + 90)/lat  + 1
  df$lons  <- (df[,2] + 180)/lon + 1

  # put the data into an nc ready array
  time_sub <- unlist(lapply(1:nyears, slice, l=sa, nsites=nsites ))
#  print(nyears)
#  print(sa)
#  print(nsites)
#  print(time_sub)
#  print(df)
  smat     <- if(!is.null(pftid)) cbind(rep(df$lons,nyears), rep(df$lats,nyears), 1, time_sub ) else 
                                  cbind(rep(df$lons,nyears), rep(df$lats,nyears), time_sub )
  rm(time_sub)
  
  # create an nc ready array 
  a_dim    <- if(!is.null(pftid)) c(360/lon,180/lat,1,nyears*sa) else 
                                  c(360/lon,180/lat,nyears*sa)
  da       <- array(mv,dim=a_dim)

#  PCM: commented out the following block
#  # remove unneeded columns from df & convert to matrix
#  # lat & lon
#  df       <- df[,-c(1:(cs-1))]
#  cdim     <- dim(df)[2] - 2
#  df       <- df[,-c((cdim+1):(cdim+2))]
#  # years of output not requested
#  #if(nyears!=cdim) df <- df[,-(1:(cdim-nyears))]
#  if(nyears!=cdim) {
#    osyr_ss <- osyr - styr
#    df <- df[,-(1:osyr_ss)]
#    if(nyears!=dim(df)[2]) df <- df[,1:nyears]
#  }
#  df <- as.matrix(df)

  # put the data into an nc ready array 
  #as.vector(as.matrix(df[,cs:(nyears+cs-1)]))
  #da[smat] 
  #print('Read file2(s):',quote=F)
  #print(fname,quote=F)
  #print(head(df),quote=F)
  #print('',quote=F)
  #print('Read file3(s):',quote=F)
  #print(fname,quote=F)
  #print(head(df[,cs:(nyears+cs-1)]),quote=F)
  #print('',quote=F)
  da[smat] <- as.vector(as.matrix(df[,cs:(nyears+cs-1)])) #PCM
  #print('Read file4(s):',quote=F)
  #print(fname,quote=F)
  #print(da[smat],quote=F)
  #print('',quote=F)
#PCM  #da[smat] <- as.vector(df)
  rm(df); rm(smat)
  
  # put the data into the ncfile
  if(!is.null(pftid)) ncvar_put(newnc, ncvar, da, start=c(1,1,pftid,1), count=a_dim )
  else                ncvar_put(newnc, ncvar, da )
  rm(da)

  newnc
}



### END ###
