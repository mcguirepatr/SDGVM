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
source('read_SDGVM_output.R')
  
  
  
########################
# Function to take concatenated SDGVM output and process to TRENDY netcdf output format
write_sdgvm_netcdf <- function(wd, afiles=NULL, mfiles=NULL, dfiles=NULL, site_vars=NULL,
                               fref=NULL,
                               mc=F, procs=4, ... ) {
  
  
  
  # Initialisation
  ########################
  setwd(wd)
 
  # scan the directory name for the simlation id, passed thru functions to generate output filename  
  if(is.null(fref)) {
    if(grepl('T1',wd)) fref   <- 'T1'
    if(grepl('T2',wd)) fref   <- 'T2'
    if(grepl('T3',wd)) fref   <- 'T3'
    if(grepl('Stest',wd)) fref <- 'Stest'
    if(grepl('spin_short',wd)) fref   <- 'spin_short'
    if(grepl('spin_accel',wd)) fref   <- 'spin_accel'
    if(grepl('S0',wd)) fref   <- 'S0'
    if(grepl('S1',wd)) fref   <- 'S1'
    if(grepl('S2',wd)) fref   <- 'S2'
    if(grepl('S3',wd)) fref   <- 'S3'
    if(grepl('S4',wd)) fref   <- 'S4'
    if(grepl('S5',wd)) fref   <- 'S5'
    if(grepl('S6',wd)) fref   <- 'S6'
    if(grepl('S7',wd)) fref   <- 'S7'
    if(grepl('S8',wd)) fref   <- 'S8'
    if(grepl('S9',wd)) fref   <- 'S9'
    if(grepl('tpu0',wd)) fref <- 'tpu0'
    if(grepl('tpu1',wd)) fref <- 'tpu1'
    if(grepl('tpu2',wd)) fref <- 'tpu2'
    if(grepl('tpu3',wd)) fref <- 'tpu3'
  }

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
  if(!is.null(site_vars)) make_netcdf_TRENDY_site(site_vars, daily=T, fref=fref, ... )        
}


slice <- function(i,l,nsites){
  s <- (i-1)*l + 1
  e <- i*l
  rep(s:e,nsites)
}


make_netcdf_TRENDY <- function(varo, 
                               annual=F, monthly=F, daily=F, fref='',
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


# test year vector for leap years
leap_year <- function(year) {
  ifelse((year%%4==0 & year%%100!=0) | year%%400==0, TRUE, FALSE )
}


# creates netcdf and call the read write function
make_netcdf_TRENDY_site <- function(vars='gpp', alma=T, pftfile=F, annual=F, monthly=F, daily=F, fref='',
                                    nsites=1, nyears=110, osyr=1901, lon=0.0, lat=0.0, mv=-99999,
                                    person='Anthony P. Walker', email='walkerap@ornl.gov', 
                                    institution='Oak Ridge National Laboratory',  ... ) {

  # read site info
  df_si <- read_SDGVM_site_info()
  lon   <- df_si$lon 
  lat   <- df_si$lat 

  # create the nc dimensions
  nclon <- ncdim_def(name='longitude', units='degrees_east', vals=lon )
  nclat <- ncdim_def(name='latitude', units='degrees_north', vals=lat )
  pft   <- if(pftfile) pftnames else NULL 
  if(!is.null(pft)) {
    ncpft <- ncdim_def(name='PFT', units='pft id, see global attributes', vals=1:length(pft) )
    print('',quote=F)
    print('',quote=F)
    print('Processing PFTs:',quote=F)
  } 
  
  # initialise time dimension 
  # - leap years not yet considered, need to be
  # annmonth is for variables that are annual output from SDGVM but are required in monthly format 
  lyears <- sum(leap_year(osyr:(osyr+nyears-1)))
  end <- 24*(365*nyears+lyears) 
  #if(annual&is.null(var$annmonth)) {
  #  st <- 24*365    ; sa <- 1
  #} else if(monthly|!is.null(var$annmonth)) {
  #  st <- 24*30     ; sa <- 12
  #} else if(daily) {
    st <- 24        ; sa <- 365   
  #}
  # provides time at mid-point of timestep 
  print('')
  print('')
  print(st)
  print(end)
  print(nyears)
  time_seq <- seq(st,end,st) - st/2
  nctime   <- ncdim_def( name='time',units=paste('hours since ',osyr,'-01-01 00:00:00',sep=''),vals=time_seq,unlim=T )
  
  
  # create the ncdf4 object(s) of the var(s)
  print('',quote=F)
  print('Processing dynamic variable:',quote=F)
  for( varo in vars ) {
    var <- get(varo)
    vname <- if(alma&!is.null(var$alma_name)) var$alma_name else var$name 
    print(vname,quote=F)

    if(!is.null(pft)) ncvar <- ncvar_def( name=vname, units=var$units, dim=list(nclon,nclat,ncpft,nctime), missval=mv, longname=var$lname )
    else              ncvar <- ncvar_def( name=vname, units=var$units, dim=list(nclon,nclat,nctime), missval=mv, longname=var$lname )

    # create ncvars vector
    ncvars <- if(varo==vars[1]) list(ncvar) else c(ncvars,list(ncvar)) 
  }


  # create additional site-level variables in netcdf
  print('',quote=F)
  print('Processing static variable:',quote=F)
  si_vars <- c('sand', 'silt', 'BD', 'carbon', 'wilt', 'field', 'sat' )
  for( varo in si_vars ) {
    var <- SDGVM_site_info_vars[[varo]]
    vname <- var$name 
    print(vname,quote=F)
    ncvar <- ncvar_def( name=vname, units=var$units, dim=list(nclon,nclat), missval=mv, longname=var$lname )

    # create ncvars vector
    ncvars <- c(ncvars,list(ncvar)) 
  }


  # create the nc file given the vars
  newnc <- nc_create(paste0(paste('SDGVM',fref,osyr,osyr+nyears-1)'.nc'), ncvars )
  
  
  # add static site variables to netcdf
  print('',quote=F)
  print(paste('Adding static variable:'),quote=F)
  siv <- 1
  for( varo in si_vars ) {
    var <- SDGVM_site_info_vars[[varo]]
    print(var$col_name,quote=F)
    ncvar_put(newnc, ncvars[[length(vars)+siv]], df_si[[var$col_name]] )
    siv <- siv + 1
  }

  
  # set attributes
  ncatt_put(newnc,'time',attname='calendar',attval='Gregorian, inc. leap years')

  # set global attributes
  conventions <- if(alma) 'ALMA' else 'CF-1.4 (or close)'
  ncatt_put(newnc,0,attname='title',attval=paste('SDGVM output for', project ))
  ncatt_put(newnc,0,attname='Conventions',attval=conventions)
  ncatt_put(newnc,0,attname='institution',attval=institution)
  ncatt_put(newnc,0,attname='history',attval=paste('created:',as.character(as.POSIXlt(Sys.time())),', by:',person,paste0('(',email,')')))
  if(!is.null(var$notes)) ncatt_put(newnc,0,attname='notes',attval=var$notes)
  if(!is.null(pft))       ncatt_put(newnc,0,attname='PFTs',attval=paste(paste(pft,collapse=' '),'. PFT distributions were derived by combining the LUH2v2h land-use and land-cover change database with the ESA GLCP 2014 data categorised according to SDGVM PFTs.',sep='') )
   
  
  # extract total veg biomass
  if('cVegd'%in%vars) {
    ss_cveg <- which(vars=='cVegd')
    vars    <- vars[-ss_cveg]
    ncvars  <- ncvars[-ss_cveg]
  }  

  # variable data reading & writing loop
  for( v in 1:length(vars) ) {
    var <- get(vars[v])
   
    print('',quote=F)
    print('',quote=F)
    print(paste('Reading variable:',var$name),quote=F)
   
    # data filename 
    fname <- var$file
    if(monthly) fname <- paste('monthly',var$file,sep='_')
    if(daily)   fname <- paste('daily',var$file,sep='_')
    cs <- 3
    if(monthly) cs <- 4

    # if by pft output required 
    if(is.null(pft)) {
      fnamefull <- paste(fname,'dat',sep='.')  
      newnc     <- readSDGVM_writeNCDF(fnamefull, var, newnc, ncvars[[v]], cs, NULL, 
                                       180, 360, nyears, 1, sa, mv, osyr, daily=daily, ... ) 
    } else { 
      for( p in pft ) {
        fnamefull <- paste(fname, '_', p, '.dat', sep='')
        newnc     <- readSDGVM_writeNCDF(fnamefull, var, newnc, ncvars[[v]], cs, which(pft==p),  
                                         180, 360, nyears, 1, sa, mv, osyr, daily=daily, ... ) 
      }
    }
  }


  # calculate total veg biomass from leaf, wood, and root mass
  if(exists('ss_cveg')) {
    if('cLeafd'%in%vars & 'cWoodd'%in%vars & 'cRootd'%in%vars ) {
      
      total_bio <- ncvar_get(newnc, 'cLeaf', start=c(1,1,1) )
      bio1      <- ncvar_get(newnc, 'cWood', start=c(1,1,1) )
      bio2      <- ncvar_get(newnc, 'cRoot', start=c(1,1,1) )
      total_bio <- total_bio + bio1 + bio2 
      ncvar_put(newnc, 'cVeg', total_bio )

    } else {
      warning('cVegd requested but not calculated as one or all of cLeafd, cWoodd, cRootd not requested.')
    }
  }


  # write the file
  print(newnc)
  nc_close(newnc)  
}



readSDGVM_writeNCDF <- function(fname, var, newnc, ncvar, cs, pftid=NULL, 
                                lat, lon, nyears, nsites, sa, mv, osyr, styr, daily=F, 
                                ... ) {

  # if SDGVM output is the same variable as required
  if(length(var$file)==1) {
    df <- if(nsites>1) read.table(fname) else read_SDGVM_daily(fname, styr, styr+nyears-1 )     
  
  # if SDGVM output needs to be combined to create a composite variable
  } else {  

    for(file in fname) {
      #df1 <- read.table(file)    
      df1 <- if(nsites>1) read.table(file) else read_SDGVM_daily(file, styr, styr+nyears-1 ) 

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
  if(fname[1]=='monthly_nep.dat') {

    #print('Made it to routine to add C fluxes to nep')

    m1 <- as.matrix(read.table('lch.dat'))
    m2 <- as.matrix(read.table('lulccc.dat'))
    m3 <- as.matrix(read.table('fcn.dat'))

    df[,cs:length(df)] <- df[,cs:length(df)] - 
				apply(m1[,(dim(m1)[2]-length(df)+cs):dim(m1)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m2[,(dim(m2)[2]-length(df)+cs):dim(m2)[2]], 2, function(v) rep(v,each=12)/12 ) - 
				apply(m3[,(dim(m3)[2]-length(df)+cs):dim(m3)[2]], 2, function(v) rep(v,each=12)/12 ) 
    rm(m1,m2,m3)
  } 

  # for biomass variables add previous year biomass to current year 'living' (i.e. respiring) biomass
  if(var$name=='cWood' | var$name=='cRoot') {
    t 
    ifile <- if(var$name=='cWood') 'stembio.dat' else 'rootbio.dat'
    v_bio  <- scan(ifile) 
    v_pbio <- scan(paste0('../ind/',ifile)) 
    v_cbio <- c(v_pbio[length(v_pbio)], v_bio[3:(length(v_bio)-1)] )

    days_in_y <- lyears <- leap_year(osyr:(osyr+nyears-1))
    days_in_y[] <- 365
    days_in_y[lyears] <- 366

    repv    <- function(v) rep(v[1],each=v[2])
    v_dbio  <- unlist(apply(cbind(v_cbio,days_in_y),1,repv))
    df[,cs] <- df[,cs] + v_dbio
  } 


  # scale variable to correct output units
  scale_day <- if(daily&!is.null(var$scale_day)) var$scale_day else 1 
  if( var$name %in% c('tas','ta','tsl') )  df[,cs:length(df)] <- df[,cs:length(df)] + var$scale
  else                                     df[,cs:length(df)] <- df[,cs:length(df)] * var$scale * scale_day

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
  if(lon==360) { 
    tsteps <- length(df[,1])  # APW: needs work, length(df) incorrect
    smat   <- if(!is.null(pftid)) cbind(1, 1, 1, 1:tsteps ) else 
                                  cbind(1, 1, 1:tsteps )
    ce     <- 3
  } else {
    tsteps   <- nyears*sa
    time_sub <- unlist(lapply(1:nyears, slice, l=sa, nsites=nsites ))
#    print(nyears)
#    print(sa)
#    print(nsites)
#    print(time_sub)
#    print(df)
    smat     <- if(!is.null(pftid)) cbind(rep(df$lons,nyears), rep(df$lats,nyears), 1, time_sub ) else 
                                    cbind(rep(df$lons,nyears), rep(df$lats,nyears), time_sub )
    rm(time_sub)
    ce <- nyears + cs - 1
  }

  # create an nc ready array
  a_dim    <- if(!is.null(pftid)) c(360/lon,180/lat,1,tsteps) else 
                                  c(360/lon,180/lat,tsteps)
  da       <- array(mv,dim=a_dim)

  da[smat] <- as.vector(as.matrix(df[,cs:ce])) 
  rm(df); rm(smat)
  
  # put the data into the ncfile
  if(!is.null(pftid)) ncvar_put(newnc, ncvar, da, start=c(1,1,pftid,1), count=a_dim )
  else                ncvar_put(newnc, ncvar, da )
  rm(da)

  newnc
}



### END ###
