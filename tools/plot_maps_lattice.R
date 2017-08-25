####################################
#
# global plotting function
#
# AWalker
# Jun 2015
#
###################################

rm(list=ls())

.libPaths('~/bin/Rlibs')
library(lattice)
library(plyr)
#library(dplyr)
library(viridis)


##################################
###user defined inputs

#directory paths
date    <- '170822'
dir     <- '/home/alp/models/SDGVM/'
rdir    <- 'run'
edir    <- 'eval_data'

# project directory
project <- 'vcmax'

# simulations
sim     <- c('Constant','Kattge','Kattge_oxisol','LUNA','Maire','vBodegom_env','vBodegom_mean','Walker_N','Walker_NP','Woodward_95')

# simulations index array (which simulations to plot from the vector 'sim', when 'sim' is arranged in alphabetical order, i.e. order(sim)[sia], simulations appear on the plot in the sia order)
sia     <- NULL

# only plot eval data, use first value of sim to stadardise by model mask 
evalonly <- F

# alternative labels for sims, should be length(sim) or if sia is not NULL length(sia) long
# - applied to the sims following order(sim)[sia]  
slabels <- NULL

# output file naming prefix
of      <- NULL 

# main name for plots - should be as long as sia (maybe sia * via) but can't be bothered to do this right now
main    <- NULL
 
# decimal places for scale legend
r       <- 1 

#  allow manually set columns (cs) and rows (rs)
cs      <- NULL 
rs      <- NULL
 
# logical vector decribing which panels in the mapplot to skip  
skip    <- NULL

# indexing list (index.cond variable in xyplot) each element of the list is an indexing vector corresponding to the order in which to plot each conditioning factor  
icon    <- NULL

# map data resolution  
deg1    <- T

# data start year
styr    <- 1901

# data end year
endyr   <- 2012

# plotting region
pregion <- 'global'

# plot year (can be a vector of years)
pyear   <- NULL 

# years over which to take a mean value, must be a two element vector 
yr_mean <- c(2001,2010)

# plot trends across whole timeseries
trends  <- T 

# write trends across whole timeseries
trend_out  <- F

# zero trend on initial year value
trend_norm <- F

# this should be a character vector. This will break up the zonal lat and lon plots into length(zon_panel_names) panels
# - the sims will be evenly distributed among the panels
zon_panel_names <- NULL

# output a table of global integrated values 
table     <- F 

# plot scatterplot/s of var v against var v-1
scatter   <- F 

# normalise the data to the norm*100 %ile, or 'sd' or 'center'
norm      <- NULL 

# output a table of global integrated NBP values for each year, for TRENDY project 
trendy    <- F 

# make a difference plot, takes the value of the simulation/dataset id string and uses it as the base for the differences
diff      <- NULL

# run EOF analysis
EOF       <- F

# load SIF data
SIF       <- F

# adjust SIF by mean of SIFadj (can be a GPP proxy or a model simulation, e.g. SIFadj <- 'MPI' )
SIFadj    <- NULL

# load SIFGPP data
SIFGPP    <- F

# load MPI GPP data
MPI       <- F

# load PR GPP data
PR        <- F

# print map with masked area (current options are: crop, bare; or any combination of these terms in a single string ) 
mask      <- NULL 
mask_perc <- 40

# line colour & type on trend and zonal plots, correspond to the labels in 'sim' when 'sim' is in alphabetical order
col_trend <- NULL 
lty       <- 1 

# plot file format
plotype   <- 'pdf'

# figure dimensions etc
width     <- 180/25.4 
width1    <-  80/25.4 
height    <- 210/25.4
pointsize <- 10
labcex    <- 1

# backgroud for plots
backg     <- 'transparent'
backg     <- 'white'

# plotting variables (these all have an associated list in 'params_map_plot.R', to add variables simply add a list in 'params_map_plot.R' with the same name as the variable).
vars <- c('npp','gpp','nbp','anlfn','antlfn',        
          'evt','trn','scn','sresp','presp',
          'mgresp','lai','anvcmax','anjmax','biot',
          'kg_beta','swr','qtotal','tmp','prc',
          'swc','field_capacity','wilting_point','cov_C3','cov_C4',
          'cov_C3crop','cov_C4crop','cov_Dc_Bl','cov_Dc_Nl','cov_Ev_Bl',
          'cov_Ev_Nl','cov_BARE','fcn','lulccc','fab')

# land-cover is fixed and so has only a single column in the output file
cov_fixed  <- T

# variable index array (which variables to plot, in the order they appear in the 'vars' vector)
via <- 1:length(vars)
#via <- 24:length(vars)
via <- 1:16



### Parse command line arguments   
##########################
# any one of the above objects can be specified as a command line argument using the syntax:
# Rscript <nameofthisscript> "<object1><-<value1>" "<object2><-<value2>"
# e.g. Rscript plot_maps_lattice.R "sty<-2000" "dir<-'/home/alp/models/SDGVM'" "vars<-c('npp','gpp')"

print('',quote=F)
print('Read command line arguments:',quote=F)
if(length(commandArgs(T))>=1) {
  for( ca in 1:length(commandArgs(T)) ) {
    eval(parse(text=commandArgs(T)[ca]))
  }
}
print(commandArgs(T))

# set default values if not specified on command line
# these are set after parsing command line arguments as they depend on other arguments that could be set on the command line
if(is.null(sia)&!evalonly) sia <- 1:length(sim)
#if(is.null(sia))           sia       <- 1:length(sim)
if(evalonly)               sia <- 1
if(is.null(of))            of  <- project
#if(is.null(lty))           lty <- rep(1:length(sia))
if(is.null(lty))           lty       <- 1
if(is.null(col_trend)) col_trend <- viridis(length(sia))
#year <- if(is.null(pyear)) 2010 else pyear
year <- if(is.null(pyear)) yr_mean[2] else pyear

# set model resolution
         mod_res <- c(3.75,2.5)
if(deg1) mod_res <- c(1,1)

# create output directory
owd <- paste(dir,rdir,project,'results',date,sep='/')
if(!file.exists(owd)) dir.create(owd)

# evaluation data directory
ewd <- paste(dir,edir,sep='')

# order simulations alphabetically
sim <- sim[order(sim)]
print('',quote=F)
print('Simulations requested:',quote=F)
print(sim[sia],quote=F)
if(evalonly) print('simulation only used to standardise eval data',quote=F)
print('',quote=F)
print('Variables requested:',quote=F)
print(vars[via],quote=F)
print('',quote=F)

# normalisation function
fnorm <- function(df,norm) {
  # takes a dataframe and normalises the 'plotdata' column
  # norm can be a character string 'sd' or 'center' or a number between 0-1 which will be the quantile to normalise on

  center <- mean(df$plotdata,na.rm=T) 
  sd     <- var(df$plotdata,na.rm=T)^0.5 
  if(norm=='sd'){ 
    df$plotdata <- (df$plotdata - center) / sd 
  } else if(norm=='center') {
    df$plotdata <- df$plotdata - center 
  } else {

    if(any(df$plotdata<0)) {
      df$plotdata <- df$plotdata / max( abs(quantile(df$plotdata,1-norm,type=8)) , abs(quantile(df$plotdata,norm,type=8)) )
    } else {
      df$plotdata <- df$plotdata / quantile(df$plotdata,norm,type=8)
    }
  }
  df
} 


###############################
### start program

setwd(paste(dir,'src/sdgvm/tools/',sep='/'))
source('params_map_plot.R')
source('functions_map_plot.R')

# load SIF data
if(SIF) {
  print('',quote=F)
  print('read SIF',quote=F)
  setwd(ewd)
  SIF_df <- read.table('SIF_l3_2007t2012mean_screened_1x1.dat',header=F)
  names(SIF_df) <- c('lat','lon','plotdata')
  SIF_df$run    <- 'SIF' 
  SIF_df$year   <- 2009.5  
}

# load SIF GPP data
if(SIFGPP) {
  print('',quote=F)
  print('read SIFGPP',quote=F)
  setwd(ewd)
  SIFGPP_df <- read.table('CASA-SIFmax_GPP_2007t2012mean_1x1.dat',header=F)
  names(SIFGPP_df) <- c('lat','lon','plotdata')
  SIFGPP_df$run    <- 'SIFGPP' 
  SIFGPP_df$year   <- 2009.5  
}

# load MPI-GPP data
if(MPI) {
  print('',quote=F)
  print('read MPI-GPP',quote=F)
  setwd(ewd)
  MPI_df <- read.table('MPI-GPP_2007t2010mean_screened_1x1.dat',header=F)
  names(MPI_df) <- c('lat','lon','plotdata')
  MPI_df$run    <- 'MPI' 
  MPI_df$year   <- 2009  
}

# load PR-GPP data
if(PR) {
  print('',quote=F)
  print('read PR-GPP',quote=F)
  setwd(ewd)
  PR_df <- read.table('PR-GPP_2007t2012mean_screened_1x1.dat',header=F)
  names(PR_df)  <- c('lat','lon','plotdata')
  PR_df$run     <- 'PR' 
  PR_df$year    <- 2009.5  
}

# load mask
kill <- F
mask_exists <- F
if(!is.null(mask)) {
  
  print('',quote=F)
  print(paste('read mask:',mask,mask_perc,'perc'),quote=F)

  wdpath  <- paste(dir,rdir,project,sim[sia[1]],'output/',sep='/')

  if(grepl('crop',mask)) {
    mydata1     <- open('cov_C3crop',c('lat','lon','cov'),wdpath)
    mydata2     <- open('cov_C4crop',c('lat','lon','cov'),wdpath)
    if(!mask_exists) mydata <- mydata1  
    mydata$cov  <- mydata$cov + mydata2$cov 
    mask_exists <- T 
    rm(mydata1)  
    rm(mydata2)  
  } 
  if(grepl('bare',mask)) {
    mydata1     <- open('cov_BARE',c('lat','lon','cov'),wdpath)
    mydata2     <- open('cov_CITY',c('lat','lon','cov'),wdpath)
    if(!mask_exists) mydata <- mydata1  
    mydata$cov  <- mydata$cov + mydata1$cov + mydata2$cov 
    mask_exists <- T 
    rm(mydata1)  
    rm(mydata2)  
  } #else {
    #print(paste('no mask method for mask type:',mask))
    #print('results will not be masked')
    #kill <- T
  #}

  # establish mask
  mydata$cov[mydata$cov > mask_perc/100] <- NA

  # subset to region
  lims <- region(pregion)
  mds  <- subset(mydata, lon>(lims[1]-mod_res[1]/2))
  mds  <- subset(mds,    lon<(lims[2]+mod_res[1]/2))
  mds  <- subset(mds,    lat>(lims[3]-mod_res[2]/2))
  mds  <- subset(mds,    lat<(lims[4]+mod_res[2]/2))      
  mask_v <- mds$cov
} 
 

# layout parameters
nruns <- length(sia) + SIF + SIFGPP + MPI - !is.null(diff) - evalonly
if(is.null(rs)){
 rs <- ceiling((nruns*length(year))^0.5)
 if((nruns*length(year)) < 4) rs <- nruns*length(year)
 if(pregion=='tropics')             rs <- rs + 2
 if(!is.null(cs)) rs <- ceiling(nruns*length(year)/cs)
}
if(is.null(cs)) cs <- ceiling(nruns*length(year)/rs)
print('',quote=F)
print(paste('rows:',rs,'cols:',cs),quote=F)
 
 
# variable plotting loop
vars <- vars[via]
# change this to an mclapply -  so far there are a few too many arguments
#mclapply(1:length(vars),,)

for( v in 1:length(vars) ) {
  print('',quote=F)
  print(vars[v],quote=F)
  
  lab <- get(vars[v])

  stripprint <- T 
  if(vars[v]=='prc'|vars[v]=='tmp'|vars[v]=='swr'|vars[v]=='qtotal') stripprint <- F

  # simulations loop
  for(m in sia) {    
    
    #open data 
    wdpath <- paste(dir,rdir,project,sim[m],'output/',sep='/')
    print(wdpath)
    if(substr(vars[v],1,3)=='cov'&cov_fixed) {
      index  <- 3
      mydata <- open(vars[v],c('lat','lon','cov'),wdpath)
      rs     <- 1 
      cs     <- 1
      stripprint <- F
    } else {
      index  <- year-styr+3
      mydata <- open(vars[v],c('lat','lon',paste('y',styr:endyr,sep='')),wdpath)
    }
    
    # subset to region
    lims <- region(pregion)
    mds  <- subset(mydata, lon>(lims[1]-mod_res[1]/2))
    mds  <- subset(mds,    lon<(lims[2]+mod_res[1]/2))
    mds  <- subset(mds,    lat>(lims[3]-mod_res[2]/2))
    mds  <- subset(mds,    lat<(lims[4]+mod_res[2]/2))      
    mydata <- mds
 
    if(!is.null(mask)) { 
      # check mask is correct length
      if( length(mydata[,1])%%length(mask_v) != 0 ) {
        print('mask is not the same length as data')
        stop 
      }

      # apply mask
      mydata <- subset(mydata,!is.na(mask_v))
    }

    # calculate mean across years
    if( !(substr(vars[v],1,3)=='cov'&cov_fixed) ) {
      ym             <- apply(as.matrix(mydata[,(yr_mean[1]-styr+3):(yr_mean[2]-styr+3)]),1,mean)
      ym             <- cbind(mydata[,1:2],ym,sim[m],mean(yr_mean))
      names(ym)[3:5] <- c('plotdata','run','year')
    }

    # calculate area integrated values
    if((lab$gsum|lab$gmean)&is.null(diff)) {

      # area integrate
      areai   <- area_integrate(mydata[,1:2],mydata[,3:length(mydata)],mod_res,lab$gmean)
      gs      <- apply(as.matrix(areai),2,sum,na.rm=T)   
      mean_gs <- mean(gs[(yr_mean[1]-styr+1):(yr_mean[2]-styr+1)])  
      global_sum      <- if(m==sia[1])      gs else rbind(global_sum,gs)
      mean_global_sum <- if(m==sia[1]) mean_gs else c(mean_global_sum,mean_gs)
      #mean_global_sum <- if(m==sia[1]) mean(gs[(yr_mean[1]-styr+1):(yr_mean[2]-styr+1)]) else 
      #                   c(mean_global_sum,mean(gs[(yr_mean[1]-styr+1):(yr_mean[2]-styr+1)]))      

      # arrange zonal data
      # - area integrated annual mean values using the year range specified in yr_mean 
      #if(substr(vars[v],1,3)!='cov') {
      if( !(substr(vars[v],1,3)=='cov'&cov_fixed) ) {
        print('Make zonal data')

        ymai  <- apply(as.matrix(areai[,(yr_mean[1]-styr+1):(yr_mean[2]-styr+1)]),1,mean)
        zmlat <- as.data.frame.table(tapply(ymai,mydata$lat,sum))
        zmlon <- as.data.frame.table(tapply(ymai,mydata$lon,sum))
      
        zmlat[,1] <- as.numeric(as.character(zmlat[,1]))
        zmlon[,1] <- as.numeric(as.character(zmlon[,1]))

        zonal_lat <- if(m==sia[1]) t(zmlat) else rbind(zonal_lat,zmlat[,2])
        zonal_lon <- if(m==sia[1]) t(zmlon) else rbind(zonal_lon,zmlon[,2])      
     
      } else { 
      
        zonal_lat  <- NULL
        zonal_lon  <- NULL
      }
      
      if(trendy&vars[v]=='nbp') { 

        trends <- F
        # zonal/regional sumns for TRENDY NBP
        trendy_df <- data.frame(year=styr:endyr,global=gs)
        # Northern extra-tropics
        sub_mydata <- subset(mydata,lat>30)
        areai      <- area_integrate(sub_mydata[,1:2],sub_mydata[,3:length(sub_mydata)],mod_res,lab$gmean)
        trendy_df$northern <- apply(as.matrix(areai),2,sum,na.rm=T)   
        # tropics
        sub_mydata <- subset(mydata,lat<=30&lat>= -30)
        areai      <- area_integrate(sub_mydata[,1:2],sub_mydata[,3:length(sub_mydata)],mod_res,lab$gmean)
        trendy_df$tropics <- apply(as.matrix(areai),2,sum,na.rm=T)   
        # Southern extra-tropics
        sub_mydata <- subset(mydata,lat< -30)
        areai      <- area_integrate(sub_mydata[,1:2],sub_mydata[,3:length(sub_mydata)],mod_res,lab$gmean)
        trendy_df$southern <- apply(as.matrix(areai),2,sum,na.rm=T)   
   
        print('Write TRENDY NBP csv')
        setwd(wdpath) 
        write.csv(trendy_df,paste('SDGVM_',sim[m],'_nbp.csv',sep=''),row.names=F,quote=F )
        setwd(owd) 
        write.csv(trendy_df,paste('SDGVM_',sim[m],'_nbp.csv',sep=''),row.names=F,quote=F )
      }  

    } else { 
      
      global_sum      <- NULL
      mean_global_sum <- NULL
      zonal_lat       <- NULL
      zonal_lon       <- NULL
        
    }
        
    # stack years for lattice plotting
    stack <- function(col,df){
      df1 <- df[,c(1:2,col)]
      df1$year <- col+styr-3
      names(df1)[3] <- 'plotdata'
      df1
    }
    #print(index)
    #print(head(mydata))
    dfl <- lapply(index,stack,df=mydata) 
    df1 <- rbind.fill(dfl)
    df1$run <- sim[m]    
  
    # normalise data
    if(!is.null(norm)) df1 <- fnorm(df1,norm) 
    if(!is.null(norm)) ym  <- fnorm(ym,norm) 
  
    # make dataframe with all simulations
    df      <- if(m==sia[1]) df1 else rbind(df,df1)
    #if(substr(vars[v],1,3)!='cov') mean_df <- if(m==sia[1])  ym else rbind(mean_df,ym)      
    if( !(substr(vars[v],1,3)=='cov'&cov_fixed) )  mean_df <- if(m==sia[1])  ym else rbind(mean_df,ym)

  # end simulations loop 
  }

  if(substr(vars[v],1,3)=='cov'&cov_fixed) mean_df <- df

  # prevent global sums and zonal data being over-written when only a single simulation called with eval datasets  
  if(!evalonly) m <- sia[1] + 1  

  # plotting order 
  df$run      <- factor(df$run,levels=sim[sia])
  mean_df$run <- factor(mean_df$run,levels=sim[sia])
  if(!is.null(slabels)) {
    levels(df$run)      <- slabels
    levels(mean_df$run) <- slabels
  }

  # append datasets with eval data
  if(SIF|SIFGPP|MPI|PR) {

    evars <- c('MPI','PR','SIFGPP','SIF')[which(c(MPI,PR,SIFGPP,SIF))]
    labs  <- c('MPI','WangFP','SIF-CASA','SIF')[which(c(MPI,PR,SIFGPP,SIF))]

    for( e in evars ) {
      # restrict dataset to pixels that occur in the model dataset only
      print(e,quote=F)
      # hijacking the now unused df1 dataframe from the simulations loop
      df1$plotdata <- NA
      df1$run      <- labs[which(evars==e)]
      edf          <- get(paste(e,'_df',sep=''))
      df1$plotdata <- edf$plotdata[match(paste(df1$lat,df1$lon),paste(edf$lat,edf$lon))]
      edf 	   <- df1
      if(!is.null(norm)) edf <- fnorm(edf,norm) 
      df      <- if(m==sia[1]) edf else rbind(df,edf)
      mean_df <- if(m==sia[1]) edf else rbind(mean_df,edf)

      # add global sums and zonal means 
      if((lab$gsum|lab$gmean)&is.null(diff)) {
  
        # area integrate
        # edf is already a mean annual value 
        areai   <- area_integrate(edf[,1:2],edf$plotdata,mod_res,lab$gmean)
        print(head(areai),quote=F)
        mean_gs <- sum(areai,na.rm=T)
        mean_global_sum <- if(m==sia[1]) mean_gs else c(mean_global_sum,mean_gs)
  
        # arrange zonal data
        # - area integrated annual mean values using the year range specified in yr_mean 
        ymai  <- areai 
        zmlat <- as.data.frame.table(tapply(ymai,edf$lat,sum))
        zmlon <- as.data.frame.table(tapply(ymai,edf$lon,sum))
      
        zmlat[,1] <- as.numeric(as.character(zmlat[,1]))
        zmlon[,1] <- as.numeric(as.character(zmlon[,1]))
  
        zonal_lat <- if(m==sia[1]) t(zmlat) else rbind(zonal_lat,zmlat[,2])
        zonal_lon <- if(m==sia[1]) t(zmlon) else rbind(zonal_lon,zmlon[,2])      
      }
      
      if(m==sia[1]) m <- m + 1 
      
    }

    # this is poor programing but will work for now to add correct labels and keep input consistent 
    if(!is.null(SIFadj)) if(SIFadj=='SIFGPP') SIFadj <- 'SIF-CASA'
    if(!is.null(diff))   if(diff=='SIFGPP')     diff <- 'SIF-CASA'

    # adjust SIF by other dataset
    if(!is.null(SIFadj)) {
      print('',quote=F)
      print(paste('Adjust SIF by:',SIFadj,'mean'),quote=F)

      if(class(df$run)=='character') {
        df$run      <- as.factor(df$run)
        mean_df$run <- as.factor(mean_df$run)
        print(levels(df$run))
      }

      sifsubs  <- which(df$run=='SIF')
      normsubs <- which(df$run==SIFadj)
      df$plotdata[sifsubs] <- ( df$plotdata[sifsubs] / mean(df$plotdata[sifsubs],na.rm=T) ) * mean(df$plotdata[normsubs],na.rm=T) 
      levels(df$run)[which(levels(df$run)=='SIF')] <- 'scaled-SIF'

      sifsubs  <- which(mean_df$run=='SIF')
      normsubs <- which(mean_df$run==SIFadj)
      mean_df$plotdata[sifsubs] <- ( mean_df$plotdata[sifsubs] / mean(mean_df$plotdata[sifsubs],na.rm=T) ) * mean(mean_df$plotdata[normsubs],na.rm=T) 
      levels(mean_df$run)[which(levels(mean_df$run)=='SIF')] <- 'scaled-SIF'
      print(length(which(is.na(mean_df$plotdata[sifsubs]))))
      if(!is.null(diff)) if(diff=='SIF') diff <- 'scaled-SIF'

      # add global sums and zonal means 
      if((lab$gsum|lab$gmean)&is.null(diff)) {
  
        # area integrate
        print(head(edf),quote=F)
        # edf is already a mean annual value 
        areai   <- area_integrate(mean_df[sifsubs,1:2],mean_df$plotdata[sifsubs],mod_res,lab$gmean)
        mean_gs <- sum(areai)
        sif_sub <- length(sia) + sum(c(MPI,PR,SIFGPP)) + !evalonly
      #  sif_sub <- length(sia) + sum(c(MPI,PR,SIFGPP)) + evalonly
        mean_global_sum[sif_sub] <-  mean_gs 

        # arrange zonal data
        # - area integrated annual mean values using the year range specified in yr_mean 
        ymai    <- areai 
        zmlat   <- as.data.frame.table(tapply(ymai,edf$lat,sum))
        zmlon   <- as.data.frame.table(tapply(ymai,edf$lon,sum))
      
        zmlat[,1]   <- as.numeric(as.character(zmlat[,1]))
        zmlon[,1]   <- as.numeric(as.character(zmlon[,1]))
  
        zonal_lat[sif_sub,] <- zmlat
        zonal_lon[sif_sub,] <- zmlon
      } 

    }

    # plotting order
    if(evalonly) {  
     df$run      <- factor(df$run,levels=c('MPI','WangFP','SIF-CASA','scaled-SIF')[which(c(MPI,PR,SIFGPP,SIF))])
     mean_df$run <- factor(mean_df$run,levels=c('MPI','WangFP','SIF-CASA','scaled-SIF')[which(c(MPI,PR,SIFGPP,SIF))])
    }
  }
 
  if(EOF) {
    # open climate data 
    for(clim_var in c('prc','tmp','swr')) {
           
      print(clim_var)
      index  <- year-styr+3
      mydata <- open(clim_var,c('lat','lon',paste('y',styr:endyr,sep='')),wdpath)
 
      # subset to region
      lims <- region(pregion)
      mds  <- subset(mydata, lon>(lims[1]-mod_res[1]/2))
      mds  <- subset(mds,    lon<(lims[2]+mod_res[1]/2))
      mds  <- subset(mds,    lat>(lims[3]-mod_res[2]/2))
      mds  <- subset(mds,    lat<(lims[4]+mod_res[2]/2))      
      mydata <- mds
   
      # apply mask
      if(!is.null(mask))  mydata <- subset(mydata,!is.na(mask_v))
      
      # calculate mean across years
      ym             <- apply(as.matrix(mydata[,(yr_mean[1]-styr+3):(yr_mean[2]-styr+3)]),1,mean)
      ym             <- cbind(mydata[,1:2],ym,clim_var,mean(yr_mean))
      names(ym)[3:5] <- c('plotdata','run','year')
  
      # normalise 
      if(!is.null(norm)) ym <- fnorm(ym,norm) 

      climdf <- if(clim_var=='prc') ym else rbind(climdf,ym)
    }

    # make data matrix: rows - lat lon vector, cols - run
    mean_mat <- matrix(as.vector(c(mean_df$plotdata,climdf$plotdata)),nrow=length(df1[,1]))
    # eigen values of covariance matrix
    #eig <- eigen(cov(mean_mat))
    
    mean_mat_noNA <- mean_mat
    mean_mat_noNA[which(is.na(mean_mat))] <- 0
 
    # use prcomp to calculate principle components
    pca <- prcomp(mean_mat_noNA, scale = F)
    # Eigenvalues
    eig <- (pca$sdev)^2
    # Proportional variances in percentage
    variance <- eig*100/sum(eig)
    # Cumulative variance
    cumvar <- cumsum(variance)
    pca$eig <- data.frame(eig = eig, variance = variance, cumvariance = cumvar)

    # add data to output
    pca$lat       <- df1$lat
    pca$lon       <- df1$lon
    pca$data      <- mean_mat
    pca$data_noNA <- mean_mat_noNA
    pca$runid     <- c(levels(mean_df$run),c('prc','tmp','swr'))

    # output PCA
    print('write PCA .rds',quote=F)
    setwd(owd)
    saveRDS(pca,paste(of,'PCA.rds',sep='_'))

  }

  # adjust plotting variables for single vs multiple simulations
  if(is.null(dim(global_sum))) {
    nsims    <- 1
    globsum  <- global_sum[index-2] 
    tglobsum <- global_sum
    if(trend_norm) tglobsum <- tglobsum - tglobsum[1]
    lt       <- length(global_sum)
  } else { 
    nsims    <- dim(global_sum)[1]
    globsum  <- global_sum[,index-2] 
    tglobsum <- t(global_sum[1:nsims,]) 
    if(trend_norm) tglobsum <- apply(tglobsum,2,function(v) v - v[1])
    lt       <- dim(global_sum)[2]
  }
  gsi <- if(is.null(icon)) 1:length(mean_global_sum) else icon[[1]]
  print(gsi) 
 
  # make difference
  if(!is.null(diff)) {
    print('') 
    print(paste('diff selected:',diff))
    # the run subsets must be identical in their lat and lon order
    df$plotdata      <- df$plotdata - df$plotdata[df$run==diff]
    mean_df$plotdata <- mean_df$plotdata - mean_df$plotdata[mean_df$run==diff]
#    rmse    <- mean_df %>%
#                 group_by(run) %>%
#                 summarise(rmse = function(x) mean(abs(x),na.rm=T)
    rmse    <- as.data.frame.table(tapply(mean_df$plotdata,mean_df$run,function(x) mean(abs(x),na.rm=T) ))
    df      <- subset(df,run!=diff)
    mean_df <- subset(mean_df,run!=diff)
  }

  # map plots
  # individual year plot
  print('make first plot')
  p1  <- plot_map_lattice(df,get(vars[v]),
                          pregion,gs=globsum[gsi],norm=norm,diff=diff,mask=mask,labcex=labcex,stripprint,
                          layout=c(cs,rs),skip=skip,index.cond=icon,main=list(main,cex=labcex))  

  # mean across years plot
  print('make second plot')
  p1x <- plot_map_lattice(mean_df,get(vars[v]),
                          pregion,gs=mean_global_sum[gsi],norm=norm,diff=diff,mask=mask,labcex=labcex,stripprint,
                          layout=c(cs,rs),skip=skip,index.cond=icon,main=list(main,cex=labcex))
  
  # make plots
  mpars <- trellis.par.get()
  mpars$strip.background$col <- c('grey90','grey80')
  mpars$fontsize$text        <- 20  
  mpars$panel.background$col <- 'white'
  res  <- 400

  setwd(owd)
  prefix     <- if(SIF)                paste(of,pregion,'SIF',sep='_')   else paste(of,pregion,sep='_')
  prefix     <- if(MPI)                paste(prefix,'MPI-GPP',sep='_')   else paste(of,pregion,sep='_')
  prefix     <- if(!is.null(norm))     paste(prefix,'norm',sep='_')      else prefix
  prefix     <- if(!is.null(diff))     paste(prefix,'_diff',diff,sep='') else prefix
  prefix     <- if(!is.null(mask))     paste(prefix,'mask:',mask,mask_perc,'perc',sep='_') else prefix
  print(paste('Prefix:',prefix))
  ofile      <- paste(prefix,'_',vars[v],sep='')
  ofile_mean <- paste(prefix,'_mean_',vars[v],sep='')
  
  if(!is.null(pyear)) {
    if(plotype=='pdf') pdf(paste(ofile,'pdf',sep='.'),width=,width,height=height,pointsize=pointsize)
    else               png(paste(ofile,'png',sep='.'),width=4*res,height=3*res,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p1)
    dev.off()  
  } 

  if(!trendy) { 
    if(plotype=='pdf') pdf(paste(ofile_mean,'pdf',sep='.'),width=width,height=height,pointsize=pointsize)
    else               png(paste(ofile_mean,'png',sep='.'),width=4*res,height=3*res,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p1x)
    dev.off()  
  }
  
  # make table
  if(table){
    if(v==1) out_table <- data.frame(sim=levels(df$run),mean_global_sum[gsi])
    else     out_table <- cbind(out_table,mean_global_sum[gsi])
    if(v==length(vars)) {
      print('')
      print('Printing table')
      names(out_table)[2:length(out_table)] <- vars
      write.csv(out_table,paste(prefix,'.csv',sep=''),quote=F,row.names=F)
    }
  }

  if(!is.null(diff)) write.csv(rmse,paste(prefix,'_rmse.csv',sep=''),quote=F,row.names=F)

  # make trend & zonal plots
  if(trends&!is.null(zonal_lon)){
    setwd(owd)
    print('')
    print('printing trends & zonal plots')

    # if requested break up below plots into panels
    zon_panels <- 1
    if(!is.null(zon_panel_names)) zon_panels <- length(zon_panel_names)
    else zon_panel_names <- lab$name

    # rename simulations
    if(is.null(slabels)) slabels <- sim[sia]

    # axis label
    axis_unit <- if(lab$gsum)     lab$sunit else lab$unit
    axis_unit <- if(lab$printsum) axis_unit else "'-'"
    axis_lab  <- eval(parse(text=paste('expression(',lab$name,"  (",axis_unit,")",')',sep='')))
 
    # key 
    #key.line  <- list(x=0.1,y=0.95,size=2,between=0.5,border=T,padding.text=0.5*labcex,
    #                  text=list(slabels,cex=labcex*0.2),lines=list(lty=lty[sia],col=col_trend[sia],lwd=labcex*0.7))
    key.line  <- list(space='top',size=2,between=0.5,border=T,padding.text=0.5*labcex,
                      #text=list(slabels,cex=labcex*0.2),lines=list(lty=lty[sia],col=col_trend[sia],lwd=labcex*0.7))
                      text=list(slabels,cex=labcex*0.2),lines=list(lty=lty,col=col_trend,lwd=labcex*0.7))
  
    print('make trend plot')
#    print(dim(tglobsum))
#    print(length(tglobsum))
#    print(c(styr,endyr,nsims,lt))
#    print(class(tglobsum))
#    print(head(tglobsum))
#    print(tglobsum)
#    print(col_trend)
#    print('')
    p2 <-
    xyplot(as.vector(tglobsum)
           ~rep(styr:endyr,nsims),
           groups=(rep(1:nsims,each=lt)),
           #ylim=c(-15,10),
           xlab=list('year',cex=labcex*0.35),
           ylab=list(axis_lab,cex=labcex*0.35),
           #type='l',lwd=labcex*1.1,lty=lty[sia],col=col_trend[sia],
           type='l',lwd=labcex*1.1,lty=lty,col=col_trend,
           scales=list(alternating=F,tck=c(-0.5,0),cex=labcex*0.35),
           strip=F,
           key=key.line
           )

    print('make zonal lat plot')
    p3 <-
    xyplot(rep(zonal_lat[1,],nsims)~as.vector(t(zonal_lat[2:(nsims+1),]))|rep(zon_panel_names,each=dim(zonal_lat)[2]*nsims/zon_panels),
           groups=(rep(1:nsims/zon_panels,each=dim(zonal_lat)[2])),
           #ylab=list('latitude',cex=labcex*0.35),type='l',lwd=labcex*1.1,lty=lty[sia],col=col_trend[sia],
           ylab=list('latitude',cex=labcex*0.35),type='l',lwd=labcex*1.1,lty=lty,col=col_trend,
           xlab=list(axis_lab,cex=labcex*0.35),
           scales=list(alternating=F,tck=c(-0.5,0),cex=labcex*0.35),
           strip=F,
           key=key.line
           )

    print('make zonal lon plot')
    p4 <-
    xyplot(t(zonal_lon[2:(nsims+1),])
           ~rep(zonal_lon[1,],nsims)|rep(zon_panel_names,each=dim(zonal_lon)[2]*nsims/zon_panels),
           groups=(rep(1:nsims,each=dim(zonal_lon)[2])),
           #xlab=list('longitude',cex=labcex*0.35),type='l',lwd=labcex*1.1,lty=lty[sia],col=col_trend[sia],
           xlab=list('longitude',cex=labcex*0.35),type='l',lwd=labcex*1.1,lty=lty,col=col_trend,
           ylab=list(axis_lab,cex=labcex*0.35),
           scales=list(alternating=F,tck=c(-0.5,0),cex=labcex*0.35),
           strip=F,
           key=key.line
           )
    
    # make plots
    ofile <- paste(of,'_',pregion,'_trend_',vars[v],sep='')
    if(plotype=='pdf') pdf(paste(ofile,'pdf',sep='.'),width=width1,height=width1,pointsize=pointsize)
    else               png(paste(ofile,'png',sep='.'),width=4*res,height=3*res,pointsize=pointsize,bg=backg)
    trellis.par.set(mpars)
    print(p2)
    dev.off()
    
    ofile <- paste(of,'_',pregion,'_zonallat_',vars[v],sep='')
    if(plotype=='pdf') pdf(paste(ofile,'pdf',sep='.'),width=width1,height=width1,pointsize=pointsize)
    else               png(paste(ofile,'png',sep='.'),width=4*res,height=3*res,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p3)
    dev.off()
    
    ofile <- paste(of,'_',pregion,'_zonallon_',vars[v],sep='')
    if(plotype=='pdf') pdf(paste(ofile,'pdf',sep='.'),width=width1,height=width1,pointsize=pointsize)
    else               png(paste(ofile,'png',sep='.'),width=4*res,height=3*res,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p4)
    dev.off()

  # end plot trends loop
  }

  if(trend_out) {
    colnames(tglobsum) <- slabels 
    saveRDS(tglobsum,paste(prefix,'.RDS',sep=''))
  }

  if(scatter) {
    if(v==1) mean_df_prev <- mean_df
    else {
      # mean across years plot
      print('make scatter plot')
      ps  <- xyplot(mean_df$plotdata~mean_df_prev$plotdata|mean_df$run,
                    ylab=list(eval(parse(text=paste('expression(',get(vars[v])$name,"   (",get(vars[v])$unit,")",')',sep=''))),
                              cex=0.5*labcex),
                    xlab=list(eval(parse(text=paste('expression(',get(vars[v-1])$name,"   (",get(vars[v-1])$unit,")",')',sep=''))),
                              cex=0.5*labcex),
                    as.table=T)
      # p1x <- plot_map_lattice(mean_df,get(vars[v]),
      #                        pregion,gs=mean_global_sum[gsi],norm=norm,diff=diff,mask=mask,labcex=labcex,stripprint,
      #                        layout=c(cs,rs),skip=skip,index.cond=icon,main=list(main,cex=labcex))
      ofile <- paste(of,'_',pregion,'_scatter_',vars[v],'-',vars[v-1],sep='')
      if(plotype=='pdf') pdf(paste(ofile,'pdf',sep='.',width=width,height=width,pointsize=pointsize))
      else               png(paste(ofile,'png',sep='.'))#,width=4*res,height=3*res,pointsize=28,bg=backg)
      trellis.par.set(mpars)
      print(ps)
      dev.off()

    }
  }


  print('Warnings:')
  print(warnings())

# end variable loop
}

print('end plot script')


