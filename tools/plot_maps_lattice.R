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



##################################
###user defined inputs

#directory paths
date    <- '151210'
dir     <- '/home/alp/models/SDGVM/'
rdir    <- 'run'
edir    <- 'eval_data'

# project directory
project <- 'vcmax'

# simulations
sim     <- c('original','orig_N','Kattge','Kattge_oxisol','Maire','vBodegom_env','vBodegom_mean','Walker_N','Walker_NP','Woodward_94','Woodward_95')

# simulations index array (which simulations to plot from the vector 'sim', when 'sim' is arranged in alphabetical order)
sia     <- NULL

# output file naming prefix
of      <- NULL 

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
year    <- 2012

# years over which to take a mean value, must be a two element vector 
yr_mean <- c(2001,2010)

# plot trends across whole timeseries
trends  <- T 

# normalise the data to the norm*100 %ile
norm    <- NULL 

# make a difference plot, takes the value of the simulation id string as the base for the differences
diff    <- NULL

# load SIF data
SIF     <- F

# line colour & type on trend and zonal plots, correspond to the labels in 'sim' when 'sim' is in alphabetical order
col_trend <- topo.colors(10)
lty       <- NULL 

# backgroud for plots
backg     <- 'transparent'
backg     <- 'white'

# plotting variables (these all have an associated list in 'params_map_plot.R', to add variables simply add a list in 'params_map_plot.R' with the same name as the variable).
vars <- c('npp','gpp','nbp',
          'anlfn','antlfn',        
          'evt','trn','scn','sresp','presp','mgresp',
          'lai','anvcmax','anjmax','biot','kg_beta',
          'swr','qtotal','tmp','prc','swc','field_capacity','wilting_point',
          'cov_C3','cov_C4','cov_C3crop','cov_C4crop','cov_Dc_Bl','cov_Dc_Nl','cov_Ev_Bl','cov_Ev_Nl','cov_BARE')

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
if(is.null(sia)) sia <- 1:length(sim)
if(is.null(lty)) lty <- rep(1:length(sia))
if(is.null(of))  of  <- project

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
print('',quote=F)


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
  sif_df <- read.table('SIF_l3_2007t2012mean_screened_1x1.dat',header=F)
  names(sif_df) <- c('lat','lon','plotdata')
  sif_df$run    <- 'SIF' 
  sif_df$year   <- 2009.5  
}

# variable plotting loop
vars <- vars[via]
# change this to an mclapply -  so far there are a few too many arguments
#mclapply(1:length(vars),,)

for( v in 1:length(vars) ){
  print('',quote=F)
  print(vars[v],quote=F)
  
  lab <- get(vars[v])

  # layout parameters
  if(!is.null(rs)){
   rs <- ceiling((length(sia)*length(year))^0.5)
   if((length(sia)*length(year)) < 4) rs <- length(sia)*length(year)
   if(pregion=='tropics')             rs <- rs + 2
   # if(rs==1)                          rs <- 2
  }
  if(!is.null(cs)) cs <- ceiling(length(sia)*length(year)/rs)
  
  # simulations loop
  for(m in sia){    
    
    #open data 
    wdpath <- paste(dir,rdir,project,sim[m],'output/',sep='/')
    print(wdpath)
    if(substr(vars[v],1,3)=='cov') {
      index  <- 3
      mydata <- open(vars[v],c('lat','lon','cov'),wdpath)
      rs     <- 3 
      cs     <- 1
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
  
    # calculate mean across years
    ym             <- apply(as.matrix(mydata[,(yr_mean[1]-styr+3):(yr_mean[2]-styr+3)]),1,mean)
    ym             <- cbind(mydata[,1:2],ym,sim[m],mean(yr_mean))
    names(ym)[3:5] <- c('plotdata','run','year')
 
    # calculate area integrated values
    if(lab$gsum|lab$gmean){

      # area integrate
      areai <- area_integrate(mydata[,1:2],mydata[,3:length(mydata)],mod_res,lab$gmean)
      gs    <- apply(as.matrix(areai),2,sum,na.rm=T)      
      global_sum      <- if(m==sia[1]) gs else rbind(global_sum,gs)
      mean_global_sum <- if(m==sia[1]) mean(gs[(yr_mean[1]-styr+1):(yr_mean[2]-styr+1)]) else 
                         c(mean_global_sum,mean(gs[(yr_mean[1]-styr+1):(yr_mean[2]-styr+1)]))      

      # arrange zonal data 
      if(substr(vars[v],1,3)!='cov') {
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
    dfl <- lapply(index,stack,df=mydata) 
    df1 <- rbind.fill(dfl)
    df1$run <- sim[m]    
  
    # normalise data
    if(!is.null(norm)) {
      # centre with mean and scale by standard deviation
      if(norm=='sd'){ 
      
        df1$plotdata <- ( df1$plotdata - mean(df1$plotdata,na.rm=T) )  / (var(df1$plotdata,na.rm=T)^0.5) 
        ym$plotdata  <- ( ym$plotdata  - mean(ym$plotdata,na.rm=T) )   / (var(ym$plotdata,na.rm=T)^0.5)  
      
      } else {
      
        if(any(df1$plotdata<0)) {
          df1$plotdata <- df1$plotdata / max( abs(quantile(df1$plotdata,1-norm,type=8)) , abs(quantile(df1$plotdata,norm,type=8)) )
          ym$plotdata  <- ym$plotdata  / max( abs(quantile(ym$plotdata,1-norm,type=8))  , abs(quantile(ym$plotdata,norm,type=8)) )
        } else {
          df1$plotdata <- df1$plotdata / quantile(df1$plotdata,norm,type=8)
          ym$plotdata  <- ym$plotdata  / quantile(ym$plotdata,norm,type=8)
        }
      
      }
    } 
  
    # make dataframe with all simulations
    df      <- if(m==sia[1]) df1 else rbind(df,df1)
    mean_df <- if(m==sia[1]) ym else rbind(mean_df,ym)      

  # end simulations loop
  }

  # append datasets with eval data
  if(SIF) {
    # restrict SIF dataset to pixels that occur in the model dataset only
    # hijacking the now unused df1 dataframe from the simulations loop
    df1$plotdata <- NA
    df1$run      <- 'SIF'
    df1$plotdata <- sif_df$plotdata[match(paste(df1$lat,df1$lon),paste(sif_df$lat,sif_df$lon))]
    sif_df 	 <- df1
    if(!is.null(norm)) {
      if(norm=='sd') sif_df$plotdata <- ( sif_df$plotdata - mean(sif_df$plotdata,na.rm=T) ) / (var(sif_df$plotdata,na.rm=T)^0.5)
      else           sif_df$plotdata <- sif_df$plotdata  / max( abs(quantile(sif_df$plotdata,1-norm,type=8,na.rm=T)) , abs(quantile(sif_df$plotdata,norm,type=8,na.rm=T)) )
    }
    df      <- rbind(df,sif_df)
    mean_df <- rbind(mean_df,sif_df)
  }
 
  # adjust plotting variables for single vs multiple simulations
  if(is.null(dim(global_sum))==1) {
    nsims    <- 1
    globsum  <- global_sum[index-2] 
    tglobsum <- global_sum
    lt       <- length(global_sum)
  } else { 
    nsims    <- dim(global_sum)[1]
    globsum  <- global_sum[,index-2] 
    tglobsum <- t(global_sum[1:nsims,]) 
    lt       <- dim(global_sum)[2]
  }
  gsi <- if(is.null(icon)) 1:length(globsum) else icon[[1]]
  
  # make difference
  if(!is.null(diff)) {
    # the run subsets must be identical in their lat and lon order
    df$plotdata      <- df$plotdata - df$plotdata[df$run==diff]
    mean_df$plotdata <- mean_df$plotdata - mean_df$plotdata[mean_df$run==diff]
  }

  # map plots
  # individual year plot
  print('make first plot')
  p1  <- plot_map_lattice(df,get(vars[v]),
                          pregion,gs=globsum[gsi],norm=norm,diff=diff,
                          layout=c(cs,rs),skip=skip,index.cond=icon)
  
  # mean across years plot
  print('make second plot')
  p1x <- plot_map_lattice(mean_df,get(vars[v]),
                          pregion,gs=mean_global_sum[gsi],norm=norm,diff=diff,
                          layout=c(cs,rs),skip=skip,index.cond=icon)

  # make plots
  mpars <- trellis.par.get()
  mpars$strip.background$col <- c('grey90','grey80')
  mpars$fontsize$text        <- 20  
  mpars$panel.background$col <- 'white'
  res  <- 400

  setwd(owd)
  prefix     <- if(!is.null(norm)) paste(of,pregion,'norm',sep='_') else paste(of,pregion,sep='_')
  prefix     <- if(!is.null(diff)) paste(prefix,'diff',sep='_')     else prefix
  print(paste('Prefix:',prefix))
  ofile      <- paste(prefix,'_',vars[v],'.png',sep='')
  ofile_mean <- paste(prefix,'_mean_',vars[v],'.png',sep='')
  pdfofile   <- paste(prefix,'_',vars[v],'.pdf',sep='')
  
  png(ofile,width=4*res,height=3*res,pointsize=28,bg=backg)
  trellis.par.set(mpars)
  print(p1)
  dev.off()  
  
  png(ofile_mean,width=4*res,height=3*res,pointsize=28,bg=backg)
  trellis.par.set(mpars)
  print(p1x)
  dev.off()  
  
  #pdf(pdfofile,width=18,height=12,pointsize=28)
  #trellis.par.set(mpars)
  #print(p1)
  #dev.off()  
  
  # trend & zonal plots
  if(trends&!is.null(zonal_lon)){
    setwd(owd)
    print('printing trends & zonal plots')

    print('make trend plot')
    p2 <-
    xyplot(tglobsum
           ~rep(styr:endyr,nsims),
           groups=(rep(1:nsims,each=lt)),
           xlab='year',
           ylab=paste(lab$name,if(lab$gsum) lab$sunit else lab$unit),
           type='l',lwd=3,lty=lty[sia],col=col_trend[sia],
           scales=list(alternating=F,tck=c(-0.5,0)),
           key=list(x=0,y=0.98,text=list(sim[sia]),lines=list(lty=lty[sia],col=col_trend[sia]))
           )
    print('make zonal lat plot')
    p3 <-
    xyplot(t(zonal_lat[2:(nsims+1),])
           ~rep(zonal_lat[1,],nsims),
           groups=(rep(1:nsims,each=dim(zonal_lat)[2])),
           xlab='latitude',type='l',lwd=3,lty=lty[sia],col=col_trend[sia],
           ylab=paste(lab$name,if(lab$gsum) lab$sunit else lab$unit),
           scales=list(alternating=F,tck=c(-0.5,0)),
           key=list(x=0,y=0.98,text=list(sim[sia]),lines=list(lty=lty[sia],col=col_trend[sia]),text=list(as.character(round(globsum,r))))
           )
    print('make zonal lon plot')
    p4 <-
    xyplot(t(zonal_lon[2:(nsims+1),])
           ~rep(zonal_lon[1,],nsims),
           groups=(rep(1:nsims,each=dim(zonal_lon)[2])),
           xlab='longitude',type='l',lwd=3,lty=lty[sia],col=col_trend[sia],
           ylab=paste(lab$name,if(lab$gsum) lab$sunit else lab$unit),
           scales=list(alternating=F,tck=c(-0.5,0)),
           key=list(x=0.65,y=0.98,text=list(sim[sia]),lines=list(lty=lty[sia],col=col_trend[sia]),text=list(as.character(round(globsum,r))))
           )
    
    # make plots
    ofile <- paste(of,'_',pregion,'_trend_',vars[v],'.png',sep='')
    png(ofile,width=1200,height=900,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p2)
    dev.off()
    
    ofile <- paste(of,'_',pregion,'_zonallat_',vars[v],'.png',sep='')
    png(ofile,width=1200,height=900,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p3)
    dev.off()
    
    ofile <- paste(of,'_',pregion,'_zonallon_',vars[v],'.png',sep='')
    png(ofile,width=1200,height=900,pointsize=28,bg=backg)
    trellis.par.set(mpars)
    print(p4)
    dev.off()
  }  
  print('Warnings:')
  print(warnings())
  print('end plot script')
}



