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
date    <- '150615'
wd      <- '/home/alp/models/SDGVM/'
rdir    <- 'run'

# project directory
project <- 'vcmax'

# simulations
sim     <- c('original','orig_N','Kattge','Kattge_oxisol','Maire','vBodegom_env','vBodegom_mean','Walker_N','Walker_NP','Woodward_94','Woodward_95')

# simulations index array (which simulations to plot)
sia     <- 1:length(sim)

# output file naming prefix
of      <- project

# decimal places for scale legend
r <- 2

# map data resolution  
deg1    <- T

# data start year
styr    <- 1901

# data end year
endyr   <- 2010

# plotting region
pregion <- 'global'

# plot year (can be a vector of years)
year    <- 2012

# plot trends across whole timeseries
trends  <- T

# plotting variables (these all have an associated list in 
vars <- c('npp','gpp','nbp',
          'anlfn','antlfn',        
          'evt','trn','scn','sresp','presp','mgresp',
          'lai','anvcmax','anjmax','biot','kg_beta','swr','qtotal',
          'tmp','prc','swc','field_capacity','wilting_point',
          'cov_C3','cov_C4','cov_C3crop','cov_C4crop','cov_Dc_Bl','cov_Dc_Nl','cov_Ev_Bl','cov_Ev_Nl','cov_BARE')

# variable index array (which variables to plot)
via <- 1:length(vars)
via <- 2


###############################
### start program

setwd(paste(wd,'src/sdgvm/tools/',sep='/'))
source('params_map_plot.R')
source('functions_map_plot.R')

         mod_res <- c(3.75,2.5)
if(deg1) mod_res <- c(1,1)

# create output directory
owd     <- paste(wd,rdir,project,'results',date,sep='/')
if(!file.exists(owd)) dir.create(owd)



# variable plotting loop
for(v in via){
  lab <- get(vars[v])

  # layout parameters
  rs  <- ceiling((length(sia)*length(year))^0.5)
  if((length(sia)*length(year)) < 4) rs <- length(sia)*length(year)
  # if(rs==1) rs <- 2
  cs  <- ceiling(length(sia)*length(year)/rs)
  
  # or manually set columns (cs) and rows (rs)
  # cs <- 2 
  # rs <- 3
  
  # data loop
  for(m in sia){    
    #open data 
    wdpath  <- paste(wd,rdir,project,sim[m],'output/',sep='/')
    print(wdpath)
    
    #mydata column to plot
    if(substr(vars[v],1,3)=='cov') {
      index   <- 3
      mydata  <- open(vars[v],c('lat','lon','cov'),wdpath)
      rs <- 3 
      cs <- 1
    } else {
      index   <- year-styr+3
      mydata  <- open(vars[v],c('lat','lon',paste('y',styr:endyr,sep='')),wdpath)
    }
    
    #get global area integrated values
    if(lab$gsum|lab$gmean){
      gs <- area_integrate(mydata[,1:2],mydata[,3:length(mydata)],mod_res,lab$gmean)
      if(m==sia[1]) { global_sum <- gs 
      } else          global_sum <- rbind(global_sum,gs)    
    } else global_sum <- NULL

    # arrange data for lattice plotting
    stack <- function(col,df){
      df1 <- df[,c(1:2,col)]
      df1$year <- col+styr-3
      names(df1)[3] <- 'plotdata'
      df1
    }
    dfl <- lapply(index,stack,df=mydata) 
    df1 <- rbind.fill(dfl)
    df1$run <- sim[m]    
    if(m==sia[1]) {
      df <- df1
    } else df <- rbind(df,df1)
    
  }
  
  # make plot
  p1 <- plot_map_lattice(df,index,get(vars[v]),
                         pregion,sim[m],mar=mar,gs=global_sum,layout=c(cs,rs))

  # plots
  mpars <- trellis.par.get()
  mpars$strip.background$col <- c('grey90','grey80')

  setwd(owd)
  ofile    <- paste(of,'_',pregion,'_',vars[v],'.png',sep='')
  pdfofile <- paste(of,'_',pregion,'_',vars[v],'.pdf',sep='')
  
  png(ofile,width=1200,height=900,pointsize=28)
  trellis.par.set(mpars)
  print(p1)
  dev.off()  
  pdf(pdfofile,width=18,height=12,pointsize=28)
  trellis.par.set(mpars)
  print(p1)
  dev.off()  
  
  
  #plot trends
  if(trends&(substr(vars[v],1,3)!='cov')&(lab$gsum|lab$gmean)){
    setwd(owd)
    lty <- 1
#     if(sum(gsum[v]+gmean[v])>0){
    if(sum(lab$gsum+lab$gmean)>0){
        ofile <- paste(of,'_trend_',vars[v],'.png',sep='')
      
      nsims <- dim(global_sum)[1]
      png(ofile,width=1200,height=900,pointsize=28)
      p2 <-
        xyplot(t(global_sum[1:nsims,])
               ~rep(styr:endyr,nsims),
               groups=(rep(1:nsims,each=dim(global_sum)[2])),
               xlab='year',ylab='',type='l',lwd=3,lty=lty,col=c('black','blue','red'),
               scales=list(alternating=F,tck=c(-0.5,0)),
               key=list(x=0,y=0.98,text=list(sim[sia]),lines=list(lty=lty,col=c('black','blue','red')[1:nsims]))
        )
      print(p2)
      dev.off()
    }  
  }
}



