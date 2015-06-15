####################################
#
# global plotting function
#
# AWalker
# Feb 2012
#
###################################

rm(list=ls())

setwd('/mnt/disk2/models/SDGVM/tools/')
source('params_map.R')
setwd('/mnt/disk2/script_library/R/map_tools/')
source('map_plot.R')

library(lattice)
library(plyr)

##################################
###user defined inputs

#directory paths
date    <- '150302'
wd      <- '/mnt/disk2/Research_Projects/leaf_trait_SDGVM/'
trends  <- T
deg1    <- T

project <- 'TERRABITES'
sim     <- c('T1_131001_GLC2000','T1_131001_GLOB2009','T1_140129','T1_noC4','T1','T2','T3')
sia     <- c(1:7)
sia     <- 5:7
# sia     <- 3:4
of      <- project

# project <- 'sub-daily_PAR_test'
# sim     <- c('one','five','ten','twenty','thirty')
# sim     <- c('one','five','ten','thirteen','sixteen','twenty','thirty')
# 
# project <- 'SDGVM_light'
# sim     <- c('orig','orig_corr','orig_PAR','orig_Vcmax','varPAR','varPAR_Cl','varPAR_Cln','xPAR','xPAR_Cl','xPAR_Cln')
# of      <- 'SDGVMdev'
# sia <- 1:4
# #of      <- 'SDGVMnew'
# #sia     <- c(8,5,9,6,10,7)
# # sia     <- c(9,6) 
# # sia     <- 5
#
# project <- 'SDGVM_light_140928'
# sim     <- c('orig','orig_corr','orig_PAR','orig_Vcmax','new','new_vm_maire','new_vm_orig1995','new_vm_P','new_C','new_vm_maire_Cl','new_vm_orig1995_Cl','new_vm_P_Cl')
# of      <- 'SDGVMdev'
# sia <- 6
# #sia     <- 5:length(sim)
# #sia     <- 1:4
 
owd     <- paste(wd,'results',project,date,sep='/')
if(!file.exists(owd)) dir.create(owd)


# year/s to plot
year  <- c(1990,2000,2010)
#year  <- c(2009)
styr  <- 1901
endyr <- 2010

# variables
mod_res <- c(3.75,2.5)
if(deg1) mod_res <- c(1,1)
pregion <- 'global'

vars<-c('npp','gpp','nbp',
        'anlfn','antlfn',        
        'evt','trn','scn','sresp','presp','mgresp',
        'lai','anvcmax','anjmax','biot','kg_beta','swr','qtotal',
        'tmp','prc','swc','field_capacity','wilting_point',
        'cov_C3','cov_C4','cov_C3crop','cov_C4crop','cov_Dc_Bl','cov_Dc_Nl','cov_Ev_Bl','cov_Ev_Nl','cov_BARE')

#variable index array
via <- 1:length(vars)
# via <- (length(vars)-8):length(vars)
# via <- 1:21
# via <- 2
# via <- 13:15

###############################
###start program

#decimal places for scale legend
r<-2
# variable loop
v   <- 2
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
  m <- 1
  for(m in sia){    
    #open data 
    wdpath  <- paste(wd,'simulations',project,sim[m],'output/',sep='/')
    if(deg1) wdpath  <- paste(wd,'simulations',project,sim[m],'output_1deg/',sep='/')
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
#   pdf(pdfofile,width=18,height=12,pointsize=28)
#   trellis.par.set(mpars)
#   print(p1)
#   dev.off()  
  
  
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



