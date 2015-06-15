####################################
#
# global plotting function
#
# AWalker
# Feb 2012
#
###################################

rm(list=ls())

setwd('/mnt/disk2/script_library/R/map_tools/')
# source('plotting_functions.R')
source('map_plot.R')

library(lattice)

##################################
###user defined inputs

#directory paths
date    <- '150223'
wd      <- '/mnt/disk2/Research_Projects/leaf_trait_SDGVM/'
trends  <- F

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
year <- 2010


# variables and assigned qualities
mod_res <- c(3.75,2.5)
pregion <- 'global'
styr    <- 1901
endyr   <- 2010

vars<-c('npp','gpp','nbp',
        'anlfn','antlfn',        
        'evt','trn','scn','sresp','presp','mgresp',
        'lai','anvcmax','anjmax','biot','kg_beta','swr','qtotal',
        'tmp','prc','swc','field_capacity','wilting_point',
        'cov_C3','cov_C4','cov_C3crop','cov_C4crop','cov_Dc_Bl','cov_Dc_Nl','cov_Ev_Bl','cov_Ev_Nl','cov_BARE')

#variable index array
via <- 1:length(vars)-2
# via <- c(1:4,7,8,12,15:20)
# via <- 19
# via <- c((length(vars)-8):length(vars))
via <- 1:21
via <- 2

#scale max min and number of divisions on scale
scaleminv<-c(0,0,-400,
             0,0,
             0,0,0,0,0,0,0,
             0,0,0,0,0,0,
             -10,0,0,0,0,
             rep(0,9))

scalemaxv<-c(1600,4000,400,
             6,3,
             1500,1000,40000,1500,400,1600,
             10,70,130,25000,1,260,20000,
             30,4000,700,0.45,0.2,
             rep(1,9))

vname<-c('NPP','GPP','NBP',
         'Leaf N','Topleaf N',
         'Evapotranspiration','Transpiration',
         'Soil Carbon','Soil Respiration','Canopy Respiration','Plant Respiration',
         'LAI',
         'V','J','Vegetation Biomass','soil water stress scalar','SWradiation','PAR',
         'Air Temperature','Precipitation','Soil Water Content','Field Capacity','Wilting Point',
         'cov_C3','cov_C4','cov_C3crop','cov_C4crop','cov_Dc_Bl','cov_Dc_Nl','cov_Ev_Bl','cov_Ev_Nl','cov_BARE')

vsub<-c(rep('',12),'cmax','max',rep('',18))
  
sunit<-c(rep('gm^-2',5),
         rep('kgm^-2',2),
         rep('gm^-2',4),
         'm^2*m^-2',
         rep('mu*mol*" "*m^-2*" "*s^-1',2),
         'gm^-2','(unitless)','unknown','mol*m^-2','.^o*C','mm','mm','(unitless)','(unitless)',
         rep('(unitless)',9))

colo<-c(1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,rep(1,9))

#add global sum or mean
gsum  <- c(rep(T,4),F,rep(T,6),rep(F,3),T,F,T,F,F,F,F,F,F,rep(F,9))
gmean <- c(rep(F,15),T,T,F,T,T,T,T,T,rep(T,9))
unit  <- c(rep('Pg',5),rep('Eg',2),rep('Pg',4),rep('',3),'Pg','','?','?','oC','mm','mm','','',rep('',9))



#setup colour schemes  
white.red    <- colorRampPalette(c('red3','lightgoldenrod'),space='Lab',bias=10)
yellow.green <- colorRampPalette(c('lightgoldenrod','darkgreen'),space='Lab',bias=10)
# green.blue<-colorRampPalette(c('darkgreen','midnightblue'),space='Lab',bias=10)
# poscolours<-c(yellow.green(9),green.blue(5)[2:5])
negcolours   <- c(white.red(7)[1:6],'lightgoldenrod',yellow.green(7)[2:7])
# drycolours<-white.red(13)
# colv<-list(poscolours,negcolours,c('blue','red'),drycolours)
# ncolours<-13

#setup colours
colv <- list(rev(topo.colors(10)),negcolours)
ncol <- 10
sfactor <- scalemaxv/ncol

###############################
###start program

#decimal places for scale legend
r<-2

# layout parameters
rs  <- ceiling((length(sia)*length(year))^0.5)
if((length(sia)*length(year)) < 4) rs <- length(sia)*length(year)
if(rs==1) rs <- 2
cs  <- ceiling(length(sia)*length(year)/rs)

# or manually set columns (cs) and rows (rs)
# cs <- 2 
# rs <- 3

trs <- rs+1
ifelse(cs==2,he<-0.8,he<-0.5)
end <- rs*cs+1
vec <- c(1:(rs*cs),rep(end,cs))
blank <- ( cs*rs - length(sim)*length(year) )

# variable loop
v   <- length(vars)
# via <- 1
for(v in via){
    
  #open plot file & layout
  setwd(owd)
  ofile <- paste(of,'_',pregion,'_',year[1],'_',vars[v],'.png',sep='')
  png(ofile,width=1200,height=900,pointsize=28)
  mar <- c(1,0,1,0)
  par(oma=c(0,1,1,1),mar=mar,cex.axis=1,tck=0.02,lwd=1,mgp=c(2,0.1,0))  
  layout(t(matrix(vec,cs,trs)), heights=c(rep(0.9,rs),0.4))
  
  # plotting loop
  m <- 1
  for(m in sia){    
    
    for(y in 1:length(year)){
      
      #open data 
      wdpath  <- paste(wd,'simulations',project,sim[m],'output/',sep='/')
      print(wdpath)
      #mydata column to plot
      if(substr(vars[v],1,3)=='cov') {
        index   <- 3
        mydata  <- open(vars[v],c('lat','lon','cov'),wdpath)
      } else {
        index   <- year[y]-styr+3
        mydata  <- open(vars[v],c('lat','lon',paste('y',styr:endyr,sep='')),wdpath)
      }
      
      #get global area integrated values
      if((sum(gsum[v]+gmean[v])>0)&y==1){
        gs <- area_integrate(mydata[,1:2],mydata[,3:length(mydata)],mod_res,gmean[v])
        if(m==sia[1]) global_sum <- gs
        else     global_sum <- rbind(global_sum,gs)    
      }

      #add polygons and area integration to data
      mydata  <- grid_calc(mydata,gsum[v],index,mod_res)
      
      #plot
      plot_map(mydata,
               index=index,scalemin=scaleminv[v],scalemax=scalemaxv[v],sfactor[v],
               vname[v],vsub[v],sunit[v],colo[v],gsum[v],unit[v],
               colv,pregion,sim[m],year[y],mar=mar)
    }
  }

  #fill blank plot layout matrix so scale can go at bottom of plot
  if(blank>0){
    for( b in 1:blank ){
      plot(1~1,axes=F,xlab='',pch=NA)
    }  
  }
  
  #plot scale bar
  plot_scale(scaleminv[v],scalemaxv[v],colv[[colo[v]]],r,sfactor[v],sunit[v],vname[v],vsub[v])
  dev.off()  
  
  #plot trends
  if(trends){
    setwd(owd)
    lty <- 1
    if(sum(gsum[v]+gmean[v])>0){
      ofile <- paste(of,'_trend_',vars[v],'.png',sep='')
      
      nsims <- dim(global_sum)[1]
      png(ofile,width=1200,height=900,pointsize=28)
      plot1 <-
        xyplot(t(global_sum[1:nsims,])
               ~rep(styr:endyr,nsims),
               groups=(rep(1:nsims,each=dim(global_sum)[2])),
               xlab='year',ylab='',type='l',lwd=3,lty=lty,col=c('black','blue','red'),
               scales=list(alternating=F,tck=c(-0.5,0)),
               key=list(x=0,y=0.98,text=list(sim[sia]),lines=list(lty=lty,col=c('black','blue','red')[1:nsims]))
        )
      print(plot1)
      dev.off()
    }  
  }
}



