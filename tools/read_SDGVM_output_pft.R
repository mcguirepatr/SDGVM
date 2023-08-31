######################
#
# Read SDGVM output
#
# AWalker
# Oct 2013
#
######################

# .libPaths('~/bin/Rlibs')
# library(plyr)
library(parallel)
library(gdata)


# Function to concatenate SDGVM output when run across multiple cores
########################
stich_sdgvm_mp_apply <- function(wd,grids=4,mc=T,
                                 annual=T,monthly=T,daily=T,pft=T,
                                 ...){
  
  # annual files just open and stich together rowise
  # monthly open, convert structure and stich rowise
  # daily open, convert structure and stich rowise
  
  # Initialisation
  ########################
  setwd(wd)
  
  #read & sort output file names
  setwd(paste(wd,'grid1/',sep=''))
  print(paste(wd,'grid1/',sep=''))
  
  files  <- list.files()
  print('FILES 1')
  print(files)
  files  <- files[-which(files=='site_info.dat')]
  files  <- files[-which(files=='simulation.dat')] #PCM
  files  <- files[-which(files=='diag.dat')]

# set the following switch if there is a subset of files that need rerunning
  #file_subset <- T
  file_subset <- F
# set the following switch if there are monthly files that need rerunning
#  monthly <- T
#  monthly <- F
#try only these files: (There has to be at least one annual file)
  #if(file_subset) files <- c('monthly_nep.dat','monthly_srp.dat','leafc.dat')
  #if(file_subset) files <- c('leafc.dat','mgresp.dat','nppstore.dat','presp.dat','rootbio.dat','snn.dat','stemper.dat','swr.dat','trn.dat','monthly_nep.dat','monthly_qdr.dat','monthly_srp.dat','monthly_vcm.dat')
  if(file_subset) files <- c('leafc.dat','mgresp.dat','nppstore.dat','presp.dat','rootbio.dat','snn.dat','stemper.dat','swr.dat','trn.dat')
  #if(file_subset) files <- c('mgresp.dat','presp.dat','snn.dat','swr.dat')
#  if(file_subset) files <- 'mgresp.dat'
  print('FILES 2')
  print(files)
  isubs  <- grep('init',files)
  files  <- files[-isubs]
  print('FILES 3')
  print(files)
  #if(file_subset) files <- c('monthly_qdf.dat','monthly_snw.dat','monthly_trn.dat')
 
  afiles <- files 
  if((!file_subset) || monthly){
    msubs  <- grep('monthly',files)
    dsubs  <- grep('daily',files)
    dfiles <- NULL
    mfiles <- NULL
    if(length(msubs) > 0 ) afiles <- files[-msubs]
    if(length(dsubs) > 0 ) afiles <- afiles[-dsubs]
    if(length(msubs) > 0 ) mfiles <- files[msubs]
    if(length(dsubs) > 0 ) dfiles <- files[dsubs]
    if(monthly) print(mfiles,quote=F)

    mpftsubs  <- grep('[A-Z]',mfiles)
    print(mpftsubs,quote=F)
    dpftsubs  <- grep('[A-Z]',dfiles)
    print(dpftsubs,quote=F)
    mpftfiles <- NULL
    dpftfiles <- NULL
    if(length(mpftsubs) >0 ) mpftfiles <- mfiles[mpftsubs]
    if(length(dpftsubs) > 0 ) dpftfiles <- dfiles[dpftsubs]
    if(length(mpftsubs) > 0 ) mfiles    <- mfiles[-mpftsubs]
    if(length(dpftsubs) > 0 ) dfiles    <- dfiles[-dpftsubs]
    print(mfiles,quote=F)
   }

  #mfiles    <- c('monthly_nep.dat','monthly_srp.dat')
  #mpftfiles <- NULL
 
  print(afiles,quote=F)
  print(paste('PFT=',pft,sep=''))
  if(annual)  print(afiles,quote=F)
  if(monthly) print(mfiles,quote=F)
  if(daily)   print(dfiles,quote=F)
  if(monthly&pft) print(mpftfiles,quote=F)
  if(daily&pft)   print(dpftfiles,quote=F)
 

  print("starting stitch") 
  # output variable loops 
  ########################
  data <- numeric(0)
  if(mc){
    if(annual)      data[] <- mclapply(afiles,stitch_annual,grids,wd,...)
    #if(monthly)     data[] <- mclapply(mfiles,stitch_subannual,grids,wd,atr=12,ad=2,...)
    if(daily)       data[] <- mclapply(dfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)    
    if(monthly&pft) data[] <- mclapply(mpftfiles,stitch_subannual,grids,wd,atr=12,ad=1,...)
    if(daily&pft)   data[] <- mclapply(dpftfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)    
    print("finished MC stitch") 
  } else {
    if(annual)      lapply(afiles,stitch_annual,grids,wd)
    #if(monthly)     lapply(mfiles,stitch_subannual,grids,wd,atr=12,ad=2,...)
    if(daily)       lapply(dfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)        
    if(monthly&pft) lapply(mpftfiles,stitch_subannual,grids,wd,atr=12,ad=1,...)
    if(daily&pft)   lapply(dpftfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)        
    print("finished non-MC stitch") 
  }
}



# Ancilliary functions
#######################

stitch_annual <- function(ifile,grids,wd,...){
 
  print(paste("starting stitch_annual: ",ifile,sep=""))
  ldf <- lapply(1:grids,read_grid,ifile,wd)
  #print(head(ldf[[9]])[,1:10])
  #df  <- rbind.fill(ldf)
  df  <- do.call('rbind',ldf)
  
  # write complete dataset to output dir
  setwd(wd)
  print(paste("starting write_sdgvm (stitch): ",ifile,sep=""))
  write_sdgvm(df,ifile)
}  

stitch_subannual <- function(ifile,grids,wd,...){
  
  print(paste("starting stitch_subannual: ",ifile,sep=""))
  ldv <- lapply(1:grids,scan_grid,ifile,wd)
  #ldv <- lapply(c(1:8,10:grids),scan_grid,ifile,wd)
  dv  <- unlist(ldv)         #PCM
  print(paste("starting process_sdgvm_matrix (stitch): ",ifile,sep=""))
  fdf <- process_sdgvm_matrix(dv,...) #PCM
  #dv  <- do.call('rbind',ldv) #PCM
  
  # write complete dataset to output dir
  setwd(wd)
  print(paste("starting write_sdgvm (stitch): ",ifile,sep=""))
  write_sdgvm(fdf,ifile) #pcm
}

process_sdgvm_matrix <- function(dv,atr,ad,nyears,...){
  
  slice <- function(i,mat,l){
    s <- (i-1)*l + 1
    e <- i*l
    mat[,s:e]
  }
  #pull lat and lon from the vector dv
  blk_l    <- nyears*(atr+ad)+2
  lat_subs <- seq(1,length(dv),blk_l)
  lon_subs <- seq(2,length(dv),blk_l)
  lat      <- dv[lat_subs]
  lon      <- dv[lon_subs]
  nsites   <- length(lat)
  latlon_subs <- c(lat_subs,lon_subs)
  dv_noll     <- dv[-latlon_subs] 

  # create a matrix of the data and re-arrange
  dmat        <- matrix(dv_noll,nrow=(atr+ad))
  if(atr==12) dmat <- dmat[1:13,]
  yrs         <- dmat[1,1:nyears]
  dmat        <- dmat[2:(atr+1),]  
  print(paste("C4c0: nyears,nsites,atr,ad ",nyears,nsites,atr,ad,sep=',')) 
  #print(paste("C4c:",dmat,sep=''))  
  dmatlist    <- lapply(1:nsites,slice,mat=dmat,l=nyears)
  dmatstack   <- do.call(rbind,dmatlist)
    
  #create dataframe for each site
  cbind( rep(lat,each=atr),
         rep(lon,each=atr),
         rep(1:atr,nsites),
         dmatstack )
}

read_grid <- function(g,ifile,wd){
  print(paste('read_grid:',ifile,'grid:',g))
  print(paste('read_grid:',wd,'grid',g,'/',sep=''))
  print(paste('read_grid:','grid',g, 'ifile= ',ifile))
  setwd(paste(wd,'grid',g,'/',sep=''))
  read.table(ifile,na.strings=c('*******','********','*********','**********','***********','************','Infinity'))
}

scan_grid <- function(g,ifile,wd){
  print(paste('scan_grid:','grid',g))
  print(paste('scan_grid:',wd,'grid',g,'/',sep=''))
  print(paste('scan_grid:','grid',g, 'ifile= ',ifile))
  setwd(paste(wd,'grid',g,'/',sep=''))
  scan(ifile,na.strings=c('*******','********','*********','**********','***********','************','Infinity'))
}

write_sdgvm <- function(df,file,w=10){
  print(paste("writing",file,sep=" "))
  write.table(format(df,width=w),file,row.names=F,col.names=F,quote=F) #PCM
  #write.fwf(df,file,width=w,justify='left',rownames=F,colnames=F,na='NA') #PCM
}



### END ###
