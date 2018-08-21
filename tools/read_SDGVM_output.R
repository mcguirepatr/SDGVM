######################
#
# Read SDGVM output
#
# AWalker
# Oct 2013
#
######################

.libPaths('~/bin/Rlibs')
library(plyr)
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
  files  <- files[-which(files=='site_info.dat')]
  files  <- files[-which(files=='simulation.dat')]
  files  <- files[-which(files=='diag.dat')]
  isubs  <- grep('init',files)
  files  <- files[-isubs]
  print(files)
  
  msubs  <- grep('monthly',files)
  dsubs  <- grep('daily',files)
  afiles <- files[-c(dsubs,msubs)]
  mfiles <- files[msubs]
  dfiles <- files[dsubs]
 
  mpftsubs  <- grep('[A-Z]',mfiles)
  dpftsubs  <- grep('[A-Z]',dfiles)
  mpftfiles <- mfiles[mpftsubs]
  dpftfiles <- dfiles[dpftsubs]
  mfiles    <- mfiles[-mpftsubs]
  dfiles    <- dfiles[-dpftsubs]
 
  print(afiles,quote=F)
  if(annual)  print(afiles,quote=F)
  if(monthly) print(mfiles,quote=F)
  if(daily)   print(dfiles,quote=F)
  if(monthly&pft) print(mpftfiles,quote=F)
  if(daily&pft)   print(dpftfiles,quote=F)
  
  # output variable loops 
  ########################
  data <- numeric(0)
  if(mc){
    if(annual)      data[] <- mclapply(afiles,stitch_annual,grids,wd,...)
    if(monthly)     data[] <- mclapply(mfiles,stitch_subannual,grids,wd,atr=12,ad=2,...)
    if(daily)       data[] <- mclapply(dfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)    
    if(monthly&pft) data[] <- mclapply(mpftfiles,stitch_subannual,grids,wd,atr=12,ad=1,...)
    if(daily&pft)   data[] <- mclapply(dpftfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)    
  } else {
    if(annual)      lapply(afiles,stitch_annual,grids,wd)
    if(monthly)     lapply(mfiles,stitch_subannual,grids,wd,atr=12,ad=2,...)
    if(daily)       lapply(dfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)        
    if(monthly&pft) lapply(mpftfiles,stitch_subannual,grids,wd,atr=12,ad=1,...)
    if(daily&pft)   lapply(dpftfiles,stitch_subannual,grids,wd,atr=360,ad=1,...)        
  }
}



# Ancilliary functions
#######################

stitch_annual <- function(ifile,grids,wd,...){
  
  print(ifile)
  ldf <- lapply(1:grids,read_grid,ifile,wd)
  #print(head(ldf[[9]])[,1:10])
  #df  <- rbind.fill(ldf)
  df  <- do.call('rbind',ldf)
  
  # write complete dataset to output dir
  setwd(wd)
  write_sdgvm(df,ifile)
}  

stitch_subannual <- function(ifile,grids,wd,...){
  
  print(ifile)
  ldv <- lapply(1:grids,scan_grid,ifile,wd)
  #ldv <- lapply(c(1:8,10:grids),scan_grid,ifile,wd)
  dv  <- unlist(ldv)
  fdf <- process_sdgvm_matrix(dv,...)
  
  # write complete dataset to output dir
  setwd(wd)
  write_sdgvm(fdf,ifile)
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
  dmatlist    <- lapply(1:nsites,slice,mat=dmat,l=nyears)
  dmatstack   <- do.call(rbind,dmatlist)
    
  #create dataframe for each site
  cbind( rep(lat,each=atr),
         rep(lon,each=atr),
         rep(1:atr,nsites),
         dmatstack )
}

read_grid <- function(g,ifile,wd){
  #print(paste(ifile,'grid:',g))
  setwd(paste(wd,'grid',g,'/',sep=''))
  read.table(ifile,na.strings=c('*******','********','*********','**********','***********','************'))
}

scan_grid <- function(g,ifile,wd){
  #print(paste('grid',g))
  setwd(paste(wd,'grid',g,'/',sep=''))
  scan(ifile,na.strings=c('*******','********','*********','**********','***********','************'))
}

write_sdgvm <- function(df,file,w=10){
  #write.table(format(df,width=w),file,row.names=F,col.names=F,quote=F)
  write.fwf(df,file,width=w,justify='left',rownames=F,colnames=F,na='NA')
}

