##############################
#
# Open a CRUNCEP netcdf files 
#
# AWalker
# Oct 2013
#
##############################

rm(list=ls())

### HYDE LUH2 v2h files genetrated by Gerge Hurtt's group, see README for details
# describing annual land cover from 850 to 2017

### ESA CCI LCP 2014 files processed by Natasha MacBean for SDGVM PFTs, see README for details
# 1   100 BARE
# 2   100 Ev_Bl
# 3   100 Dc_Bl
# 4   100 Ev_Nl
# 5   100 Dc_Nl
# 6   100 Shrub
# 7   100 C3 
# 8   100 C4
# 9   100 C3crop
# 10  100 C4crop



### Functions
#######################################

library(gdata)

join_hyde <- function(v) {
  # takes a 15 element vector, 1:10 - SDGVM PFTs from ESA data, 11:15 HYDE aggregated land cover 
  # and combines into a 10 element vector of SDGVM PFTs
  #
  # SDGVM PFTs: 1 BARE, 2 Ev_Bl, 3 Dc_Bl, 4 Ev_Nl, 5 Dc_Nl, 6 Shrub, 7 C3, 8 C4, 9 C3crop, 10 C4crop
  # HYDE aggregated land cover: 11 forest, 12 non-forest, 13 C3 crop, 14 C4 crop, 15 pasture   
  
  # first forest area and potential forest area from HYDE is used to adjust forest area in ESA
  # then C3 and C4 croplands in HYDE are used to adjust ESA
  # C3 and C4 grassland area is then adjusted to account for changes in pasture and natural grassland cover
  # rangelands in hyde and shrubs in esa are not explicitly included in the calculation of change 
  # though a proportion of them could be accounted for implicitly if rangelands habour some forest vegetation
  
  # print(v)
  
  if(is.finite(sum(v[1:10])) & !is.finite(sum(v[11:15])) ) {
    # when ESA has landcover but hyde not
    # arises from different land masks
    return(v[1:10])
    
  } else if(is.finite(sum(v))) {
    
    # create output vector 
    # print(v)
    ov <- v[1:10]
    
    # subscripts for ESA forest PFTs
    fid <- 2:6
    
    # stores for increase in ESA but when either no forest or grass cover in ESA 
    nofcov  <- 0
    nogcov  <- 0
    
    # change forest cover
    # ESA forest cover
    for_esa <- sum(v[fid])
    # If forest cover is greater in hyde
    if(v[11]>for_esa) {
      diff   <- v[11] - for_esa
      if(for_esa>0) {
        # increase forest cover proportionally
        for(f in fid) ov[f]  <- ov[f] + diff * (v[f] / for_esa)
      } else {
        # no forest cover in ESA - forest cover stored for post-processing
        # write(paste('no ESA forest','BARE:',ov[1],'C3 CROP:',ov[9],'C4 CROP:',ov[10]),'1nofor_error.txt',append=T)
        nofcov <- diff
      }
    } else if(for_esa>v[11]) {
      # forest cover greater in ESA
      diff   <- for_esa - v[11]
      # decrease forest cover proportionally
      for(f in fid) ov[f]  <- ov[f] - diff * (v[f] / for_esa)
    }
    
    
    # increase or decrease crop cover 
    # C3 crops
    if(v[9]!=v[13])  ov[9]  <- v[13]
    
    # C4 crops
    if(v[10]!=v[14]) ov[10] <- v[14]
    
    
    # increase or decrease grass (pasture) cover 
    # total cover so far in hyde ESA combined dataset - accounts for forest cover as yet unassigned to a PFT
    esahyde_cov <- sum(ov) + nofcov
    
    # total cover in hyde
    hyde_cov    <- sum(v[11:15])
    
    # ESA grasslands (potential pasture)
    past_esa    <- v[7] + v[8]
    
    # if new cover is different from hyde cover
    past_err <- F
    if(esahyde_cov!=hyde_cov) {
      diff   <- hyde_cov - esahyde_cov
      if(diff<0) {
        if(abs(diff)<past_esa) {
          
          # change C3/C4 grass cover proportionally
          ov[7]  <- ov[7] + diff * (v[7] / past_esa)
          ov[8]  <- ov[8] + diff * (v[8] / past_esa)
          
        } else {
          
          # ESA pasture needs reduction but reduction is more than pasture cover
          ov[7] <- 0
          ov[8] <- 0
          
          # assume remaining reduction from bare ground
          diff <- abs(diff) - past_esa
          ov[1] <- if(ov[1]>=diff) ov[1] - diff else 0
        }   
        
      } else if(diff>0) {
        if(past_esa>0) {
          
          # change C3/C4 grass cover proportionally
          ov[7]  <- ov[7] + diff * (v[7] / past_esa)
          ov[8]  <- ov[8] + diff * (v[8] / past_esa)
          
        } else {
          
          # no grass cover in ESA - grass cover stored for post-processing
          # write(paste('no ESA grass','BARE:',ov[1],'C3 CROP:',ov[9],'C4 CROP:',ov[10]),'1nopast_error.txt',append=T)
          nogcov <- diff
        }
      }
    }
    
    # test for cover sum >100%
    if( sum(ov[1:10])-100 > 1e-3 ) write(paste('Error,cov:',sum(ov[1:10]),hyde_cov,'. Bare:',v[1]),'1cov_error.txt',append=T)
    
    # assign as yet unassigned forest cover as negative value to DcBl
    if(nofcov>0) ov[3]  <- -nofcov
    # assign as yet unassigned grass cover as negative value to C3 grass
    if(nogcov>0) ov[7]  <- -nogcov
    
    return(ov[1:10])
    
  } else {
    return(rep(NA,10))
  }
}


# Assigns forest or grass pfts to ESA HYDE according to latitude when no forest or grass PFTs exist in ESA but do in HYDE 
f_lat_assignPFT <- function(j,m,loff) {
  # this function assigns forest and grassland cover to an appropriate PFT in the combined HYDE ESA dataset
  # when there was no forest or grass cover in the ESA data.
  # in this case the forest cover is assigned as a negative value to deciduous braodleaved PFT m[3,] 
  # in this case the grass  cover is assigned as a negative value to C3 grass PFT m[7,]
  # if in temperate latitudes this negative value is simply switched to positive
  # if in tropical latitudes this negative value is switched to positive and assigned to evergreen broadleaved PFT or C4 grass 
  
  # m is a pft x lat matrix
  # loff is the change in subscript by lat - 0 for temperate (DcBL & C3), 1 for tropical (EvBl & C4)
  
  # if there was no forest or grass cover in the ESA data 
  if(any(m[,j]<0)) {
    # switch DcBl and EvBl PFT to allow loff to work for both forest and grass
    m[2:3,j] <- m[3:2,j]
    # subscrips of negative cover
    sub      <- which(m[,j]<0)
    # take absolute value of negative covers and assign to appropriate PFT
    m[sub+loff[j],j] <- abs(m[sub,j])
    # if tropical PFT is appropriate zero negative cover in temperate PFT
    if(loff[j]!=0) m[sub,j] <- 0
    # switch DcBl and EvBl PFT back to original placement in vector
    m[2:3,j] <- m[3:2,j]
    # return vector
    m[,j]
  } else {
    m[,j]
  }
}


# write functions - intended to be called from lapply
write_landcover <- function(i,mat) {
  write.fwf(matrix(round(mat[i,],0),ncol=1),width=3,nsmall=0,paste(paste(fname_s,i,time[t],sep='-'),send,sep=''),colnames=F,sep='')
}

write_landcover_array <- function(i,arr) {
  v <- as.vector(arr[,,i])
  write.fwf(matrix(round(v,0),ncol=1),width=3,nsmall=0,paste(paste(fname_s,i,time[t],sep='-'),send,sep=''),colnames=F,sep='')
}



### Initialise
#######################################

library(lattice)
library(grid)
library(parallel)
library(ncdf4)
library(viridis)

#srcdir <- '/mnt/disk2/script_library/R/map_tools/'
srcdir  <- '/gws/nopw/j04/nexcs/pmcguire/TRENDYv8/scripts/'
source(paste(srcdir,'regrid.R',sep=''))

# paths
date <- Sys.Date()
#dir  <- '/mnt/disk3/'
# for speed: cp -pr  /gws/nopw/j04/nexcs/pmcguire/sdgvmD/data/land_use/global/ESACCILCP2014 /work/scratch-pw2/pmcguire/TRENDYv12/db/landcover/ESACCILCP2014
dir  <- '/work/scratch-pw2/pmcguire/TRENDYv12/'
# wdh  <- paste(dir,'/Databases/landcover/hyde/3.2/hyde32_baseline_for_gcp2016/',sep='')
# wdh  <- paste(dir,'/Databases/landcover/hyde/LUH2/v2h/',sep='')
# wdh  <- paste(dir,'/Databases/landcover/hyde/LUH2/v2h.1/',sep='')
# wdh  <- paste(dir,'/Databases/landcover/hyde/LUH2/v2h_2020/',sep='')
# wdh  <- paste(dir,'/Databases/landcover/hyde/LUH2/v2h_2022/',sep='')
# wdg  <- paste(dir,'/Databases/landcover/ESACCILCP2014/30min',sep='')
wdh  <- paste(dir,'db/landcover/hyde/LUH2/v2h_2023/',sep='')
wdg  <- paste(dir,'db/landcover/ESACCILCP2014/30min',sep='')
# wdo  <- paste(dir,'/Databases/SDGVMdata/land_use/global/ESACCILCP2014_LUH2v2h/30min',sep='')
# wdo  <- paste(dir,'/Databases/SDGVMdata/land_use/global/ESACCILCP2014_LUH2v2h.1/30min',sep='')
# wdo  <- paste(dir,'/Databases/SDGVMdata/land_use/global/ESACCILCP2014_LUH2v2h_2022/30min',sep='')
wdo  <- paste(dir,'db/SDGVMdata/land_use/global/ESACCILCP2014_LUH2v2h_2023/30min',sep='')

# filename variables
fname    <- 'states'
fname_s  <- 'cont_lu'
esadate  <- 2009
time     <- 850:2022
fend     <- '.nc'
send     <- '.dat'
lon_res  <- 1440
lat_res  <- 720
lon_ress <- 720
lat_ress <- 360

# time variables  
syear    <- 1700
# syear    <- 1894

# LUH2 variables
ncvars   <- c('primf','secdf','secdn','primn',
              'c3ann','c3per','c3nfx','c4ann','c4per',
              'pastr','range','urban')
misval   <- 1e19

# processing parameter
# latitudinal offset for assigning correct forest or grass PFT when none exist in ESA
lat_offset <- c(rep(0,132),rep(1,96),rep(0,132))

# ESA PFTs
esa_pfts <- c('BARE','Ev_Bl','Dc_Bl','Ev_Nl','Dc_Nl','Shrub','C3','C4','C3crop','C4crop')

# plot parameters
col <- c('grey85',rev(viridis(50)))
hyde_agg_names <- c('forest','nat non-forest','c3 crop','c4 crop','pasture','urban')



### Program proper
####################################

# setup ESA arrays
esaarray <- array(0,dim=c(lon_ress,lat_ress,15)) 

# Open ESA CCILCP 2014 dataset 
setwd(wdg)
for(i in 1:length(esa_pfts)) {
  esav          <- scan( paste(paste(fname_s,i,esadate,sep='-'),send,sep='') )
  # esamat <- if(i==1) esav else cbind(esamat,esav)
  esaarray[,,i] <- as.matrix(esav,nrow=lon_ress)
}
rm(esav)

# convert missing values      
esaarray[which(esaarray==255)] <- NA      

# # plot ESA
# setwd('../plots/')
# p1 <- levelplot(esaarray[,,1:10],xlim=c(1,720),ylim=c(360,1),col.regions=col,contour=F,at=seq(0,100,2),strip=strip.custom(factor.levels=esa_pfts),as.table=T)
# pdf(paste('ESA.pdf',sep='')); print(p1); dev.off()


## loop through Hyde data years (that are pertinent to the project, starting 1700)
tia <- which(time==syear):length(time)
time[tia]

# year loop
t <- which(time==syear)
for(t in tia) {
  
  # Open HYDE LUH2 dataset
  setwd(wdh)
  mynetcdf <- nc_open(paste(fname,fend,sep=''))

  # setup arrays
  dummy1   <- array(0,dim=c(lon_res,lat_res)) 
  dummya   <- array(0,dim=c(lon_ress,lat_ress)) 
  dummya1  <- array(0,dim=c(lon_ress,lat_ress,length(ncvars)+1)) 
  dummya2  <- array(0,dim=c(lon_ress,lat_ress,length(hyde_agg_names)+1)) 

  # extract data
  setwd('./plots/')
  for(v in 1:length(ncvars)) {
    
    # read hyde variable
    dummy1 <- ncvar_get(mynetcdf,ncvars[v],start=c(1,1,t),count=c(lon_res,lat_res,1))
    
    # regrid from 0.25 to 0.5
    dummya <- reduce_res(dummy1)

    # convert missing values      
    dummya[which(dummya > misval)] <- NA      

    # all hyde default classes array 
    dummya1[,,v] <- dummya
    
    # aggregate Hyde landcover types 
    # primary and secondary forest
    if(v<=2)       dummya2[,,1] <- dummya2[,,1] + dummya
    # primary and secondary non-forest
    else if(v<=4)  dummya2[,,2] <- dummya2[,,2] + dummya
    # c3 crops
    else if(v<=7)  dummya2[,,3] <- dummya2[,,3] + dummya
    # c4 crops
    else if(v<=9)  dummya2[,,4] <- dummya2[,,4] + dummya
    # # pasture and rangelands
    # else if(v<=11) dummya2[,,5] <- dummya2[,,5] + dummya
    # pasture and not rangelands
    else if(v==10) dummya2[,,5] <- dummya2[,,5] + dummya
    # rangelands - add to non-forest
    else if(v==11) dummya2[,,2] <- dummya2[,,2] + dummya
    # urban
    else if(v==12) dummya2[,,6] <- dummya
  } # data extraction loop
  
  # Close HYDE LUH2 dataset
  rm(mynetcdf)
  
  # plot hyde default classes
  dummya1[,,length(ncvars)+1] <- apply(dummya1,1:2,sum)
  p1 <- levelplot(dummya1,xlim=c(1,720),ylim=c(360,1),
                  col.regions=col,contour=F,at=seq(0,1.02,0.02),as.table=T,
                  #panel=panel.levelplot.raster)#,
                  #useRaster=T,
                  strip=strip.custom(factor.levels=c(ncvars,'sum')))
  pdf(paste('Hyde_',time[t],'.pdf',sep='')); print(p1); dev.off()
  
  # change to percent      
  dummya2 <- dummya2 * 100      
  
  # plot hyde aggregate classes
  dummya2[,,length(hyde_agg_names)+1] <- apply(dummya2,1:2,sum)
  p1 <- levelplot(dummya2,xlim=c(1,720),ylim=c(360,1),
                  col.regions=col,contour=F,at=seq(0,102,2),as.table=T,
                  #useRaster=T,
                  strip=strip.custom(factor.levels=c(hyde_agg_names,'sum')))
  pdf(paste('Hyde_agg_',time[t],'.pdf',sep='')); print(p1); dev.off()
  
  # create processing matrix
  esaarray[,,11:15] <- dummya2[,,1:5]
  
  # process, returns an array of PFT, lon, lat  
  out   <- apply(esaarray,c(1,2),join_hyde)
  
  # convert missing values to SDGVM missing value
  out[is.na(out)] <- 255

  # deal with forest or grass cover where no forest or grass cover existed in ESA
  # permute out array to PFT, lat, lon
  out <- aperm(out,c(1,3,2))
  out <- 
    sapply(1:dim(out)[3],
           function(i,a) sapply(1:dim(a)[2],f_lat_assignPFT,m=a[,,i],loff=lat_offset) ,
           a=out,simplify='array')
  # permute to lon, lat, PFT
  out <- aperm(out,c(3,2,1))
  
  # write PFT year data
  setwd(wdo)
  lapply(1:10,write_landcover_array,arr=out)

  # print combined ESA HYDE
  setwd('../plots/')
  p1 <- levelplot(out,xlim=c(1,720),ylim=c(360,1),
                  col.regions=col,contour=F,at=seq(0,102,2),
                  strip=strip.custom(factor.levels=esa_pfts),
                  #useRaster=T,
                  as.table=T)
  # p2 <- levelplot(out,xlim=c(1,720),ylim=c(360,1),col.regions=c('red','grey40','grey80'),contour=F,at=c(-100,-1,0,100),strip=strip.custom(factor.levels=esa_pfts),as.table=T)
  # pdf(paste('HydeESA_',time[t],'.pdf',sep='')); print(p1); print(p2); dev.off()
  pdf(paste('HydeESA_',time[t],'.pdf',sep='')); print(p1); dev.off()
  
} ## end loop



### END ###
