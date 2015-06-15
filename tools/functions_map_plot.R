####################################
#
# global plotting function
#
# AWalker
# Feb 2012
#
###################################

library(maps)
library(fBasics)

###############################
#open IMOGEN-SDGVM / annual data SDGVM files
open <- function(var,mnames,wdpath){
  ifile <- paste(var,'.dat',sep='')
  print(ifile)
  x     <- read.table(paste(wdpath,ifile,sep=''),header=F,na.strings='********')
  names(x)<-mnames  
  x
}


###############################
grid_calc <- function(x,gsum,index,res){        
  
  # res 1 - longitude resolution
  # res 2 - latitude resolution
  
  #add lat and lon boundaries
  x$w_bound <- x$lon - res[1]/2
  x$e_bound <- x$lon + res[1]/2
  x$s_bound <- x$lat - res[2]/2
  x$n_bound <- x$lat + res[2]/2
  
  #add polygons to data
  poly_p <- cbind(x$w_bound, x$w_bound, x$e_bound, x$e_bound, x$s_bound, x$n_bound, x$n_bound, x$s_bound)
  x      <- cbind(x,poly_p)
  
  #parameters for summing global values 
  circ  <- 40008.0
  rad   <- circ/(2.0*pi)
  ydist <- circ*res[2]/360.0
  
  #take the global sum/mean values i.e. convert m-2 values to absolutes
  x$xdist<-2.0*pi*rad*cos(x$lat*pi/180.0)*res[1]/360.0

  x$area<-x$xdist*ydist
  
  ifelse(gsum, 
         x$sum<-x[,index]*x$area,
         x$sum<-(x[,index]*x$area)/sum(x$area)
  )
  
  #return processed dataframe
  x
  
}

###############################
scalef <- function(x,scalemin,scalemax,ncolours){
  #calculate subscript for colour array 
  
  range <- scalemax-scalemin
  xcorr <- x-scalemin+0.0001
  i     <- xcorr*ncolours/range
  i     <- ceiling(i)+1
  adj_ncolours <- ncolours+2
  
  #adjust for out of bounds values
  ifelse(is.na(x)|x<=0,i<-1,
         ifelse(i>ncolours+1,i<-adj_ncolours,ifelse(i<1,
                                                    i<-1,i<-i)))
  
  #returns
  i  
}

###############################
write_title  <- function(mod,vn,vs) as.expression(substitute(list(mod*' '*vname[vsub]),list(mod=mod,vname=vn,vsub=vs)))

###############################
write_var   <- function(vn,vs) as.expression(substitute(list(vname[vsub]),list(vname=vn,vsub=vs)))

###############################
###scale plotting function
# colours<-colv[[colo]]
# scalemin<-scaleminv[v]
# scalemax<-scalemaxv[v]

plot_scale<-function(scalemin,scalemax,colours,r,sfactor,sunit,vname,vsub){
  #plots a horizontal scale at the bottom of the figure 
  
  cex <- 0.7
  
  #add out of range high colour
  colours  <- c(colours,'red')
  ncolours <- length(colours)
  
  range0   <- scalemax-scalemin
  
  scalemax0<- scalemax
  scalemax <- scalemax+range0/(ncolours-1)
  range    <- scalemax-scalemin
  
  par(mar=c(0,7,0,7))
  plot(c(scalemin,scalemax),c(-2,4),type="n", ann=FALSE, axes=FALSE)
  
  y1<-c(0,1,1,0)
  y1<-y1+0.2
  xs<-scalemin
  xi<-range/ncolours
  xa<-xs
  xb<-xa+xi
  
  for(c in 1:length(colours)){
    x1<-c(xa,xa,xb,xb)
    polygon(x1,y1,col=colours[c],border='black')
    xa<-xb
    xb<-xa+xi
  }
  
  text(scalemin,-1,as.character(round(scalemin,r)),cex=cex)
  text(scalemax0,-1,paste('>',as.character(round(scalemax0,r)),sep=''),cex=cex)
   
  nticks<-(scalemax0-scalemin)/sfactor-1
  
  if(scalemin>=0){
    for(si in 1:nticks){
      x<-xi*si+scalemin
      text(x,-1,as.character(round(x,r)),cex=cex)
    }}
  
  text(scalemin,2.4,write_var(vname,vsub),pos=4,cex=cex*1.5)
  text((scalemin+range*0.2),2.4,parse(text=sunit),pos=4,cex=cex*1.5)
  
  #end function
}   

###############################
region<-function(x,region){
  #subset dataset 'x' by region 
  
  # returns 'xylims' a 6 element vector composing of:
  # (W lon boundary, E lon boundary, S lat boundary, N lat boundary, lon tick res, lat tick res)
  
  if(region=='global') {
    xylims <- c(-180,180,-60,90,90,30)
  }
  
  else if(region=='europe'){
    xylims <- c(-10,80,30,70,15,10)
  }
  
  else if(region=='tropics'){
    xylims<-c(-90,180,-30,30,90,15)
  }
  
  else if(region=='boreal NA'){
    xylims <- c(-180,-45,45,70,15,10)
  }
  
  else if(region=='siberia'){
    xylims <- c(80,180,45,70,20,15)
  }
  
  else if(region=='north'){
    xylims <- c(-180,180,-45,70,90,15)    
  }
  
  else if(region=='amazonia'){
    xylims <- c(-90,-30,-30,15,30,15)
  }
  
  #returns
  xylims
}


###############################   
#integrate per m values to the total amount for the globe/region
# x <- mydata[,1:2]
# y <- mydata[,3:length(mydata)]
# res<-mod_res
# gm<-gmean[v]
area_integrate <- function(x,y,res,gm){
  # x is a dataframe with variables 'lat' , 'lon'
  # y is a dataframe of the actual data, n columns the same as the columns of the output  
  
  circ  <- 40008.0
  rad   <- circ/(2.0*pi)
  ydist <- circ*res[2]/360.0
  
  x$xdist <- 2.0*pi*rad*cos(x$lat*pi/180.0)*res[1]/360.0
  x$area  <- x$xdist*ydist
  areai   <- y*x$area 
  if(gm) areai       <- areai/sum(x$area)
  global_sum         <- apply(as.matrix(areai),2,sum,na.rm=T)
  if(!gm) global_sum <- global_sum*1e-9
  global_sum
}



###############################   
#map plotting function
# v<-1
# variable<-vars[v]
# scalemin<-scaleminv[v]
# scalemax<-scalemaxv[v]
# sfactor<-sfactor[v]
# vname<-vname[v]
# varsub<-vsub[v]
# sunit<-sunit[v]
# colo<-colo[v]
# gsum<-gsum[v]
# unit<-unit[v]
# colv<-colv
# x<-subdata
plot_map <- function(x,index,scalemin,scalemax,sfactor,vname,varsub,sunit,colo,gsum,unit,colv,pregion='global',model,year,mar){
#rl is the last column in the data frame before the ploygon columns
#could add polygon calculation to this function 
#index is the column with the data to be plotted

  #mar=c(0.5,2,0.5,0)
    
  rl      <- length(x) - 11
  print(rl)
  cl      <- length(x[,1])  
  ncol    <- length(colv[[colo]])
  colours <- c('white',colv[[colo]],'red')
  
  ifelse(gsum, 
         global_sum<-round(sum(x$sum,na.rm=T)*1e-9,1),
         global_sum<-round(sum(x$sum,na.rm=T),2)
  )
  
  #write lat lon -90-90 & -180-180
  date_line<-min(x$e_bound[which(abs(x$e_bound-180)==min(abs(x$e_bound-180)))])
  x$e_bound[x$e_bound>date_line]  <- x$e_bound[x$e_bound>date_line]-360
  x$w_bound[x$w_bound>=date_line] <- x$w_bound[x$w_bound>=date_line]-360

  #subset dataset 
  #remove nas
  x <- x[is.finite(x[,index]),]  
  #by region
  lims  <- region(x,pregion)
  #subset data to region
  x <- subset(x,w_bound>lims[1])
  x <- subset(x,e_bound<lims[2])
  x <- subset(x,s_bound>lims[3])
  x <- subset(x,n_bound<lims[4])
#   print(head(x))
#   lims <- c(-90,-30,-30,15)
#   lims <- c(-90,180,-30,30)
    
  #plot
#   plot(getMap('high'),lwd=0.5,xlim=c(-180,180),ylim=c(-60,90))
  map("world", xlim=c(lims[1],lims[2]), ylim=c(lims[3],lims[4]), interior=T, bg="white", col="black",cex=1,mar=mar)
  axis(2,at=seq(lims[3],lims[4],lims[6]),pos=lims[1],las=TRUE) 
  axis(1,at=seq(lims[1],lims[2],lims[5]),pos=lims[3]) 
  axis(3,labels=F,tck=0,pos=lims[4],at=c(lims[1],lims[2]))
  axis(4,labels=F,tck=0,pos=lims[2],at=c(lims[3],lims[4]))

  title(main=write_title(paste(model,year),vname,varsub))
  
  print(rl);print(length(x));print(cl)
  for (i in (1:cl)){
    polygon(x[i,(rl+1):(rl+4)], x[i,(rl+5):(rl+8)], col=colours[scalef(x[i,index],scalemin,scalemax,ncol)],border=NA)
  }
  
#   plot(getMap('high'),add=T,lwd=0.5,xlim=c(-180,180),ylim=c(-60,90))
  map("world",  xlim=c(lims[1],lims[2]), ylim=c(lims[3],lims[4]), interior=T, bg="white", col="black",add=T,mar=mar)
  text(x=-135,y=-45,labels=ifelse(gsum|gmean,paste(global_sum,unit),''))
  
  abline(h=0,lty=2,lwd=0.6)
  abline(h=c(-23.3,23.3,66.5),lty=3,lwd=0.6)
}   
   


### Lattice map plotting functions
#rworldmap mods
###############################
# dev.off()
# mapGriddedData(mapRegion='Africa')
# mapGriddedData
# rmw <- rworldmap:::rwmNewMapPlot
# rworldmap:::setMapExtents
 
library(lattice)
library(latticeExtra)
library(rworldmap)
library(rworldxtra)

plotmap <- function(i,map) {
  pm <- function(j,i,pg){
    panel.polygon(x=pg[[j]]@coords[,1],y=pg[[j]]@coords[,2],border='black',lwd=0.5)                
  }
  polys <- map@polygons[[i]]@Polygons
  lapply(1:length(polys),pm,i=i,pg=polys)
}

plot_map_lattice <- function(x,index,lab,pregion='global',model,
                             year,mar,gs=NULL,...){
  #initialise
  world <- getMap()    
  
  ncol    <- length(lab$cols)
  colours <- lab$cols
  
  if(!is.null(gs)){
    if(lab$gsum) {
           global_sum <- round(gs[,index-2],1)
    } else global_sum <- round(gs[,index-2],2)
  } else global_sum <- gs
  
  #subset dataset 
  #remove nas
  x <- x[is.finite(x$plotdata),]  
  #by region
  lims  <- region(x,pregion)
  #subset data to region
  x <- subset(x,lon>lims[1]-2)
  x <- subset(x,lon<lims[2]+2)
  x <- subset(x,lat>lims[3]-1.25)
  x <- subset(x,lat<lims[4]+1.25)
  
  #   print(head(x))
  
  #plot
  levelplot(plotdata~lon*lat|as.factor(year)*run,x,...,
            as.table=T,
            main=write_var(lab$name,lab$nsub),
            # scale bar
            col.regions=colours,cuts=ncol-2,
            at=lab$at,
            colorkey=list(labels=list(labels=lab$at,at=lab$at)),
            # axes
            aspect='iso',xlim=c(lims[1:2]),ylim=c(lims[3:4]),xlab=NULL,ylab=NULL,
            scales=list(alternating=F),
            xscale.components=function(...,top=F,lim=lims[1:2]){
              axis.ticks(lim=lim,...)},
            yscale.components=function(...,lim=lims[3:4]){
              axis.ticks(lim=lim,y=T,...)},
            # panel function
            panel=function(...){
              panel.levelplot(...)
              lapply(1:length(world@polygons),plotmap,map=world)
              panel.abline(h=c(0,-23.3,23.3,66.5,-66.5),lty=c(2,3,3,3,3,3),lwd=0.6)
              panel.text(x=-160,y=-35,pos=4,labels=ifelse(lab$gsum|lab$gmean,t(global_sum)[panel.number()],''))      
              panel.text(x=-160,y=-45,pos=4,labels=ifelse(lab$gsum|lab$gmean,parse(text=lab$sunit),''))      
            }) #+ layer_(panel.2dsmoother(..., n = 200))
}   

axis.ticks <- function (...,y=F,ticks = seq(-180,180,20), ticks2 = seq(-180,180,5)) 
{
  ans    <- xscale.components.default(...)
  ticks2 <- ticks2[!(ticks2 %in% ticks)]
  ans$bottom$ticks$at  <- c(ticks, ticks2)
  ans$bottom$ticks$tck <- c(rep(-1,length(ticks)), rep(-0.5,length(ticks2)))
  ans$bottom$labels$at <- ticks
  ans$bottom$labels$labels <- ticks
  if(y){
    ans    <- yscale.components.default(...)
    ticks  <- c(seq(-90,90,20))
    ticks2 <- seq(-90,90,5)
    ticks2 <- ticks2[!(ticks2 %in% ticks)]
    ans$left$ticks$at  <- c(ticks, ticks2)
    ans$left$ticks$tck <- c(rep(-1,length(ticks)), rep(-0.5,length(ticks2)))
    ans$left$labels$at <- ticks
    ans$left$labels$labels <- ticks
  }
  ans
}

# xyplot(1:10~1:10,
#        panel=function(...){
#          panel.xyplot(...)
#          panel.text(x=5,y=5,'text')
#          panel.text(x=8,y=5,parse(text='mu*mol'))
#        })


