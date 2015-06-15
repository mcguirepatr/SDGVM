####################################
#
# global plotting function
#
# AWalker
# Jun 2015
#
###################################

#library(maps)
#library(fBasics)

.libPaths('~/bin/Rlibs')
library(lattice)
library(latticeExtra)
library(rworldmap)
library(rworldxtra)

  
  
###############################
#open annual data SDGVM files

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

write_title <- function(mod,vn,vs) as.expression(substitute(list(mod*' '*vname[vsub]),list(mod=mod,vname=vn,vsub=vs)))
write_var   <- function(vn,vs) as.expression(substitute(list(vname[vsub]),list(vname=vn,vsub=vs)))



###############################

region <- function(x,region){
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



### Lattice map plotting functions
# rworldmap mods
###############################

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



axis.ticks <- function (...,y=F,ticks = seq(-180,180,20), ticks2 = seq(-180,180,5)){
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

