####################################
#
# global plotting function
#
# AWalker
# Jun 2015
#
###################################

.libPaths('~/bin/Rlibs')
#library(lattice)
#library(latticeExtra)
#library(rworldmap)
#library(rworldxtra)

  
  
###############################
#open annual data SDGVM files

open <- function(var,mnames,wdpath){
  ifile    <- paste(var,'.dat',sep='')
  x        <- read.table(paste(wdpath,ifile,sep=''),header=F,na.strings='********')
  x        <- as.data.frame(apply(x,2,as.numeric))
  #print(head(x))
  names(x) <- mnames  
  x
}

  

###############################

write_title  <- function(mod,vn,vs) as.expression(substitute(list(mod*' '*vname[vsub]),list(mod=mod,vname=vn,vsub=vs)))
write_title1 <- function(vn,vs,txt) as.expression(substitute(list(vname[vsub]*' '*text),list(vname=vn,vsub=vs,text=txt)))
write_var    <- function(vn,vs)     as.expression(substitute(list(vname[vsub]),list(vname=vn,vsub=vs)))



###############################

region <- function(region){
  # returns 'xylims' an 8 element vector composing of:
  # (W lon boundary, E lon boundary, S lat boundary, N lat boundary, )
  
  if(region=='global') {
    xylims <- c(-180,180,-60,90,90,30,30,-43)
  }
  
  else if(region=='europe'){
    xylims <- c(-10,80,30,70,15,10,0,0)
  }
  
  else if(region=='tropics'){
    #xylims <- c(-120,180,-23.5,23.5,90,15,55,-15)
    xylims <- c(-120,180,-30,30,90,15,55,-15)
  }
  
  else if(region=='neotropics'){
    xylims <- c(-120,-30,-23.5,23.5,90,15,-110,-15)
  }
  
  else if(region=='afrotropics'){
    xylims <- c(-30,60,-23.5,23.5,90,15,-25,-15)
  }
  
  else if(region=='asiantropics'){
    xylims <- c(60,180,-23.5,23.5,90,15,65,-15)
  }
  
  else if(region=='boreal NA'){
    xylims <- c(-180,-45,45,70,15,10,0,0)
  }
  
  else if(region=='siberia'){
    xylims <- c(80,180,45,70,20,15,0,0)
  }
  
  else if(region=='north'){
    xylims <- c(-180,180,-45,70,90,15,0,0)    
  }
  
  else if(region=='amazonia'){
    xylims <- c(-90,-30,-23.5,15,30,15,0,0)
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
  # convert to mean value per global area, or convert to human readable units (e.g. for C in gm-2 converts to Pg)
  if(gm) areai/sum(x$area) else areai*1e-9
}



### Lattice map plotting functions
# rworldmap mods
###############################

plotmap <- function(i,map) {
  pm <- function(j,i,pg){
    panel.polygon(x=pg[[j]]@coords[,1],y=pg[[j]]@coords[,2],border='black',lwd=0.1)                
  }
  polys <- map@polygons[[i]]@Polygons
  lapply(1:length(polys),pm,i=i,pg=polys)
}



plot_map_lattice <- function(x,lab,pregion='global',gs=NULL,norm=norm,diff=diff,mask=mask,labcex=1,stripprint,...){
  # This function plots world (or subsetted) maps of gridded data

  # expects a dataframe, 'x', with the columns 'lat', 'lon', 'plotdata', 'year', and 'run'
  # 'run' is a factorial column with the label for the simulation that produced 'plotdata'  
  # 'year' is also a factorial column, I'm not sure yet that this is compatible with 'gs'
  # 'gs' is a vector or matrix of globally area integrated values for the data in 'x$plotdata'

  #initialise
  world <- getMap()    
  
  if(!is.null(gs)){
    if(lab$gsum) {
           global_sum <- round(gs,1)
    } else global_sum <- round(gs,2)
  } else global_sum <- gs
  
  #subset dataset 
  #remove nas
  x <- x[is.finite(x$plotdata),]  
  #by region
  lims  <- region(pregion)
  #subset data to region
  x <- subset(x,lon>lims[1]-2)
  x <- subset(x,lon<lims[2]+2)
  x <- subset(x,lat>lims[3]-1.25)
  x <- subset(x,lat<lims[4]+1.25)
  
  #   print(head(x))
   
  # modify plotting parameters if norm or diff is used normalise data 
  text <- NULL
  if(!is.null(norm)) {
    
    text <- paste('normalised:',norm)

    if(norm=='sd') {
      #lab$cols <- col.neg
      lab$at   <- c(seq(-4.0,-0.5,0.5),seq(0.5,4.0,0.5))
      lab$at   <- c(seq(-3,3,0.5)) 
      #lab$at   <- c(seq(-2,2,0.5)) 
      lab$unit <- 'sigma' 
      lab$cols <- viridis(length(lab$at)-1)
    } else if(any(x$plotdata<0)) {
      lab$cols <- col.neg
      lab$at   <- c(seq(-1.3,0.1,0.2),seq(0.1,1.3,0.2)) 
    } else {
      lab$cols <- col.inc.gpp
      lab$at   <- seq(0,1.3,0.1) 
    }
  }

  if(!is.null(diff)) {
    text      <- paste(text,'. Difference plot',sep='')
    lab$cols  <- rev(col.neg)
#    lab$at    <- lab$at / 2 
    if(is.null(norm))  lab$at <- c(seq(-600,-100,100),seq(100,600,100)) else lab$at <- lab$at * 0.5 
    lab$at    <- round(lab$at,2)  
  }

  lab$lab <- lab$at

  # expand range of labels based on data
  if(F){
    lab$at_orig <- lab$at
    lla         <- length(lab$at)
    if(round(min(x$plotdata),1)<lab$at[1])   lab$at[1]   <- round(min(x$plotdata),1)
    if(round(max(x$plotdata),1)>lab$at[lla]) lab$at[lla] <- round(max(x$plotdata),1)
    # rescale range based on expanded max/min
    if(F) lab$at <- seq(lab$at[1],lab$at[lla],(lab$at[lla]-lab$at[1])/(lla-1))
    lab$at      <- round(lab$at,2)  
  }
  if(T){
    lab$at_orig <- lab$at
    lla         <- length(lab$at)
    if(round(min(x$plotdata),1)<lab$at[1]) { 
      x$plotdata[x$plotdata < lab$at[1]] <- lab$at[1] + 1e-3
      lab$lab[1] <- paste('<',lab$lab[1],sep='') 
    }
    if(round(max(x$plotdata),1)>lab$at[lla]) { 
      x$plotdata[x$plotdata > lab$at[lla]] <- lab$at[lla] - 1e-3
      lab$lab[lla] <- paste('>',lab$lab[lla],sep='')
    }
  }

  ncol    <- length(lab$cols)
  colours <- lab$cols
  print(norm)
 
  #plot
  #levelplot(plotdata~lon*lat|as.factor(year)*run,x,...,
  levelplot(plotdata~lon*lat|run,x,...,
  #contourplot(plotdata~lon*lat|run,x,...,
  #          region=T,labels=F,contour=F,
            as.table=T,cex=labcex*0.5,
            # scale bar
            col.regions=colours,cuts=ncol-2,
            at=lab$at,
            colorkey=list(space='bottom',height=0.8,width=labcex*0.7,labels=list(labels=lab$lab,at=lab$at,cex=labcex*0.4),tck=labcex*0.5),
            sub=list(eval(parse(text=paste('expression(',lab$name,"   (",lab$unit,")",')',sep=''))),cex=0.5*labcex),
            # axes
            aspect='iso',xlim=c(lims[1:2]),ylim=c(lims[3:4]),xlab=NULL,ylab=NULL,
            scales=list(alternating=F,ce=labcex*0.3),
            xscale.components=function(...,top=F,lim=lims[1:2]){
              axis.ticks(lim=lim,...)},
            yscale.components=function(...,lim=lims[3:4]){
              axis.ticks(lim=lim,y=T,...)},
            strip=stripprint,
            par.strip.text=list(cex=labcex*0.35,lines=1),
            # panel function
            panel=function(...){
              panel.levelplot(...)
              # block weird smearing of levelplot in coastal areas
              # atlantic brasil 
              panel.polygon(x=c(-35,-25,-25,-35),y=c(-10,-10,-4,-4),col='white',border=F)
              # atlantic africa 
              panel.polygon(x=c(-27,-17,-17,-27),y=c(10,10,24,24),col='white',border=F)
              #  madagascar
              panel.polygon(x=c(60,50,50,60),y=c(-10,-10,-20,-20),col='white',border=F)
              # Indus valley 
              if(!is.null(mask)) panel.polygon(x=c(70,55,55,70),y=c(22,22,30,30),col='white',border=F)
              # add world map
              lapply(1:length(world@polygons),plotmap,map=world)
              panel.abline(h=c(0,-23.3,23.3,66.5,-66.5),lty=c(2,3,3,3,3,3),lwd=0.5)
              if(is.null(diff)&is.null(norm)&lab$printsum) {
                panel.text(x=lims[7],y=lims[8],pos=4,labels=ifelse(lab$gsum|lab$gmean,t(global_sum)[panel.number()],''),cex=labcex*0.35)      
                panel.text(x=lims[7]+cs*((lims[2]-lims[1])/25),y=lims[8],pos=4,labels=ifelse(lab$gsum|lab$gmean,parse(text=lab$sunit),''),cex=labcex*0.35)      
              }   
            }) 
}   



axis.ticks <- function ( ... , y=F , ticks=seq(-180,180,20) , ticks2=seq(-180,180,5) ) {
  
  ans    <- xscale.components.default(...)
  ticks2 <- ticks2[!(ticks2 %in% ticks)]
  ans$bottom$ticks$at      <- c(ticks, ticks2)
  ans$bottom$ticks$tck     <- c(rep(-0.5,length(ticks)), rep(-0.2,length(ticks2)))
  lab_ticks                <- seq(-160,160,40) 
  ans$bottom$labels$at     <- lab_ticks
  ans$bottom$labels$labels <- lab_ticks
  ans$bottom$labels$cex    <- labcex
  
  if(y){
    ans    <- yscale.components.default(...)
    ticks  <- c(seq(-80,80,20))
    ticks2 <- seq(-90,90,5)
    ticks2 <- ticks2[!(ticks2 %in% ticks)]
    ans$left$ticks$at      <- c(ticks, ticks2)
    ans$left$ticks$tck     <- c(rep(-0.5,length(ticks)), rep(-0.25,length(ticks2)))
    lab_ticks              <- seq(-80,80,20) 
    ans$left$labels$at     <- lab_ticks
    ans$left$labels$labels <- lab_ticks
    ans$left$labels$cex    <- labcex
  }
  
  ans
}



### END ###
