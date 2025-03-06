 library(lattice)
 lat_res <- 360
 time    <- 1901:2022
 #setwd('/work/scratch-pw2/pmcguire/TRENDYv12/db/CRUJRA2022/SDGVMdata/monthly/0.5deg')
 testdf <- as.matrix(read.table('/work/scratch-pw2/pmcguire/TRENDYv12/db/clim/global/crujra_2023/30min/tmp_byyear.dat'))
 
 class(testdf)
 dim(testdf)
 dim(testdf)[1]/length(time)

 m1 <- matrix(testdf[,7],nrow=67420)
 m1[12000,]
 m1[22000,]
 m1[32000,]
 m1[42000,]
 m1[52000,]
# m1 <- matrix(testdf[,1],nrow=67420)
#levelplot(m1,ylim=c(lat_res,1),col.regions=rev(heat.colors(50)),contour=F,
#          panel=function(...){
#            panel.fill(col='black')
#            panel.levelplot(...)
#          })
#matrix(1:20,nrow=4)
