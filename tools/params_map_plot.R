##########################
#
# Map plotting variables and parameters
#
# AWalker
# Feb 2015
#
#########################

# Colours
red.white   <- colorRampPalette(c('red3','grey90'),space='Lab',bias=10)
blue.white  <- colorRampPalette(c('midnightblue','grey90'),space='Lab',bias=10)
green.white <- colorRampPalette(c('darkgreen','grey90'),space='Lab',bias=10)
brown.white <- colorRampPalette(c('darkorange3','grey90'),space='Lab',bias=10)
nc <- 7 
col.neg     <- c(blue.white(nc+1)[1:nc],'grey90',rev(red.white(nc+1)[1:(nc)]))
col.nbp     <- c(brown.white(nc+1)[1:nc],'grey90',rev(green.white(nc+1)[1:(nc)]))
col.inc     <- c('grey90',rev(topo.colors(10)),'darkred')
col.inc.gpp <- c(col.inc[1:4],'greenyellow',col.inc[6],'green4',col.inc[8:10],'navy',rev(c('orange','darkorange2','red','darkred')))


# Variable parameters
npp <- list(
  name  = 'NPP',
  nsub  = '',
  unit  = 'gCm^-2*y^-1',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,50,125,250,375,500,750,1000,1250,1500,1750,2000),
  cols  = col.inc.gpp
)  
gpp <- list(
  name  = 'GPP',
  nsub  = '',
  unit  = 'gCm^-2*y^-1',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000),
  #at    = c(0,100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500),
  #cols  = viridis(14)
  cols  = col.inc.gpp
)  
nbp <- list(
  name  = 'NBP',
  nsub  = '',
  unit  = 'gCm^-2*y^-1',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  #at    = c(-400,-300,-200,-100,-50,-10,10,50,100,200,300,400),
  at    = c(-300,-250,-200,-150,-100,-50,50,100,150,200,250,300),
  cols  = col.nbp
)  
anlfn <- list(
  name  = 'LeafN',
  nsub  = '',
  unit  = 'gNm^-2',
  sunit = 'PgN',
  gsum  = F,
  gmean = T,
  printsum = F,
  #at    = c(0,0.25,0.5,1,1.5,2,2.5,3,3.5,4,5,6),
  #cols  = col.inc
  at    = c(0,0.5,0.75,1,1.25,1.5,1.75,2,2.250,2.5,2.75,3.0,3.5,4,4.5,5),
  cols  = col.inc.gpp
) 
antlfn <- list(
  name  = 'TopleafN',
  nsub  = '',
  unit  = 'gNm^-2',
  sunit = 'PgN',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2.,2.250,2.500,2.750,3.000,3.500,4.000),
  cols  = col.inc.gpp
)         
evt <- list(
  name  = 'Total_Evapotranspiration',
  nsub  = '',
  unit  = 'kgm^-2*y^-1',
  sunit = 'Eg',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000),
  cols  = col.inc.gpp
) 
trn <- list(
  name  = 'Transpiration',
  nsub  = '',
  unit  = 'kgm^-2*y^-1',
  sunit = 'Eg',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,50,100,200,300,400,500,600,700,800,900,1000,1100,1200),
  cols  = col.inc.gpp
)
scn  <- list(
  name  = '"Total_Soil_Carbon"',
  nsub  = '',
  unit  = 'gCm^-2',
  sunit = 'PCg',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,1000,2500,5000,7500,10000,15000,20000,25000,30000,35000,40000,45000,50000),
  #cols  = col.inc
  cols  = c( col.inc[1] , rev(viridis(11)) )
)
sresp <- list(
  name  = 'Heterotrophic_Respiration',
  nsub  = '',
  unit  = 'gCm^-2*y^-1',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000),
  cols  = col.inc
)
presp <- list(
  name  = 'Canopy_Respiration',
  nsub  = '',
  unit  = 'gCm^-2*y^-1',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,10,25,50,75,100,150,200,250,300,350,400),
  cols  = col.inc
)
mgresp <- list(
  name  = 'Plant_Respiration',
  nsub  = '',
  unit  = 'gCm^-2*y^-1',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,100,250,500,750,1000,1500,2000,2500,3000,3500,4500),
  cols  = col.inc
)
lai <- list(
  name  = 'LAI',
  nsub  = '',
  unit  = 'm^2*" "*m^-2',
  sunit = 'm^2*" "*m^-2',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0,0.5,1:8,9.5,11.5),
  cols  = col.inc
)
anvcmax <- list(
  name  = 'V[cmax]',
  nsub  = 'cmax',
  unit  = 'mu*mol*" "*m^-2*" "*s^-1',
  sunit = 'mu*mol*" "*m^-2*" "*s^-1',
  gsum  = F,
  gmean = T,
  printsum = F,
#  at    = c(0,10,15,20,25,30,35,40,50,60,70,80 ,100,120,140,160),
  at    = c(0,5 ,10,15,20,30,40,50,60,75,90,105,120,140,160,180),
#  at    = c(0,5,10,15,20,25,30,35,40,45,50,55,60,70,80),
#  cols  = col.inc.gpp
  cols  = viridis(15)
)

anjmax <- list(
  name  = 'J[max]',
  nsub  = 'max',
  unit  = 'mu*mol*" "*m^-2*" "*s^-1',
  sunit = 'mu*mol*" "*m^-2*" "*s^-1',
  gsum  = F,
  gmean = T,
  printsum = F,
#  at    = c(0,20,30,40,50,60,70,80,100,120,140,160,200,240,280,320),
  at    = c(0,10,20,30,40,60,80,100,120,150,180,210,240,280,320,360),
  cols  = col.inc.gpp
)
biot <- list(
  name  = '"Vegetation_Carbon"',
  nsub  = '',
  unit  = 'gCm^-2',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,1000,2000,3000,5000,7500,10000,15000,20000,25000,30000,35000),
  cols  = col.inc
)
kg_beta<- list(
  name  = 'Soil_Water_Limitation_Scalar',
  nsub  = '',
  unit  = '',
  sunit = '',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.3),
  cols  = col.inc
)
swr<- list(
  name  = 'SWR',
  nsub  = '',
  unit  = 'Wm^-2',
  sunit = 'Wm^-2',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0,10,20,50,75,100,125,150,175,200,230,260,290),
  cols  = col.inc
)
qtotal<- list(
  name  = 'PAR',
  nsub  = '',
  unit  = 'mol*m^-2*y^-1',
  sunit = 'mol*m^-2*y^-1',
  gsum  = F,
  gmean = T,
  printsum = F,
  #at    = c(0,500,1000,2000,3000,5000,7500,10000,12500,15000,17500,20000),
  #cols  = col.inc
  at    = c(5000,6200,7400,8600,9800,11000,12200,13400,14600,15800,17000,18200,19400,20600),
  cols  = col.inc.gpp
)
tmp <- list(
  name  = 'Temperature',
  nsub  = '',
  unit  = 'degree*C',
  sunit = 'degree*C',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(-20,-15,-10,-5,0,5,10,15,20,25,30,35),
  cols  = c(rev(heat.colors(10))[1:9],'darkred','purple4')
)
prc <- list(
  name  = 'Precipitaton',
  nsub  = '',
  unit  = 'kgm^-2',
  sunit = 'Eg',
  gsum  = T,
  gmean = F,
  printsum = F,
  at    = c(0,100,250,500,750,1000,1500,2000,2500,3000,3500,4000,4500,5000),
  cols  = col.inc.gpp
)
swc <- list(
  name  = 'Soil_Water',
  nsub  = '',
  unit  = 'mm',
  sunit = 'Eg',
  gsum  = T,
  gmean = F,
  printsum = F,
  at    = c(0,100,150,200,225,250,275,300,325,350,375,400),
  cols  = col.inc
)
field_capacity <- list(
  name  = 'Field_Capacity',
  nsub  = '',
  unit  = 'proportion',
  sunit = '',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.4),
  cols  = col.inc
)
wilting_point<- list(
  name  = 'Wilting_Point',
  nsub  = '',
  unit  = 'proportion',
  sunit = '',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3),
  cols  = col.inc
)
cov_C3 <- list(
  name  = 'C3Grass_Cover',
  nsub  = '',
  unit  = 'proportion',
  sunit = 'proportion',
  gsum  = F,
  gmean = T,
  printsum = F,
  at    = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
  cols  = col.inc
)
cov_C4 <- cov_C3
cov_C4$name <- 'C4Grass_Cover'
cov_C3crop <- cov_C3
cov_C3crop$name <- 'C3Crop_Cover'
cov_C4crop <- cov_C3
cov_C4crop$name <- 'C4_Crop_Cover'
cov_Dc_Bl <- cov_C3
cov_Dc_Bl$name <- 'Deciduous_Broadleaf_Cover'
cov_Dc_Nl <- cov_C3
cov_Dc_Nl$name <- 'Deciduous_Needleleaf_Cover'
cov_Ev_Bl <- cov_C3
cov_Ev_Bl$name <- 'Evergreen_Broadleaf_Cover'
cov_Ev_Nl <- cov_C3
cov_Ev_Nl$name <- 'Evergreen_Needleleaf_Cover'
cov_BARE <- cov_C3
cov_BARE$name <- 'Bareground'
fab <- cov_C3
fab$name <- 'Fire_Burnt_Area'


fcn <- list(
  name  = 'Fire_C_losses',
  nsub  = '',
  unit  = 'gCm^-2',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000) /10 ,
  cols  = col.inc.gpp
)  
lulccc <- list(
  name  = 'Land-use_Change_C_losses',
  nsub  = '',
  unit  = 'gCm^-2',
  sunit = 'PgC',
  gsum  = T,
  gmean = F,
  printsum = T,
  at    = c(0,100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000) /10 ,
  cols  = col.inc.gpp
)

