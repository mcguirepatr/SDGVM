##############################
#
# Process CO2 data for SDGVM 
#
# AWalker
# Oct 2013
#
##############################

wd <- '/mnt/disk3/Databases/CO2/'
setwd(wd)

conv_ppm_to_Pa <- function(ppm, press=1.013e5 ) {
  ppm * press/1e6
}

# Spline interpolated SCRIPPS dataset
ifile  <- 'SCRIPPS_spline_merged_ice_core_yearly'
mydata <- read.csv(paste(ifile,'csv',sep='.'),skip=28,header=F)
dim(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[2]])
write.table(mydata[,c(1,3)],paste(ifile,'dat',sep='.'),quote=F,row.names=F,col.names=F)


# RCP 8.5 Meinshusen from PIK
ifile <- 'RCP85_MIDYEAR_CONCENTRATIONS.DAT'
mydata<- read.table(ifile,skip=38,header=T)
head(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[4]])
names(mydata)
write.table(mydata[,c(1,37)],'RCP85_MIDYEAR_CONCENTRATIONS_Pa.DAT',quote=F,row.names=F,col.names=F)


# TRENDY 2016
ifile <- 'global_co2_ann_1860_2015_TRENDY.txt'
mydata<- read.table(ifile,header=F)
head(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[2]])
names(mydata)
write.table(mydata[,c(1,3)],'global_co2_ann_1860_2015_TRENDY.dat',quote=F,row.names=F,col.names=F)

# TRENDY 2016 ENSO
ifile <- 'global_co2_ann_1860_2016_TRENDY.txt'
mydata<- read.table(ifile,header=F)
head(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[2]])
names(mydata)
write.table(mydata[,c(1,3)],'global_co2_ann_1860_2016_TRENDY.dat',quote=F,row.names=F,col.names=F)

# TRENDY 2017
ifile <- 'global_co2_ann_1700_2017_TRENDY.txt'
mydata<- read.table(ifile,header=F)
head(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[2]])
names(mydata)
mydata1 <- mydata
mydata1[,3] <- mydata1[1,3] 
write.table(mydata[,c(1,3)],'global_co2_ann_1700_2017_TRENDY.dat',quote=F,row.names=F,col.names=F)
write.table(mydata1[,c(1,3)],'global_co2_ann_1700_2017_TRENDY_1700.dat',quote=F,row.names=F,col.names=F)


# TRENDY 2019
ifile  <- 'global_co2_ann_1700_2019_TRENDY.txt'
mydata <- read.table(ifile,header=F)
head(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[2]])
names(mydata)
mydata1 <- mydata
mydata1[,3] <- mydata1[1,3] 
write.table(mydata[,c(1,3)],'global_co2_ann_1700_2019_TRENDY.dat',quote=F,row.names=F,col.names=F)


# TRENDY 2022 (2021 dataset)
ifile  <- 'global_co2_ann_1700_2021_TRENDY.txt'
mydata <- read.table(ifile,header=F)
head(mydata)
mydata$Pa <- conv_ppm_to_Pa(mydata[[2]])
names(mydata)
mydata1 <- mydata
mydata1[,3] <- mydata1[1,3] 
write.table(mydata[,c(1,3)],'global_co2_ann_1700_2021_TRENDY.dat',quote=F,row.names=F,col.names=F)
write.table(mydata1[,c(1,3)],'global_co2_ann_1700_2021_TRENDY_1700.dat',quote=F,row.names=F,col.names=F)



### END ###