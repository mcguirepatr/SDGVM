library(gdata)

neglatlon <- as.matrix(read.table('lat_lon_formatted.dat'))

neglatlon[,1] <- -neglatlon[,1] #get the negative sign correct for latitude; keep the same sign for longitude
write.fwf(neglatlon, width=c(7,9), nsmall=3, 'lat_lon_formatted2.dat', colnames=F, sep='' )

neglatlon <- as.matrix(read.table('lat_lon_formatted_1deg.dat'))

neglatlon[,1] <- -neglatlon[,1] #get the negative sign correct for latitude; keep the same sign for longitude
write.fwf(neglatlon, width=c(7,9), nsmall=3, 'lat_lon_formatted_1deg2.dat', colnames=F, sep='' )
