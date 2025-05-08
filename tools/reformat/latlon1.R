library(gdata)

neglatlon <- as.matrix(read.table('/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/CRUJRA2021/SDGVMdata/monthly/0.5deg_v2a5/lat_lon_formatted.dat'))

neglatlon[,1] <- -neglatlon[,1] #get the negative sign correct for latitude; keep the same sign for longitude
write.fwf(neglatlon, width=c(7,9), nsmall=3, '/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/CRUJRA2021/SDGVMdata/monthly/0.5deg_v2a10/lat_lon_formatted.dat', colnames=F, sep='' )

neglatlon <- as.matrix(read.table('/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/CRUJRA2021/SDGVMdata/monthly/0.5deg_v2a5/lat_lon_formatted_1deg.dat'))

neglatlon[,1] <- -neglatlon[,1] #get the negative sign correct for latitude; keep the same sign for longitude
write.fwf(neglatlon, width=c(7,9), nsmall=3, '/gws/nopw/j04/nexcs/pmcguire/TRENDYv10/db/CRUJRA2021/SDGVMdata/monthly/0.5deg_v2a10/lat_lon_formatted_1deg.dat', colnames=F, sep='' )
