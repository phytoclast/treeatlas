library(sf)
library(raster)
library(rgdal)
library(velox)
library(lattice)
library(dplyr)
treelist <- read.delim('treefiles.txt')
Tclx <- raster('nam/Tclx.tif')
Tgs <- raster('nam/Tgs.tif')
Tc <- raster('nam/Tc.tif')
MAP <- raster('nam/MAP.tif')
Deficit <- raster('nam/Deficit.tif')
Surplus <- raster('nam/Surplus.tif')
pAET <- raster('nam/pAET.tif')
M <- MAP/(Deficit + MAP - Surplus)
n <- rownames(treelist[grepl('cyrirace', treelist$name),])
spp<-read_sf(treelist[n,1])
spp2 <- 	st_transform(spp, st_crs(Tclx))
spp2 <- spp2[spp2$CODE %in% '1',]

plot(Tclx, main=treelist[n,4])
plot(st_geometry(spp2), add=T)
top = 0.98
bot = 0.02

#Tclx----
vx <- velox(Tclx)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tclx', ylab='Density', main='')
Tclx_max <- quantile(e[,2], top)
Tclx_min <- quantile(e[,2], bot)

#Tgs----
vx <- velox(Tgs)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tgs', ylab='Density', main='')
Tgs_max <- quantile(e[,2], top)
Tgs_min <- quantile(e[,2], bot)

#Tc----
vx <- velox(Tc)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tc', ylab='Density', main='')
Tc_max <- quantile(e[,2], top)
Tc_min <- quantile(e[,2], bot)

#M----
vx <- velox(M)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='M', ylab='Density', main='')
M_max <- quantile(e[,2], top)
M_min <- quantile(e[,2], bot)

#Surplus----
vx <- velox(Surplus)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Surplus', ylab='Density', main='')
Surplus_max <- quantile(e[,2], top)
Surplus_min <- quantile(e[,2], bot)

#Deficit----
vx <- velox(Deficit)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Deficit', ylab='Density', main='')
Deficit_max <- quantile(e[,2], top)
Deficit_min <- quantile(e[,2], bot)

#pAET----
vx <- velox(pAET)
e <- vx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

pAET_max <- quantile(e[,2], top)
pAET_min <- quantile(e[,2], bot)


envelop <- (Tclx <= Tclx_max)*(Tclx >= Tclx_min)*
  (Tgs <= Tgs_max)*(Tgs >= Tgs_min)*
  (M <= M_max)*(M >= M_min)*
  (pAET <= pAET_max)*(pAET >= pAET_min)*
  (Tc <= Tc_max)*(Tc >= Tc_min)*
  (Deficit <= Deficit_max)*(Deficit >= Deficit_min)*
  (Surplus <= Surplus_max)*(Surplus >= Surplus_min)
plot(envelop)
plot(st_geometry(spp2), add=T)