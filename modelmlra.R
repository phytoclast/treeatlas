library(sf)
library(raster)
library(rgdal)
library(velox)
library(lattice)
library(dplyr)
library(ggplot2)

#path = 'C:/a/geo/model/global/'
path = 'nam5k/'
treelist <- read.delim('treefiles.txt')
bt_frost <- read.delim('bt_frost2.txt')
mlramap <-read_sf('mlra/mlra2018.shp')

Tclx <- raster(paste0(path, 'Tclx.tif'))
Tgs <- raster(paste0(path, 'Tgs.tif'))
Tc <- raster(paste0(path, 'Tc.tif'))
M <- raster(paste0(path, 'M.tif'))
MAP <- raster(paste0(path, 'MAP.tif'))
Deficit <- raster(paste0(path, 'Deficit.tif'))
Surplus <- raster(paste0(path, 'Surplus.tif'))
pAET <- raster(paste0(path, 'pAET.tif'))
slope <- raster(paste0(path, 'slope.tif'))
sand <- raster(paste0(path, 'sand.tif'))
SoilpH <- raster(paste0(path, 'SoilpH.tif'))
hydric <- raster(paste0(path, 'hydric.tif'))
salids <- raster(paste0(path, 'salids.tif'))
water100k <- raster(paste0(path, 'water100k.tif'))
sealevel <- raster(paste0(path, 'sealevel.tif'))
gdem <- raster(paste0(path, 'gdem.tif'))
#sealevel <- (gdem*(gdem > 0)/(gdem+5)-1)*-1*water100k


#M <- MAP/(Deficit + MAP - Surplus)
#writeRaster(sealevel,paste0(path, 'sealevel.tif'), COMPRESS='LZW', overwrite=T)
#writeRaster(M,paste0(path, 'M.tif'), COMPRESS='LZW', overwrite=T)
mlramap <- st_transform(mlramap, st_crs(SoilpH))
rasters<-stack(Tgs,Tc,Tclx,M,Surplus,Deficit,pAET,slope,sand,SoilpH,hydric,salids,sealevel)

mlrasample <- st_sample(mlramap, 50000, type = "random", exact = FALSE)
mlrasample <- st_sf(mlrasample)
mlrajoin <- st_join(st_sf(mlrasample),mlramap[,c('MLRA')])
rasterjoin <- raster::extract(rasters, mlrajoin)
mlrasample <- cbind(mlrajoin, rasterjoin)
#plot(gdem, axes=FALSE, legend=FALSE, main='Random Sampling with NA Removal')
#plot(st_geometry(mlrasample), add=T)
     
#classification tree packages
library(randomForest)
library(rpart)
library(rpart.plot)

#randomforest
selectBiome<-subset(mlrasample, MLRA !='' & !is.na(sand) & !is.na(M))
Spwts<-aggregate(x=selectBiome$MLRA, by=list(selectBiome$MLRA), FUN=length)
names(Spwts)<-c("MLRA","x")
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by='MLRA')
rf <- randomForest(as.factor(MLRA) ~  Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, data=selectBiome, classwt=Spwts$myweights, importance=TRUE, ntree=500, na.action=na.omit, maxnodes = 512)
# Make plot  other params to try: maxnodes=64,mtry=10,
rf#statistical summary
varImpPlot(rf)

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/mlra512.tif") 

plot(vegmaprf)   
#rpart
selectBiome<- Biome_extract
Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$synbiome), FUN=length)
names(Spwts)<-c("synbiome","x")
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by='synbiome')

model <- rpart(as.factor(synbiome)  ~ Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, data= selectBiome, weights = selectBiome$wt, method = "class", cp = 0.01,
               maxdepth = 6)

png(filename="output/kuchlerclass_NorthAmerica.png",width = 10, height = 3, units = 'in', res = 600)

rpart.plot(model, extra=108) # Make plot

dev.off()

#randomforest2
selectBiome<- subset(Biome_extract, !is.na(synth)& synth != '')
Spwts<-aggregate(x=selectBiome$synth, by=list(selectBiome$synth), FUN=length)
names(Spwts)<-c("synth","x")
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by='synth')
rf <- randomForest(as.factor(as.character(synth)) ~  Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, data=selectBiome, classwt=selectBiome$wt, importance=TRUE, ntree=200, na.action=na.omit )
# Make plot
rf#statistical summary
varImpPlot(rf)

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/kuchlermodelrf2.tif") 

plot(vegmaprf)   
               

               