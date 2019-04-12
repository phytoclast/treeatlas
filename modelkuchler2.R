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
#kuchlermap <-read_sf('data/kuchler.shp')
biomemap <- read_sf('C:/a/geo/ecoregion/biomkuchler.shp')
biomlist <- read.csv('data/biomlist.csv')
biomlist$synbiome2 <- biomlist$synbiome

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

rasters<-stack(Tgs,Tc,Tclx,M,Surplus,Deficit,pAET,slope,sand,SoilpH,hydric,salids,sealevel)

biomemap <- st_transform(biomemap, st_crs(SoilpH))

biomemapmerge <- merge(biomemap,unique(biomlist[,c("TYPE","synbiome")]),by='TYPE', all.x=T)
biomemapmerge <- merge(biomemapmerge,unique(biomlist[biomlist$TYPE %in% 'unk',c("ECO_NAME","synbiome2")]),by='ECO_NAME', all.x=T)
biomemapmerge[is.na(biomemapmerge$synbiome),]$synbiome <- biomemapmerge[is.na(biomemapmerge$synbiome),]$synbiome2
biomemapmerge <- subset(biomemapmerge, !is.na(synbiome))
#saveRDS(biomemapmerge,'data/biomemapmerge.RDS')
biomemapmerge$samples <- floor(biomemapmerge$AREA/(50000^2)+2)
kuchsample <- st_sample(biomemapmerge, biomemapmerge$samples, type = "random", exact = FALSE)
kuchsample <- st_sf(kuchsample)
kuchjoin <- st_join(st_sf(kuchsample),biomemapmerge[,c('synbiome')])
rasterjoin <- raster::extract(rasters, kuchjoin)
kuchsample <- cbind(kuchjoin, rasterjoin)
#plot(gdem, axes=FALSE, legend=FALSE, main='Random Sampling with NA Removal')
#plot(st_geometry(mlrasample), add=T)

#classification tree packages
library(randomForest)
library(rpart)
library(rpart.plot)

#randomforest
selectBiome<-subset(kuchsample, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) )
#saveRDS(selectBiome, 'data/selectBiome.RDS')
Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$synbiome), FUN=length)
names(Spwts)<-c("synbiome","x")
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by='synbiome')
rf <- randomForest(as.factor(synbiome) ~  Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, 
                   data=selectBiome, importance=TRUE, ntree=500, na.action=na.omit)
# Make plot  other params to try: maxnodes=64,mtry=10, classwt=Spwts$myweights, ,
rf#statistical summary
varImpPlot(rf)

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/synbiome.tif") 

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
               

               