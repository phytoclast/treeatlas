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
biomemap <- read_sf('C:/a/geo/ecoregion/biomkuchler2.shp')
biomlist <- read.csv('data/biomlist.csv')
biomsort <- read.csv('data/biomsort.csv')
biomlist$synbiome2 <- biomlist$synbiome

Tw <- raster(paste0(path, 'Tw.tif'))
Twh <- raster(paste0(path, 'Twh.tif'))
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
bedrock <- raster(paste0(path, 'bedrock.tif'))


rasters<-stack(Tw,Twh,Tgs,Tc,Tclx,M,Surplus,Deficit,pAET,slope,sand,SoilpH,hydric,salids,sealevel,bedrock)

biomemap <- st_transform(biomemap, st_crs(SoilpH))
#biomemapmerge <- readRDS('data/biomemapmerge.RDS')
biomemap <- subset(biomemapmerge, select=c("ECO_NAME",  "TYPE",      "AREA",      "BIOME",     "ECO_NUM",   "ECO_ID",    "ECO_SYM",   "CODE","geometry"))
biomemapmerge <- merge(biomemap,unique(biomlist[,c("TYPE","ECO_NAME","synbiome")]),by=c("TYPE","ECO_NAME"), all.x=T)
biomemapmerge <- merge(biomemapmerge,unique(biomlist[biomlist$TYPE %in% 'unk',c("ECO_NAME","synbiome2")]),by='ECO_NAME', all.x=T)
biomemapmerge[is.na(biomemapmerge$synbiome),]$synbiome <- biomemapmerge[is.na(biomemapmerge$synbiome),]$synbiome2
biomemapmerge <- subset(biomemapmerge, !is.na(synbiome))
biomemapmerge <- merge(biomemapmerge, biomsort, by='synbiome')
#saveRDS(biomemapmerge,'data/biomemapmerge.RDS')

biomemapmerge$samples <-  floor(biomemapmerge$AREA/(50000^2)+2)
#kuchsample1 <- st_sample(biomemapmerge, biomemapmerge$samples, type = "random", exact = FALSE)
#kuchsample2 <- st_sample(biomemapmerge, 30000, type = "random", exact = FALSE)
kuchsample <- c(kuchsample1,kuchsample2)
kuchsample <- readRDS('data/kuchsample.RDS')
#saveRDS(kuchsample, 'data/kuchsample.RDS')

kuchsample <- st_sf(kuchsample)


kuchjoin <- st_join(st_sf(kuchsample),biomemapmerge[,c('sort','synbiome')])
rasterjoin <- raster::extract(rasters, kuchjoin)
kuchsample <- cbind(kuchjoin, rasterjoin)
#plot(gdem, axes=FALSE, legend=FALSE, main='Random Sampling with NA Removal')
#plot(st_geometry(mlrasample), add=T)

#classification tree packages
library(randomForest)
library(rpart)
library(rpart.plot)

#randomforest
selectBiome<-subset(kuchsample, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) &!is.na(slope) 
                    &(synbiome %in% c('Rock and Ice') | Tgs > 0) 
                    &(!synbiome %in% c('Rock and Ice', 'Tundra') | Tgs < 9) 
                    &(!synbiome %in% c('Warm Xeric Woodland','Mediterranean Scrub') | M < 1.41) 
                    &(!synbiome %in% c('Subtropical Mixed Forest') | M > 0.71) 
                    &(synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh','Tundra') | hydric < 67) 
                    &(!synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh') | hydric > 5) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland') | Tgs > 3) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland', 'Boreal/Subalpine Forest') | Tgs > 9|Tc > 0))

#saveRDS(selectBiome, 'data/selectBiome.RDS')
Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$sort, selectBiome$synbiome), FUN=length)
names(Spwts)<-c('sort','synbiome','x')
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by=c('sort','synbiome'))
rf <- randomForest(as.factor(sort) ~  Tw+Twh+Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel+bedrock, 
                   data=selectBiome, importance=TRUE, ntree=500,na.action=na.omit)
# Make plot  other params to try: maxnodes=512,mtry=10, classwt=Spwts$myweights, ,
#statistical summary
varImpPlot(rf)

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/synbiomebedrock.tif") 

plot(vegmaprf)   
#rpart
selectBiome<-subset(kuchsample, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) &!is.na(slope) 
                      )

Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$sort, selectBiome$synbiome), FUN=length)
names(Spwts)<-c('sort','synbiome','x')
Spwts$myweights<- (10000/(Spwts$x/1+100))/10
selectBiome<-merge(selectBiome,Spwts,by=c('sort','synbiome'))

model <- rpart(as.factor(synbiome)  ~ Tw+Twh+Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, data= selectBiome, weights = selectBiome$myweights, method = "class", cp = 0.01,
               maxdepth = 5)

png(filename="output/kuchlerclass_NorthAmerica.png",width = 10, height = 3, units = 'in', res = 600)

rpart.plot(model, extra=108) # Make plot

dev.off()

