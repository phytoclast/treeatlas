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
altveg <- read_sf('data/biotic_comm_la.shp')
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
clay <- raster(paste0(path, 'clay.tif'))


rasters<-stack(Tw,Twh,Tgs,Tc,Tclx,M,Surplus,Deficit,pAET,slope,sand,SoilpH,hydric,salids,sealevel,bedrock,clay)

altveg <- st_transform(altveg, st_crs(SoilpH))
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
#kuchsample2 <- st_sample(biomemapmerge, 35000, type = "random", exact = FALSE)
#savesample <- c(kuchsample1,kuchsample2)
#saveRDS(savesample, 'data/kuchsample2.RDS')
#kuchsample <- c(kuchsample1,kuchsample2)
#kuchsampletest <- readRDS('data/kuchsample.RDS')
kuchsample <- readRDS('data/kuchsample2.RDS')
#saveRDS(kuchsample, 'data/kuchsample.RDS')

kuchsample <- st_sf(kuchsample)


kuchjoin <- st_join(st_sf(kuchsample),biomemapmerge[,c('synbiome')])
kuchjoin <- st_join(st_sf(kuchjoin),altveg[,c('biome')])
kuchjoin$synbiome2 <- as.character(kuchjoin$synbiome)
kuchjoin[kuchjoin$biome %in% c('subalpine forest','taiga') &
           kuchjoin$synbiome %in% c('Temperate Coniferous Forest','Temperate Mixed Forest','Temperate Deciduous Forest'),]$synbiome2 <- 'Boreal/Subalpine Forest'
kuchjoin[kuchjoin$biome %in% c('tundra')&
           kuchjoin$synbiome %in% c('Temperate Coniferous Forest','Temperate Mixed Forest'),]$synbiome2 <- 'Tundra'
kuchjoin[kuchjoin$biome %in% c('ice'),]$synbiome2 <- 'Rock and Ice'
kuchjoin[kuchjoin$biome %in% c('paramo'),]$synbiome2 <- 'Montane Grassland'
kuchjoin[kuchjoin$biome %in% c('cloud forest'),]$synbiome2 <- 'Tropical Montane Forest'
kuchjoin[kuchjoin$biome %in% c('thornscrub', 'semi-evergreen forest','semi-deciduous forest','dry deciduous forest') &
           kuchjoin$synbiome %in% c('Tropical Moist Forest'),]$synbiome2 <- 'Tropical Dry Forest'
#kuchjoin[kuchjoin$biome %in% c('semidesert grassland','shortgrass prairie') &
#           kuchjoin$synbiome %in% c('Warm Desert'),]$synbiome2 <- 'Semidesert Grassland'
kuchjoin[kuchjoin$biome %in% c('evergreen woodland','conifer woodland') &
           kuchjoin$synbiome %in% c('Warm Desert'),]$synbiome2 <- 'Warm Xeric Woodland'
kuchjoin[kuchjoin$biome %in% c('thornscrub') &
           kuchjoin$synbiome %in% c('Warm Desert','Tropical Dry Forest','Tropical/Subtropical Savanna'),]$synbiome2 <- 'Tropical/Subtropical Thornscrub'
kuchjoin[kuchjoin$biome %in% c('dry deciduous forest')&
           kuchjoin$synbiome %in% c('Warm Desert'),]$synbiome2 <- 'Tropical/Subtropical Thornscrub'
kuchjoin[kuchjoin$biome %in% c('chaparral','coastal scrub'),]$synbiome2 <- 'Chaparral'
kuchjoin[kuchjoin$biome %in% c('warm swamp'),]$synbiome2 <- 'Warm Swamp'
kuchjoin[kuchjoin$synbiome %in% c('Warm Xeric Woodland') &
           kuchjoin$biome %in% c('montane conifer forest'),]$synbiome2 <- 'Subtropical Mixed Forest'
kuchjoin[kuchjoin$synbiome %in% c('Subtropical Mixed Forest') &
           kuchjoin$biome %in% c('shortgrass prairie', 'semidesert grassland','desertscrub','evergreen woodland','dry deciduous forest')
         ,]$synbiome2 <- 'Warm Xeric Woodland'
kuchjoin$synbiome <- kuchjoin$synbiome2
kuchjoin <- merge(kuchjoin,biomsort[,c('sort','synbiome')], by='synbiome')
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
                    #&(!synbiome %in% c('Warm Xeric Woodland','Mediterranean Scrub') | M < 1.41) 
                    #&(!synbiome %in% c('Subtropical Mixed Forest') | M > 0.71) 
                    &(synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh','Tundra') | hydric < 67) 
                    &(!synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh') | hydric > 15) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland') | Tgs > 3) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland', 'Boreal/Subalpine Forest') | Tgs > 9|Tc > 0))
#selectBiome<-subset(kuchsample, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) &!is.na(slope) 
#                    &(synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh','Tundra') | hydric < 67) 
#                    &(!synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh') | hydric > 5) 
#)

#saveRDS(selectBiome, 'data/selectBiome.RDS')
Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$sort, selectBiome$synbiome), FUN=length)
names(Spwts)<-c('sort','synbiome','x')
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by=c('sort','synbiome'))
rf <- randomForest(as.factor(sort) ~  Tw+Twh+Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+
                     hydric+salids+sealevel+bedrock+clay, 
                   data=selectBiome, importance=TRUE, ntree=200,na.action=na.omit)
# Make plot  other params to try: maxnodes=512,mtry=10, classwt=Spwts$myweights, ,
#statistical summary
varImpPlot(rf)

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/synbiomeclay2.tif") 




















######################
 
#rpart
selectBiome<-subset(kuchsample, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) &!is.na(slope) 
                      )

Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$sort, selectBiome$synbiome), FUN=length)
names(Spwts)<-c('sort','synbiome','x')
Spwts$myweights<- (10000/(Spwts$x/1+100))/10
selectBiome<-merge(selectBiome,Spwts,by=c('sort','synbiome'))
selectBiome$Cindex <- pmin(selectBiome$Tc, selectBiome$Tclx+15)
model <- rpart(as.factor(synbiome)  ~ Tgs+Cindex+M+Surplus+Deficit+pAET+slope+sand+hydric+salids+sealevel, data= selectBiome, weights = selectBiome$myweights, method = "class", cp = 0.01,
               maxdepth = 5)

png(filename="output/kuchlerclass_NorthAmerica.png",width = 10, height = 3, units = 'in', res = 600)

rpart.plot(model, extra=108) # Make plot

dev.off()


#test randomforest

kuchsampletest <- st_sf(kuchsampletest)


kuchjoin <- st_join(st_sf(kuchsampletest),biomemapmerge[,c('sort','synbiome')])
rasterjoin <- raster::extract(rasters, kuchjoin)
kuchsampletest <- cbind(kuchjoin, rasterjoin)

selectBiometest<-st_drop_geometry(subset(kuchsampletest, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) &!is.na(slope) 
                    &(synbiome %in% c('Rock and Ice') | Tgs > 0) 
                    &(!synbiome %in% c('Rock and Ice', 'Tundra') | Tgs < 9) 
                    &(!synbiome %in% c('Warm Xeric Woodland','Mediterranean Scrub') | M < 1.41) 
                    &(!synbiome %in% c('Subtropical Mixed Forest') | M > 0.71) 
                    &(synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh','Tundra') | hydric < 67) 
                    &(!synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh') | hydric > 5) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland') | Tgs > 3) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland', 'Boreal/Subalpine Forest') | Tgs > 9|Tc > 0)))
selectBiometest$predicted <- predict(rf,selectBiometest) 
library(caret) #confusion matrix
confusionMatrix(as.factor(selectBiometest$predicted), as.factor(selectBiometest$sort))
plot(rf)