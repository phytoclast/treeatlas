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
BPS <- raster(paste0(path, 'BPS.tif'))


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

#import tree ranges to filter subtropical forest
magnolia <-read_sf('tree_atlas_dd83/l-o/magngran.shp')
liveoak <-read_sf('tree_atlas_dd83/q/quervirg.shp')
magnolia <- 	st_transform(magnolia, st_crs(Tclx))
magnolia$tree1 <- magnolia$CODE
liveoak <- 	st_transform(liveoak, st_crs(Tclx))
liveoak$tree2 <- liveoak$CODE

#bring in fresh copy of vegetation crosswalks
biomlist <- read.csv('data/biomlist.csv')
biomsort <- read.csv('data/biomsort.csv')
biomlist$synbiome2 <- biomlist$synbiome
bpsgroups <- read.csv('data/bpsgroups.csv')

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
kuchjoin <- st_join(st_sf(kuchjoin),magnolia[,c('tree1')])
kuchjoin <- st_join(st_sf(kuchjoin),liveoak[,c('tree2')])
#add layer to make raster stack for proper extract
BPS$BPS2 <- 1

bpsval <- (raster::extract(BPS, kuchjoin))
kuchsample <- cbind(kuchjoin, bpsval)
kuchsample <- subset(kuchsample, select = -BPS2)
biomsort2 <- biomsort
colnames(biomsort2)[1] <- 'sort2'
kuchsample <- merge(kuchsample, bpsgroups, by.x = 'BPS', by.y = 'VALUE', all.x=TRUE)
#kuchsample <- merge(kuchsample, biomsort2, by.x = 'legend', by.y = 'synbiome', all.x=TRUE)

kuchsample$synbiome2 <- as.character(kuchsample$synbiome)
kuchsample[kuchsample$biome %in% c('subalpine forest','taiga') &
           kuchsample$synbiome %in% c('Temperate Coniferous Forest','Temperate Mixed Forest','Temperate Deciduous Forest'),]$synbiome2 <- 'Boreal/Subalpine Forest'
kuchsample[kuchsample$biome %in% c('tundra')&
           kuchsample$synbiome %in% c('Temperate Coniferous Forest','Temperate Mixed Forest'),]$synbiome2 <- 'Tundra'
kuchsample[kuchsample$biome %in% c('ice'),]$synbiome2 <- 'Rock and Ice'
kuchsample[kuchsample$biome %in% c('paramo'),]$synbiome2 <- 'Montane Grassland'
kuchsample[kuchsample$biome %in% c('cloud forest'),]$synbiome2 <- 'Tropical Montane Forest'
kuchsample[kuchsample$biome %in% c('thornscrub', 'semi-evergreen forest','semi-deciduous forest','dry deciduous forest') &
           kuchsample$synbiome %in% c('Tropical Moist Forest'),]$synbiome2 <- 'Tropical Dry Forest'
#kuchsample[kuchsample$biome %in% c('semidesert grassland','shortgrass prairie') &
#           kuchsample$synbiome %in% c('Warm Desert'),]$synbiome2 <- 'Semidesert Grassland'
kuchsample[kuchsample$biome %in% c('evergreen woodland','conifer woodland') &
           kuchsample$synbiome %in% c('Warm Desert'),]$synbiome2 <- 'Warm Xeric Woodland'
kuchsample[kuchsample$biome %in% c('thornscrub') &
           kuchsample$synbiome %in% c('Warm Desert','Tropical Dry Forest','Tropical/Subtropical Savanna'),]$synbiome2 <- 'Tropical/Subtropical Thornscrub'
kuchsample[kuchsample$biome %in% c('dry deciduous forest')&
           kuchsample$synbiome %in% c('Warm Desert'),]$synbiome2 <- 'Tropical/Subtropical Thornscrub'
kuchsample[kuchsample$biome %in% c('chaparral','coastal scrub'),]$synbiome2 <- 'Chaparral'
kuchsample[kuchsample$biome %in% c('warm swamp'),]$synbiome2 <- 'Warm Swamp'
kuchsample[kuchsample$synbiome %in% c('Warm Xeric Woodland') &
           kuchsample$biome %in% c('montane conifer forest'),]$synbiome2 <- 'Subtropical Mixed Forest'
kuchsample[kuchsample$synbiome %in% c('Subtropical Mixed Forest') &
           kuchsample$biome %in% c('shortgrass prairie', 'semidesert grassland','desertscrub','evergreen woodland','dry deciduous forest')
         ,]$synbiome2 <- 'Warm Xeric Woodland'
kuchsample$synbiome <- kuchsample$synbiome2
#savekuchsample <- kuchsample

#separate the BPS points and add them as new points to hedge the model
kuchsample$legend <- as.character(kuchsample$legend)
bpspoints <- kuchsample[!is.na(kuchsample$legend),]
#kuchsample <- kuchsample[is.na(kuchsample$legend),]
#try to amplify subtropical forests near riparian areas
bpspoints[bpspoints$legend %in% 'Tropical Pine Forest' & bpspoints$synbiome2 %in% 'Warm Swamp',]$legend <-'Subtropical Mixed Forest'
bpspoints[is.na(bpspoints$tree1),]$tree1 <-0
bpspoints[is.na(bpspoints$tree2),]$tree2 <-0
selection <- bpspoints[(bpspoints$legend %in% 'Tropical Pine Forest' | bpspoints$synbiome2 %in% 'Subtropical Mixed Forest'),]
bpspoints[(bpspoints$legend %in% 'Temperate Deciduous Forest' | bpspoints$synbiome2 %in% 'Subtropical Mixed Forest')&
            (bpspoints$tree1 >=1 | bpspoints$tree2 >=1),]$legend <-'Subtropical Mixed Forest'
bpspoints[(bpspoints$legend %in% 'Subtropical Mixed Forest' | bpspoints$synbiome2 %in% 'Temperate Oak Forest')&
            (bpspoints$tree1 == 0 & bpspoints$tree2 == 0) & !grepl('California',bpspoints$BPS_NAME),]$legend <-'Temperate Oak Forest'


#copy BPS to main vegetation field
bpspoints$synbiome <- bpspoints$legend

#fold points back in
kuchsample <- rbind(kuchsample,bpspoints)
kuchsample <- merge(kuchsample,biomsort[,c('sort','synbiome')], by='synbiome')
#original replacement with BPS
#kuchsample[!is.na(kuchsample$legend),]$synbiome <- kuchsample[!is.na(kuchsample$legend),]$legend



rasterjoin <- raster::extract(rasters, kuchsample)



kuchsample1 <- cbind(kuchsample, rasterjoin)


#plot(gdem, axes=FALSE, legend=FALSE, main='Random Sampling with NA Removal')
#plot(st_geometry(mlrasample), add=T)

#classification tree packages
library(randomForest)
library(rpart)
library(rpart.plot)

#randomforest
selectBiome<-subset(kuchsample1, synbiome !='' & !is.na(sand) & !is.na(M) &!is.na(salids) &!is.na(slope) 
                    &(synbiome %in% c('Rock and Ice') | Tgs > 0) 
                    &(!synbiome %in% c('Rock and Ice', 'Tundra') | Tgs < 9) 
                    #&(!synbiome %in% c('Warm Xeric Woodland','Mediterranean Scrub') | M < 1.41) 
                    #&(!synbiome %in% c('Subtropical Mixed Forest') | M > 0.71) 
                    &(synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh','Tundra') | hydric < 67) 
                    &(!synbiome %in% c('Cool Swamp','Warm Swamp','Freshwater Marsh','Salt Marsh') | hydric > 15) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland') | Tgs > 3) 
                    & (synbiome %in% c('Rock and Ice', 'Tundra', 'Montane Grassland', 'Boreal/Subalpine Forest') | Tgs > 9|Tc > 0))


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

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/synbiomeBPS.tif") 













##############################################
plot(vegmaprf)   
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