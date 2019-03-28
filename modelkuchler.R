library(sf)
library(raster)
library(rgdal)
library(velox)
library(lattice)
library(dplyr)
library(ggplot2)

#path = 'C:/a/geo/model/global/'
path = 'nam/'
treelist <- read.delim('treefiles.txt')
bt_frost <- read.delim('bt_frost2.txt')
kuchlermap <-read_sf('data/kuchler.shp')

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
plot(sealevel)
Cindex <- min(Tclx+15, Tc)
#M <- MAP/(Deficit + MAP - Surplus)
#writeRaster(sealevel,paste0(path, 'sealevel.tif'), COMPRESS='LZW', overwrite=T)

Biomeclimate <- readRDS('C:/workspace2/PrepareClimateData/data/bigRadBiomeclimate.RDS')
Biomeclimate <- readRDS('C:/workspace2/modelmap/data/bigRadBiomeclimate.RDS')
Biomeclimate <-  subset(Biomeclimate, Norm %in% '1990')
Biomeclimate$MAP <- apply(Biomeclimate[,c('p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12')], MARGIN = 1, FUN='sum')
Biomeclimate$MAAT <- apply(Biomeclimate[,c('t01','t02','t03','t04','t05','t06','t07','t08','t09','t10','t11','t12')], MARGIN = 1, FUN='mean')
Biomeclimate$M <- Biomeclimate$MAP / (Biomeclimate$Deficit + Biomeclimate$MAP - Biomeclimate$Surplus)+0.0001

Biomeclimate$Cindex <- pmin(Biomeclimate$Tc , Biomeclimate$Tclx + 15)
Biomeclimate$alpine <- ifelse(Biomeclimate$BIOME %in% c('11'), 1, 
                              ifelse(Biomeclimate$BIOME %in% c('6'),2,0 ))

cor(Biomeclimate[Biomeclimate$alpine > 0 & !is.na(Biomeclimate$alpine),c(
  'alpine', 'Tw', 'Twh', 'Tg', 'MAAT', 'Tc', 'Tcl', 'Tclx', 'Cindex'
)])


Biomeclimate$boreal <- ifelse(Biomeclimate$BIOME %in% c('6'), 1, 
                              ifelse(Biomeclimate$BIOME %in% c('4')&Biomeclimate$Latitude >= 45,2,0 ))
cor(Biomeclimate[Biomeclimate$boreal > 0 & !is.na(Biomeclimate$boreal),c(
  'boreal', 'Tw', 'Twh', 'Tg', 'MAAT', 'Tc', 'Tcl', 'Tclx', 'Cindex'
)])

Biomeclimate$temperate <- ifelse(Biomeclimate$BIOME %in% c('4') & !grepl('evergreen',Biomeclimate$ECO_NAME) & Biomeclimate$Latitude < 45 & Biomeclimate$BIOME > -45, 1, 
                                 ifelse(Biomeclimate$BIOME %in% c('1')&!grepl('subtropical',Biomeclimate$ECO_NAME),2,0 ))

cor(Biomeclimate[Biomeclimate$temperate > 0 & !is.na(Biomeclimate$temperate),c(
  'temperate', 'Tw', 'Twh', 'Tg', 'MAAT', 'Tc', 'Tcl', 'Tclx', 'Cindex'
)])


Biomeclimate$y <- Biomeclimate$Latitude
Biomeclimate$x <- Biomeclimate$Longitude
Biome_sf <- st_as_sf(Biomeclimate, coords = c("x", "y"), crs = st_crs(kuchlermap))

Biome_sf <- st_join(Biome_sf,kuchlermap[,c('TYPE','CODE')])

Biome_sftrans <- st_transform(Biome_sf, st_crs(SoilpH))


xslope <- raster::extract(slope, Biome_sftrans)
Biome_extract <- cbind(Biome_sftrans, xslope)
xsand <- raster::extract(sand, Biome_sftrans)
Biome_extract <- cbind(Biome_extract, xsand)
xSoilpH <- raster::extract(SoilpH, Biome_sftrans)
Biome_extract <- cbind(Biome_extract, xSoilpH)
xhydric <- raster::extract(hydric, Biome_sftrans)
Biome_extract <- cbind(Biome_extract, xhydric)
xsalids <- raster::extract(salids, Biome_sftrans)
Biome_extract <- cbind(Biome_extract, xsalids)
xsealevel <- raster::extract(sealevel , Biome_sftrans)
Biome_extract <- cbind(Biome_extract, xsealevel)
xwater100k <- raster::extract(water100k, Biome_sftrans)
Biome_extract <- cbind(Biome_extract, xwater100k)
Biome_extract$Cindex <- pmin(Biome_extract$Tc, Biome_extract$Tclx+15)
Biome_extract$gdem <- Biome_extract$Elevation

colnames(Biome_extract)[which(names(Biome_extract) %in% c('Tg','xslope','xsand', 'xSoilpH','xhydric','xsalids','xsealevel','xwater100k'))] <- 
  c('Tgs','slope','sand','SoilpH','hydric','salids','sealevel','water100k')
Biome_extract <- subset(Biome_extract, !is.na(SoilpH)& !is.na(slope))
Biome_extract$TYPE[is.na(Biome_extract$TYPE)] <- 'unk'
Biome_extract$CODE[is.na(Biome_extract$CODE)] <- 0
#biomlist <- st_set_geometry(Biome_extract, NULL)
#biomlist <- aggregate(biomlist[,(c('BIOME'))], by=list(biomlist$BIOME, biomlist$biomname, biomlist$ECO_NAME, biomlist$CODE, biomlist$TYPE), FUN='length')
#colnames(biomlist) <- c('BIOME', 'biomname', 'ECO_NAME', 'CODE', 'TYPE', 'count')
#write.csv(biomlist, 'output/biomlist.csv')
biomlist <- read.csv('data/biomlist.csv')
Biome_extract <- merge(Biome_extract, biomlist[,c('ECO_NAME','TYPE','synbiome', 'synth')], by= c('ECO_NAME','TYPE'))
plot(Biome_extract[,'synbiome'])
Biome_extract <- subset(Biome_extract, !synbiome %in% 'Mangroves' | Elevation <= 5)
rasters<-stack(Tgs,Tc,Tclx,M,Surplus,Deficit,pAET,slope,sand,SoilpH,hydric,salids,sealevel)
     
#classification tree packages
library(randomForest)
library(rpart)
library(rpart.plot)
#Discriminate Analysis
library(MASS)           
#LDA___________________________________________________________________________
selectBiome<-Biome_extract
selectBiome$synbiome <- as.character(selectBiome$synbiome)
selectBiome[selectBiome$synbiome %in% c('Cool Wetland','Warm Wetland'),]$synbiome <- 'Wetland'
Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$synbiome), FUN=length)
names(Spwts)<-c("synbiome","x")
Spwts$myweights<- (10000/(Spwts$x/1+10000))/10
Spwts$effective<-Spwts$myweights*Spwts$x
Spwtsum<-sum(Spwts$effective)
Spwts$prior<-Spwts$effective/Spwtsum
Spwtsum<-sum(Spwts$prior)
#selectBiome2<-subset(selectBiome, (grepl('Warm Desert',selectBiome$synbiome)))
vegmodel<-qda(synbiome ~   Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, data=selectBiome, prior=Spwts$prior)



#impliment model__________________________________________________________________________
#predict(rasters,biomeclass,progress="window",overwrite=TRUE,filename="kuchlermodel.tif") 
#vegmap2<-predict(rasters,biomeclass,progress="window",filename="kuchlermodel.tif") 

vegmap<-predict(rasters,vegmodel,progress="window",overwrite=TRUE, filename="output/kuchlermodelqda.tif") 

plot(vegmap)      
               
               
#randomforest
selectBiome<-Biome_extract
Spwts<-aggregate(x=selectBiome$synbiome, by=list(selectBiome$synbiome), FUN=length)
names(Spwts)<-c("synbiome","x")
Spwts$myweights<- (10000/(Spwts$x/1+1000))/10
selectBiome<-merge(selectBiome,Spwts,by='synbiome')
rf <- randomForest(as.factor(synbiome) ~  Tgs+Tc+Tclx+M+Surplus+Deficit+pAET+slope+sand+SoilpH+hydric+salids+sealevel, data=selectBiome, classwt=selectBiome$wt, importance=TRUE, ntree=200, na.action=na.omit )
# Make plot
rf#statistical summary
varImpPlot(rf)

vegmaprf<-predict(rasters,rf,progress="window",overwrite=TRUE, filename="output/kuchlermodelrf.tif") 

plot(vegmaprf)   

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
               

               