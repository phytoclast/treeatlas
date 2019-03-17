library(sf)
library(raster)
library(rgdal)
library(velox)
library(lattice)
library(dplyr)
library(ggplot2)
treelist <- read.delim('treefiles.txt')
forms <- read.delim('data/NorthAmericaGrowthHabits20180924.txt', encoding="Latin1")
syns <- read.delim('data/BONAPGRIN_NAMES.txt', encoding="UTF-8")
colnames(syns) <- c("Binomial","Taxon","B_Taxon","B_Binomial","G_Taxon","G_Binomial") 
formsyns <- merge(syns, forms, by.x='Binomial')
formsyns <- unique(formsyns[!is.na(formsyns$B_Binomial) & !is.na(formsyns$HabitSymbol)& !formsyns$B_Binomial %in% "", c('HabitSymbol','B_Binomial')])
treelistsyns <- merge(treelist, syns, by.x='Latin.Name', by.y='Binomial', all.x=TRUE)
treelistsyns <- unique(treelistsyns[!is.na(treelistsyns$B_Binomial),c('name','B_Binomial')])
treelistforms <- merge(treelistsyns, formsyns, by='B_Binomial', all.x=TRUE)
BE <- as.character(unique(treelistforms[grepl('BE', treelistforms$HabitSymbol),'name']))
BD <- as.character(unique(treelistforms[grepl('BD', treelistforms$HabitSymbol),'name']))
NE <- as.character(unique(treelistforms[grepl('NE', treelistforms$HabitSymbol)|grepl('ND', treelistforms$HabitSymbol),'name']))
unique(treelistforms$HabitSymbol)

Tclx <- raster('nam/Tclx.tif')
Tgs <- raster('nam/Tgs.tif')
Tc <- raster('nam/Tc.tif')
MAP <- raster('nam/MAP.tif')
Deficit <- raster('nam/Deficit.tif')
Surplus <- raster('nam/Surplus.tif')
pAET <- raster('nam/pAET.tif')
slope <- raster('nam/slope.tif')
sand <- raster('nam/sand.tif')
SoilpH <- raster('nam/SoilpH.tif')
hydric <- raster('nam/hydric.tif')
salids <- raster('nam/salids.tif')

M <- MAP/(Deficit + MAP - Surplus)
vTclx <- velox(Tclx)
vTgs <- velox(Tgs)
vTc <- velox(Tc)
vM <- velox(M)
vDeficit <- velox(Deficit)
vSurplus <- velox(Surplus)
vpAET <- velox(pAET)
vslope <- velox(slope)
vsand <- velox(sand)
vSoilpH<- velox(SoilpH)
vhydric<- velox(hydric)
vsalids<- velox(salids)
line <- cbind(file = 'name',taxon = 'name',name = 'name', Tclx_max=0, Tclx_min=0, Tgs_max=0, Tgs_min=0)

top = 0.98
bot = 0.02



n <- rownames(treelist[grepl('liqustyr', treelist$name),])

spp<-read_sf(treelist[n,1])
spp2 <- 	st_transform(spp, st_crs(Tclx))
spp2 <- spp2[spp2$CODE %in% '1',]

#plot(Tclx, main=treelist[n,4])
#plot(st_geometry(spp2), add=T)


#Tclx----

e <- vTclx$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tclx', ylab='Density', main='')
Tclx_max <- quantile(e[,2], top)
Tclx_min <- quantile(e[,2], bot)

#Tgs----

e <- vTgs$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tgs', ylab='Density', main='')
Tgs_max <- quantile(e[,2], top)
Tgs_min <- quantile(e[,2], bot)

#Tc----

e <- vTc$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tc', ylab='Density', main='')
Tc_max <- quantile(e[,2], top)
Tc_min <- quantile(e[,2], bot)

#M----

e <- vM$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='M', ylab='Density', main='')
M_max <- quantile(e[,2], top)
M_min <- quantile(e[,2], bot)

#Surplus----

e <- vSurplus$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Surplus', ylab='Density', main='')
Surplus_max <- quantile(e[,2], top)
Surplus_min <- quantile(e[,2], bot)

#Deficit----

e <- vDeficit$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Deficit', ylab='Density', main='')
Deficit_max <- quantile(e[,2], top)
Deficit_min <- quantile(e[,2], bot)

#pAET----

e <- vpAET$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

pAET_max <- quantile(e[,2], top)
pAET_min <- quantile(e[,2], bot)


#slope----

e <- vslope$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

slope_max <- quantile(e[,2], top)
slope_min <- quantile(e[,2], bot)

#sand----

e <- vsand$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

sand_max <- quantile(e[,2], top)
sand_min <- quantile(e[,2], bot)

#SoilpH----

e <- vSoilpH$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

SoilpH_max <- quantile(e[,2], top)
SoilpH_min <- quantile(e[,2], bot)

#hydric----

e <- vhydric$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

hydric_max <- quantile(e[,2], top)
hydric_min <- quantile(e[,2], bot)

#salids----

e <- vsalids$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='pAET', ylab='Density', main='')

salids_max <- quantile(e[,2], top)
salids_min <- quantile(e[,2], bot)


envelopT <- (Tclx <= Tclx_max)*(Tclx >= Tclx_min)*
  (Tgs <= Tgs_max)*(Tgs >= Tgs_min)*
  (Tc <= Tc_max)*(Tc >= Tc_min)

envelopM <- (M <= M_max)*(M >= M_min)*
  (pAET <= pAET_max)*(pAET >= pAET_min)*
  (Deficit <= Deficit_max)*(Deficit >= Deficit_min)*
  (Surplus <= Surplus_max)*(Surplus >= Surplus_min)

envelopS <- (slope <= slope_max)*(slope >= slope_min)*
  (sand <= sand_max)*(sand >= sand_min)*
  (SoilpH <= SoilpH_max)*(SoilpH >= SoilpH_min)*
  (hydric <= hydric_max)*(hydric >= hydric_min)*
  (salids <= salids_max)*(salids >= salids_min)
envelop3 <- envelopT*envelopM*envelopS



minusTclx <- (slope <= slope_max)*(slope >= slope_min)*
  (sand <= sand_max)*(sand >= sand_min)*
  (SoilpH <= SoilpH_max)*(SoilpH >= SoilpH_min)*
  (hydric <= hydric_max)*(hydric >= hydric_min)*
  (salids <= salids_max)*(salids >= salids_min)*
  (Tclx <= Tclx_max)*(Tclx >= Tclx_min)*
  (Tgs <= Tgs_max)*(Tgs >= Tgs_min)*
  (Tc <= Tc_max)*(Tc >= Tc_min)*
  (M <= M_max)*(M >= M_min)*
  (pAET <= pAET_max)*(pAET >= pAET_min)*
  (Deficit <= Deficit_max)*(Deficit >= Deficit_min)*
  (Surplus <= Surplus_max)*(Surplus >= Surplus_min)

diff <- minusTclx-envelop3
cellStats(diff, stat='mean')
plot(diff, main=treelist[n,3])

plot(envelop3, main=treelist[n,3])
plot(st_geometry(spp2), add=T)

