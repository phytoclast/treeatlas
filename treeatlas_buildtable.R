library(sf)
library(raster)
library(rgdal)
library(velox)
library(lattice)
library(dplyr)
library(ggplot2)
treelist <- read.delim('treefiles.txt')
bt_frost <- read.delim('bt_frost2.txt')
kuchlermap <-read_sf('data/kuchler.shp')
forms <- read.delim('data/NorthAmericaGrowthHabits20180924.txt', encoding="Latin1")
syns <- read.delim('data/BONAPGRIN_NAMES.txt', encoding="UTF-8")
colnames(syns) <- c("Binomial","Taxon","B_Taxon","B_Binomial","G_Taxon","G_Binomial") 
formsyns <- merge(syns, forms, by.x='Taxon', by.y='Scientific.Name')
formsyns <- unique(formsyns[!is.na(formsyns$B_Binomial) & !is.na(formsyns$HabitSymbol)& !formsyns$B_Binomial %in% "", c('HabitSymbol','B_Binomial')])
treelistsyns <- merge(treelist, syns, by.x='Latin.Name', by.y='Taxon', all.x=TRUE)
treelistsyns <- unique(treelistsyns[!is.na(treelistsyns$B_Binomial),c('name','B_Binomial')])
treelistforms <- merge(treelistsyns, formsyns, by='B_Binomial', all.x=TRUE)
BE <- as.character(unique(treelistforms[grepl('BE', treelistforms$HabitSymbol),'name']))
PALM <- as.character(unique(treelistforms[grepl('P', treelistforms$HabitSymbol),'name']))
BD <- as.character(unique(treelistforms[grepl('BD', treelistforms$HabitSymbol),'name']))
NE <- as.character(unique(treelistforms[grepl('NE', treelistforms$HabitSymbol)|grepl('ND', treelistforms$HabitSymbol),'name']))
SU <- as.character(unique(treelistforms[grepl('U', treelistforms$HabitSymbol),'name']))

unique(treelistforms$HabitSymbol)
bigTgs <- raster('C:/a/Ecological_Sites/GIS/model/global/Tgs.tif')
bigTclx <- raster('C:/a/Ecological_Sites/GIS/model/global/Tclx.tif')
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
Cindex <- min(Tclx+15, Tc)

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
vCindex<- velox(Cindex)

vbigTgs <- velox(bigTgs)
vbigTclx <- velox(bigTclx)

line <- cbind(file = 'name',taxon = 'name',name = 'name', Tclx_max=0, Tclx_min=0, Tgs_max=0, Tgs_min=0, Cindex_max = 0, Cindex_min=0)

top = 0.98
bot = 0.02

for(i in 1:nrow(treelist)){

#n <- rownames(treelist[grepl('larrdiva', treelist$name),])
n <- i
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

#Cindex----

e <- vCindex$extract(spp2, df=T, small = T)
e <- e[!is.na(e[,2]),]
#densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tgs', ylab='Density', main='')
Cindex_max <- quantile(e[,2], top)
Cindex_min <- quantile(e[,2], bot)

line2 <- cbind(file = as.character(treelist[n,2]), taxon = as.character(treelist[n,3]),
               name = as.character(treelist[n,4]), Tclx_max, Tclx_min, Tgs_max, Tgs_min, Cindex_max, Cindex_min)
line <- rbind(line, line2)

}

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
hist(envelopM)
envelop3 <- envelopT*envelopM*envelopS
plot(envelop3)
plot(st_geometry(spp2), add=T)

TreeClim <- as.data.frame(line)
TreeClim$Tclx_min <- as.numeric(as.character(TreeClim$Tclx_min))
TreeClim$Tclx_max <- as.numeric(as.character(TreeClim$Tclx_max))
TreeClim$Tgs_min <- as.numeric(as.character(TreeClim$Tgs_min))
TreeClim$Tgs_max <- as.numeric(as.character(TreeClim$Tgs_max))
TreeClim$Cindex_max <- as.numeric(as.character(TreeClim$Cindex_max))
TreeClim$Cindex_min <- as.numeric(as.character(TreeClim$Cindex_min))
TreeClim <- unique(TreeClim[2:nrow(TreeClim),])
rownames(TreeClim) <- TreeClim$file
#saveRDS(TreeClim, 'data/TreeClim.RDS')
#TreeClim <- readRDS('data/TreeClim.RDS')

BroadEvergreens <- unique(TreeClim[TreeClim$file %in% BE,])
Monocot <- unique(TreeClim[TreeClim$file %in% PALM,])
Succulent <- unique(TreeClim[TreeClim$file %in% SU,])
BroadDeciduous <- unique(TreeClim[TreeClim$file %in% BD,])
Needle <- unique(TreeClim[TreeClim$file %in% NE,])
TreeClim$file <- as.character(TreeClim$file)
densityplot(Needle$Tclx_min)
densityplot(BroadDeciduous$Tclx_min, add=T)
densityplot(BroadEvergreens$Tclx_min, add=T)
TreeClim$diff <- TreeClim$Cindex_min - TreeClim$Tclx_min
ggplot() + 
  geom_density(data = BroadEvergreens, 
               mapping = aes(x = Cindex_min, fill='Broadleaf Evergreen', color='Broadleaf Evergreen'),size=1, alpha=0.1)+
  geom_density(data = BroadDeciduous, 
               mapping = aes(x = Cindex_min, fill='Broadleaf Deciduous', color='Broadleaf Deciduous'), size=1, alpha=0.1)+
  geom_density(data = Monocot, 
               mapping = aes(x = Cindex_min, fill='Monocot', color='Monocot'), size=1, alpha=0.1)+
  geom_density(data = Succulent, 
               mapping = aes(x = Cindex_min, fill='Succulent', color='Succulent'), size=1, alpha=0.1)+
  geom_density(data = Needle, 
               mapping = aes(x = Cindex_min, fill='Needleleaf', color='Needleleaf'), size=1, alpha=0.1)+
  scale_fill_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                         "Broadleaf Deciduous" = "green",
                                         "Monocot" = "yellow",
                                         "Succulent" = "purple",
                                         "Needleleaf" = "blue"))+
  scale_color_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                          "Broadleaf Deciduous" = "dark green",
                                          "Monocot" = "orange",
                                          "Succulent" = "purple",
                                          "Needleleaf" = "blue"))+
  scale_x_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65+15, 5+15, by = 5))+
  scale_y_continuous(name= "Tree Species Density")

ggplot() + 
  geom_smooth(data = BroadEvergreens, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Broadleaf Evergreen', color='Broadleaf Evergreen'),size=1, alpha=0.1)+
  geom_smooth(data = BroadDeciduous, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Broadleaf Deciduous', color='Broadleaf Deciduous'), size=1, alpha=0.1)+
  geom_smooth(data = Monocot, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Monocot', color='Monocot'), size=1, alpha=0.1)+
  geom_smooth(data = Succulent, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Succulent', color='Succulent'), size=1, alpha=0.1)+
  geom_smooth(data = Needle, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Needleleaf', color='Needleleaf'), size=1, alpha=0.1)+
  geom_point(data = BroadEvergreens, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Broadleaf Evergreen', color='Broadleaf Evergreen'),size=1, alpha=0.1)+
  geom_point(data = BroadDeciduous, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Broadleaf Deciduous', color='Broadleaf Deciduous'), size=1, alpha=0.1)+
  geom_point(data = Monocot, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Monocot', color='Monocot'), size=1, alpha=0.1)+
  geom_point(data = Succulent, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Succulent', color='Succulent'), size=1, alpha=0.1)+
  geom_point(data = Needle, 
               mapping = aes(x = Cindex_min, y = Tgs_min, fill='Needleleaf', color='Needleleaf'), size=1, alpha=0.1)+
  scale_fill_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                         "Broadleaf Deciduous" = "green",
                                         "Monocot" = "yellow",
                                         "Succulent" = "purple",
                                         "Needleleaf" = "blue"))+
  scale_color_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                          "Broadleaf Deciduous" = "dark green",
                                          "Monocot" = "orange",
                                          "Succulent" = "purple",
                                          "Needleleaf" = "blue"))+
  scale_x_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65+15, 5+15, by = 5))+
  scale_y_continuous(name= "Maximum Growing Season Temperature", breaks=seq(0, 36, by = 2))

ggplot() + 
  geom_density(data = BroadEvergreens, 
               mapping = aes(x = Tgs_min, fill='Broadleaf Evergreen', color='Broadleaf Evergreen'),size=1, alpha=0.1)+
  geom_density(data = BroadDeciduous, 
               mapping = aes(x = Tgs_min, fill='Broadleaf Deciduous', color='Broadleaf Deciduous'), size=1, alpha=0.1)+
  geom_density(data = Monocot, 
               mapping = aes(x = Tgs_min, fill='Monocot', color='Monocot'), size=1, alpha=0.1)+
  geom_density(data = Succulent, 
               mapping = aes(x = Tgs_min, fill='Succulent', color='Succulent'), size=1, alpha=0.1)+
  geom_density(data = Needle, 
               mapping = aes(x = Tgs_min, fill='Needleleaf', color='Needleleaf'), size=1, alpha=0.1)+
  scale_fill_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                         "Broadleaf Deciduous" = "green",
                                         "Monocot" = "yellow",
                                         "Succulent" = "purple",
                                         "Needleleaf" = "blue"))+
  scale_color_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                          "Broadleaf Deciduous" = "dark green",
                                          "Monocot" = "orange",
                                          "Succulent" = "purple",
                                          "Needleleaf" = "blue"))+
  scale_x_continuous(name= "Minimum Growing Season Temperature", breaks=seq(0, 36, by = 2))+
  scale_y_continuous(name= "Tree Species Density")

ggplot() + 
  geom_density(data = BroadEvergreens, 
               mapping = aes(x = Tgs_max, fill='Broadleaf Evergreen', color='Broadleaf Evergreen'),size=1, alpha=0.1)+
  geom_density(data = BroadDeciduous, 
               mapping = aes(x = Tgs_max, fill='Broadleaf Deciduous', color='Broadleaf Deciduous'), size=1, alpha=0.1)+
  geom_density(data = Monocot, 
               mapping = aes(x = Tgs_max, fill='Monocot', color='Monocot'), size=1, alpha=0.1)+
  geom_density(data = Succulent, 
               mapping = aes(x = Tgs_max, fill='Succulent', color='Succulent'), size=1, alpha=0.1)+
  geom_density(data = Needle, 
               mapping = aes(x = Tgs_max, fill='Needleleaf', color='Needleleaf'), size=1, alpha=0.1)+
  scale_fill_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                         "Broadleaf Deciduous" = "green",
                                         "Monocot" = "yellow",
                                         "Succulent" = "purple",
                                         "Needleleaf" = "blue"))+
  scale_color_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                          "Broadleaf Deciduous" = "dark green",
                                          "Monocot" = "orange",
                                          "Succulent" = "purple",
                                          "Needleleaf" = "blue"))+
  scale_x_continuous(name= "Maximum Growing Season Temperature", breaks=seq(0, 36, by = 2))+
  scale_y_continuous(name= "Tree Species Density")



ggplot() + 
  geom_density2d(data = BroadEvergreens, 
               mapping = aes(x = Tgs_min, y = Tclx_min, fill='Broadleaf Evergreen', color='Broadleaf Evergreen'),size=1, alpha=0.8)+
  geom_density2d(data = BroadDeciduous, 
               mapping = aes(x = Tgs_min, y = Tclx_min, fill='Broadleaf Deciduous', color='Broadleaf Deciduous'), size=1, alpha=0.8)+
  geom_density2d(data = Monocot, 
               mapping = aes(x = Tgs_min, y = Tclx_min, fill='Monocot', color='Monocot'), size=1, alpha=0.8)+
  geom_density2d(data = Succulent, 
               mapping = aes(x = Tgs_min, y = Tclx_min, fill='Succulent', color='Succulent'), size=1, alpha=0.8)+
  geom_density2d(data = Needle, 
               mapping = aes(x = Tgs_min, y = Tclx_min, fill='Needleleaf', color='Needleleaf'), size=1, alpha=0.8)+
  scale_fill_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                         "Broadleaf Deciduous" = "green",
                                         "Monocot" = "yellow",
                                         "Succulent" = "purple",
                                         "Needleleaf" = "blue"))+
  scale_color_manual("Legend", values = c("Broadleaf Evergreen" = "red",
                                          "Broadleaf Deciduous" = "dark green",
                                          "Monocot" = "orange",
                                          "Succulent" = "purple",
                                          "Needleleaf" = "blue"))+
  scale_x_continuous(name= "Minimum Growing Season Temperature", breaks=seq(0, 36, by = 2))+
  scale_y_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65, 5, by = 5))

library(BiodiversityR)
library(ape)
library(cluster)

dist <- vegdist(TreeClim[,4:7], method = 'gower')
clust <- hclust(dist)
plot(clust, cex = 0.1)
cuttree <- as.data.frame(cutree(clust,8))
cuttree$file <- rownames(cuttree)
colnames(cuttree) <- c('cluster', 'file')
TreeClimcut <- merge(TreeClim, cuttree, by='file')
TreeClimcut$cluster <- as.character(TreeClimcut$cluster)
ggplot() + 
  geom_density2d(data = TreeClimcut, 
                 mapping = aes(x = Tgs_min, y = Tclx_min, fill=cluster, color=cluster),size=1, alpha=0.8)+
  geom_point(data = TreeClimcut, 
                 mapping = aes(x = Tgs_min, y = Tclx_min, fill=cluster, color=cluster),size=1, alpha=0.8)+
  scale_fill_manual("Legend", values = c("1" = "red",
                                         "2" = "green",
                                         "3" = "blue",
                                         "4" = "orange",
                                         "5" = "yellow",
                                         "6" = "cyan",
                                         "7" = "purple",
                                         "8" = "pink"))+
  scale_color_manual("Legend", values = c("1" = "red",
                                          "2" = "green",
                                          "3" = "blue",
                                          "4" = "orange",
                                          "5" = "yellow",
                                          "6" = "cyan",
                                          "7" = "purple",
                                          "8" = "pink"))+
  scale_x_continuous(name= "Minimum Growing Season Temperature", breaks=seq(0, 36, by = 2))+
  scale_y_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65, 5, by = 5))


ggplot() + 
  geom_boxplot(data = TreeClimcut, 
               mapping = aes(x = cluster, y = Tgs_min, fill=cluster, color=cluster),size=1, alpha=0.5)+
  scale_fill_manual("Legend", values = c("1" = "red",
                                         "2" = "green",
                                         "3" = "blue",
                                         "4" = "orange",
                                         "5" = "yellow",
                                         "6" = "cyan",
                                         "7" = "purple",
                                         "8" = "pink"))+
  scale_color_manual("Legend", values = c("1" = "red",
                                          "2" = "green",
                                          "3" = "blue",
                                          "4" = "orange",
                                          "5" = "yellow",
                                          "6" = "cyan",
                                          "7" = "purple",
                                          "8" = "pink"))+
  scale_y_continuous(name= "Minimum Growing Season Temperature", breaks=seq(0, 36, by = 2))
#scale_y_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65, 5, by = 5))
ggplot() + 
  geom_boxplot(data = TreeClimcut, 
               mapping = aes(x = cluster, y = Tgs_max, fill=cluster, color=cluster),size=1, alpha=0.5)+
  geom_boxplot(data = TreeClimcut, 
               mapping = aes(x = cluster, y = Tgs_min, fill=cluster, color=cluster),size=1, alpha=0.5)+
  scale_fill_manual("Legend", values = c("1" = "red",
                                         "2" = "green",
                                         "3" = "blue",
                                         "4" = "orange",
                                         "5" = "yellow",
                                         "6" = "cyan",
                                         "7" = "purple",
                                         "8" = "pink"))+
  scale_color_manual("Legend", values = c("1" = "red",
                                          "2" = "green",
                                          "3" = "blue",
                                          "4" = "orange",
                                          "5" = "yellow",
                                          "6" = "cyan",
                                          "7" = "purple",
                                          "8" = "pink"))+
  scale_y_continuous(name= "Maximum Growing Season Temperature", breaks=seq(0, 36, by = 2))
#scale_y_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65, 5, by = 5))

ggplot() + 
  geom_boxplot(data = TreeClimcut, 
               mapping = aes(x = cluster, y = Tclx_min, fill=cluster, color=cluster),size=1, alpha=0.5)+
  scale_fill_manual("Legend", values = c("1" = "red",
                                         "2" = "green",
                                         "3" = "blue",
                                         "4" = "orange",
                                         "5" = "yellow",
                                         "6" = "cyan",
                                         "7" = "purple",
                                         "8" = "pink"))+
  scale_color_manual("Legend", values = c("1" = "red",
                                          "2" = "green",
                                          "3" = "blue",
                                          "4" = "orange",
                                          "5" = "yellow",
                                          "6" = "cyan",
                                          "7" = "purple",
                                          "8" = "pink"))+
  #scale_y_continuous(name= "Minimum Growing Season Temperature", breaks=seq(0, 36, by = 2))
scale_y_continuous(name= "Minnimum Annual Extreme Low", breaks=seq(-65, 5, by = 5))
aggregate(TreeClimcut[,c('Tclx_max','Tclx_min','Tgs_max','Tgs_min')], by=list(TreeClimcut$cluster), FUN='median')
aggregate(TreeClimcut[,c('Tclx_max','Tclx_min','Tgs_max','Tgs_min')], by=list(TreeClimcut$cluster), FUN='sd')


library(rpart)
library(rpart.plot)
model <- rpart(cluster~ Tclx_min +Tclx_max +Tgs_min +Tgs_max, data= TreeClimcut, method = "class", maxdepth = 4)
rpart.plot(model, extra=108)

bt_frost$b01 <- 0
bt_frost$b02 <- 0
bt_frost$b03 <- 0
bt_frost$b04 <- 0
bt_frost$b05 <- 0
bt_frost$b06 <- 0
bt_frost$b07 <- 0
bt_frost$b08 <- 0
bt_frost$b09 <- 0
bt_frost$b10 <- 0
bt_frost$b11 <- 0
bt_frost$b12 <- 0
for (i in 0:11){
  bt_frost[,which(colnames(bt_frost)=='b01')+i]  <- bt_frost[,which(colnames(bt_frost)=='t01')+i]*
    (bt_frost[,which(colnames(bt_frost)=='t01')+i]>0)*1
}

bt_frost$bl01 <- 0
bt_frost$bl02 <- 0
bt_frost$bl03 <- 0
bt_frost$bl04 <- 0
bt_frost$bl05 <- 0
bt_frost$bl06 <- 0
bt_frost$bl07 <- 0
bt_frost$bl08 <- 0
bt_frost$bl09 <- 0
bt_frost$bl10 <- 0
bt_frost$bl11 <- 0
bt_frost$bl12 <- 0
for (i in 0:11){
  bt_frost[,which(colnames(bt_frost)=='bl01')+i]  <- bt_frost[,which(colnames(bt_frost)=='tl01')+i]*
    (bt_frost[,which(colnames(bt_frost)=='tl01')+i]>0)*1
}
bt_frost$gl01 <- 0
bt_frost$gl02 <- 0
bt_frost$gl03 <- 0
bt_frost$gl04 <- 0
bt_frost$gl05 <- 0
bt_frost$gl06 <- 0
bt_frost$gl07 <- 0
bt_frost$gl08 <- 0
bt_frost$gl09 <- 0
bt_frost$gl10 <- 0
bt_frost$gl11 <- 0
bt_frost$gl12 <- 0
for (i in 0:11){
  bt_frost[,which(colnames(bt_frost)=='gl01')+i]  <-( bt_frost[,which(colnames(bt_frost)=='tl01')+i]-10)*
    (bt_frost[,which(colnames(bt_frost)=='tl01')+i]>10)*1
}
colnames(bt_frost)
bt_frost$bt <- apply(bt_frost[,c('b01', 'b02', 'b04', 'b04', 'b05', 'b06', 'b07', 'b08', 'b09', 'b10', 'b11', 'b12')], 1, FUN = mean)
bt_frost$btl <- apply(bt_frost[,c('bl01', 'bl02', 'bl04', 'bl04', 'bl05', 'bl06', 'bl07', 'bl08', 'bl09', 'bl10', 'bl11', 'bl12')], 1, FUN = mean)
bt_frost$gtl <- apply(bt_frost[,c('gl01', 'gl02', 'gl04', 'gl04', 'gl05', 'gl06', 'gl07', 'gl08', 'gl09', 'gl10', 'gl11', 'gl12')], 1, FUN = mean)
bt_frost$gdif <- bt_frost$btl - bt_frost$gtl 
cor(bt_frost[,c('Frost50', 'bt0510','tmin','t01','bt','btl','gtl','gdif')])
model <- glm(Frost50  ~ bt0510+bt+gdif+Elevation, data=bt_frost)
summary(model)
bt_frost$fit <- predict(model, bt_frost)
model2 <- glm(Frost50  ~ bt0510+tmin+t01+bt+btl+gdif, data=bt_frost)
summary(model2)
bt_frost$fit2 <- predict(model2, bt_frost)
bt_frostx <- subset(bt_frost, select= c("Station_Name","State","Elevation","bt0510","tmin","Frost50","Freeze50","t01","fit","fit2"))



#kuchler----
top = 0.98
bot = 0.02

kuchlertab <- unique(subset(as.data.frame(kuchlermap[,c('TYPE', 'CODE')]), select= -c(geometry)))
rownames(kuchlertab) <- kuchlertab$CODE
for(i in 1:nrow(kuchlertab)){
  #i <-4
  #n <- kuchlertab[grepl('Southeastn spruce-fir forest', kuchlertab$TYPE),'CODE']
  n <- kuchlertab[i,'CODE']
  veg <- kuchlermap[kuchlermap$CODE %in% n,]
  plot(veg)
  veg2 <- 	st_transform(veg, st_crs(Tgs))
  
  
  #plot(Tclx, main=treelist[n,4])
  #plot(st_geometry(spp2), add=T)
  
  
  #Tclx----
  
  e <- vTclx$extract(veg2, df=T, small = T)
  e <- e[!is.na(e[,2]),]
  #densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tclx', ylab='Density', main='')
  Tclx_max <- quantile(e[,2], top)
  Tclx_min <- quantile(e[,2], bot)
  
  #Tgs----
  
  e <- vTgs$extract(veg2, df=T, small = T)
  e <- e[!is.na(e[,2]),]
  #densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tgs', ylab='Density', main='')
  Tgs_max <- quantile(e[,2], 1)
  Tgs_min <- quantile(e[,2], 0)
  
  #Cindex----
  
  e <- vCindex$extract(spp2, df=T, small = T)
  e <- e[!is.na(e[,2]),]
  #densityplot(e[,2], plot.points=FALSE, bw=2, lwd=2, col='RoyalBlue', xlab='Tgs', ylab='Density', main='')
  Cindex_max <- quantile(e[,2], top)
  Cindex_min <- quantile(e[,2], bot)
  
  line2 <- cbind(file = as.character(treelist[n,2]), taxon = as.character(treelist[n,3]),
                 name = as.character(treelist[n,4]), Tclx_max, Tclx_min, Tgs_max, Tgs_min, Cindex_max, Cindex_min)
  line <- rbind(line, line2)
  
}

