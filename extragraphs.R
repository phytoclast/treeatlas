#https://rkabacoff.github.io/datavis/
# create a biplot
# load data
KuchlerClim <- readRDS('output/KuchlerClim2.RDS')
rownames(KuchlerClim) <- paste(KuchlerClim$CODE,KuchlerClim$TYPE)
# fit a principal components model
fit <- prcomp(x = KuchlerClim[,3:26], 
              center = TRUE, 
              scale = TRUE)

# plot the results
library(factoextra)
fviz_pca(fit, 
         repel = TRUE, 
         labelsize = 3) + 
  theme_bw() +
  labs(title = "Biplot of Kuchler data")

# create a heatmap


superheat(KuchlerClim[,3:26],
          scale = TRUE,
          left.label.text.size=3,
          bottom.label.text.size=3,
          bottom.label.size = .05,
          row.dendrogram = TRUE )