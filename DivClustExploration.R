### Divclust Exploration


install.packages("remotes")
remotes::install_github("chavent/divclust")
library(divclust)

## Vignette
data(protein) # pure quantitatives data
tree <- divclust(protein) # full clustering
plot(tree)
plot(1:(tree$kmax-1),tree$height,xlab="number of cluster",ylab="height",main="Split levels")
c_5 <- divclust(protein, K=3) # stops clustering to 5 clusters
plot(c_5,nqbin=4)
c_5$B*100 #explained inertia
c_5$clusters  # retrieve the list of observations in each cluster
c_5$description # and their monothetic description


par(mfrow = c(1,2))
c_3 <- divclust(protein, K=3) # stops clustering to 5 clusters
plot(c_3, nqbin = 1)

mc_tree <- MonoClust(protein, nclusters = 3)
plot(mc_tree)

#### 
dc1 <- divclust(final_dat, K = 3)

sparsedc1 <- SparseMonoClust(final_dat, nclusters = 3, wbound = 3, divclust = TRUE)
sparsedc1$clustob$clusters
sparsedc1$clustob$description

