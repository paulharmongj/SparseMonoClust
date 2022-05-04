## Assess Weights from Sparcl and Distance Matrix
library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(sparcl)
#test dataset
nba <- read.csv("C:/Users/paulh/OneDrive/Documents/Utah Jazz Simulations/UtahJazzPredictions/Player_Stats.csv")
head(nba)


## Consider a dataframe
nba1 <- select(nba, FG, ORB, PTS) %>% apply(2,scale)

hc1 <- HierarchicalSparseCluster(as.matrix(nba1), wbound = 1.01)
hc1$ws
hc1$dists #gives the feature-wise dissimilarity matrix (nxn)*p
hc1$u  %>% round(3)#this is the 'sparse' distance matrix

#heatmaps
heatmap(hc1$u, main = "Heatmap of Sparse Distance Matrix")
image(hc1$u, main = "Heatmap of Sparse Distance Matrix") ## Target

##############################
#### Try with Multiplication #
##############################
nbamat <- as.matrix(nba1)
nba_new <- tibble(x1 = hc1$ws[1] * nbamat[,1], x2 = hc1$ws[2] * nbamat[,2], x3 = hc1$ws[3] * nbamat[,3])

d1 <- dist(nba_new) %>% as.matrix()
image(d1, main = "Start with Weights")

#### Verify that we can calculate u from the original weights and data
#sum_j wj sum_ii' d_i,i' 





#### Ignore ####
# instead, try this
# take your n = 15 observations and generate n*n)xp = 225*3 distance matrix

a = dist(scale(nbamat[,1])) #15x15
b = dist(scale(nbamat[,2]))
c = dist(scale(nbamat[,3]))

# sum over j the weights within each feature
mat2 <- as.matrix(sqrt(hc1$ws[1]) * a + sqrt(hc1$ws[2]) * b + sqrt(hc1$ws[3]) * c)
image(mat2, main = "Summed version")
hc2 <- hclust(as.dist(mat2), method = "complete")

## vectorize and plot sets of points vs eachother

## temp calc

temp <- matrix(nrow = 3, c(1,2,0,3,2,1,0,1,2))
dist(temp) %>% as.matrix()

hc2 <- HierarchicalSparseCluster(temp, wbound = 1.2)
hc2$ws
hc2$u










