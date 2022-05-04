## Code to step through the Cluster analysis on the golf data

# ---------------------------------------------------------------------------------------------------#
#### Libarary packages and read in Monoclust                                                     #####
#  --------------------------------------------------------------------------------------------------# 
library(readr);library(ggplot2);library(dplyr)
library(magrittr); library(GGally)
#ClusGap example with Monothetic clustering
library(cluster)
## Source the MonoClust Function
source('~/Doctoral Work/monoClust/R/MonoClust.R')

# --------------------------------------------------------------------------------------------------#
#### Read in some Data                                                                          #####
# --------------------------------------------------------------------------------------------------#
golf <- read_csv("~/Doctoral Work/Golf Example/Round3.csv")
#table(golf$back_nine)
golf2 <- dplyr::filter(golf, back_nine == FALSE)
dim(golf2)

golfmat <- as.data.frame(golf[,1:18])
ggparcoord(golf,columns = 1:18, groupColumn = 25, splineFactor = 0, alpha = .6) + ggtitle("US Open Round 3 Trajectories (Smoothed)") + theme_bw() + scale_color_viridis_d(option = "D", begin = 0.2, end = 0.9)

golfmat2 <- as.data.frame(golf2[,1:18])
#mc1<- MonoClust(golfmat, nclusters = 2)

### Need function to interface with the Gap Statisic calculator 
## argument 1 is the dataset, argument 2 is K
## output is a list
MC1 <- function(x, k){
  xnew <- as.data.frame(x)
  temp <- MonoClust(xnew, nclusters = k)
  
return(list(cluster = temp$Membership))}

#Try it on the ruspini dataset
gs_mc_RU <- clusGap(ruspini, FUN = MC1, K.max = 8, B = 60)
gs_mc_RU
plot(gs_mc_RU, main = "Gap statistic for the 'ruspini' data")
mtext("k = 4 is best .. and  k = 5  pretty close")

## ggplot method for this: (from https://joey711.github.io/phyloseq/gap-statistic.html)
plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
  library(ggplot2)
  gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
  p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
  p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim), color = "tomato2", alpha = .6)
  p = p + ggtitle(title) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

plot_clusgap(gs_mc_RU, title = "Gap Statistic on Ruspini (Monothetic)")



## Now we can try this on the Golf data
gs_mc_golf <- clusGap(golfmat, FUN = MC1, K.max = 8, B = 60, d.power = 2, spaceH0 = "original")
gs_mc_golf
plot_clusgap(gs_mc_golf) #this is a problem



## try with removal of 2, 4, 7, 11, 13, 17

gs_mc_golf2 <- clusGap(golfmat2[,-c(2,4,7,11,13,17)], FUN = MC1, K.max = 8, B = 60)
gs_mc_golf2
plot_clusgap(gs_mc_golf2)


#### Other clustering methods: ####

## partitioning around medoids
gs_mc_golf3 <- clusGap(golfmat2, FUN = pam, K.max = 8, B = 60)
gs_mc_golf3
plot_clusgap(gs_mc_golf3)

## Hierarchical clustering (standard, with squared Ward's distance)
hclusCut <- function(x, k, d.meth = "euclidean", ...)
  list(cluster = cutree(hclust(dist(x, method=d.meth), method = "ward.D2", ...), k=k))
gs_mc_golf2 <- clusGap(golfmat2, FUN = hclusCut, K.max = 20, B = 100, spaceH0 = "original", d.power = 2)
gs_mc_golf2
plot_clusgap(gs_mc_golf2)



#### Alternative Methods ####

#reads in all 4 days
g1 <- read_csv("~/Doctoral Work/Golf Example/Round1.csv")
g2 <- read_csv("~/Doctoral Work/Golf Example/Round2.csv")
g3 <- read_csv("~/Doctoral Work/Golf Example/Round3.csv")
g4 <- read_csv("~/Doctoral Work/Golf Example/Round4.csv")
g1$Day <- rep("Day1", nrow(g1))
g2$Day <- rep("Day2", nrow(g2))
g3$Day <- rep("Day3",nrow(g3))
g4$Day <- rep("Day4", nrow(g4))

golf_full <- rbind(g1,g2,g3,g4)
golf_full_mat <- golf_full[,1:18] %>% as.data.frame()

#Look at the gap statistic (try with increasing B, use MC and Wards, go either squared or not, use original)
gs_mc_golf4 <- clusGap(golf_full_mat, FUN = hclusCut, K.max = 10, B = 60)
gs_mc_golf4
plot_clusgap(gs_mc_golf4)



## Silhouette statistic
data(ruspini)
pr4 <- pam(golf_full_mat, 4)
str(si <- silhouette(pr4))
(ssi <- summary(si))
plot(si) # silhouette plot
plot(si, col = c("red", "green", "blue", "purple"))# with cluster-wise coloring






