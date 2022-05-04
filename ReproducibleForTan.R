write.csv(g3[,1:18], "golfexample.csv")

g4 <- read.csv("golfexample.csv")
mc2 <- MonoClust(as.data.frame(g4[,2:19]), nclusters = 3)
g4$Cluster <- mc2$membership %>% as.factor()

colors = c(viridis::inferno(3, begin = 0, end = .7))

#par(mar = c(2,2,2,2))
#this is where we'd put in the dendrogram
plot(mc2, main = "3-Cluster Dendrogram", cols = colors, col.type = "l")

