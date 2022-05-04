#### Simulation in Line with Witten et al. 
source("C:/Users/paulh/Documents/Doctoral Work/Sparcl Witten/sparcl/sparclproject/SparseClusteringFunctions.R")
#simulate an example with 6 equally sized classes of n = 120 
# p = 2000 features, and 200 features differ between classes


# Compare Sparse Monothetic Clustering to Regular Monothetic and Hierarchical Clustering (Wards Method)
total_noise_features <- 180  #1800 in high-dim setting
num_noise_features <- 180
tsf <- 20
#### SIMULATE DATA #####################################################################
#No structure: Simluates 20 rv's with mean 50 and variance 25 (independent)
no_structure <- rmvnorm(n = 120, mean = rep(46,total_noise_features), sigma = diag(30,total_noise_features)) %>% as_tibble()

###Structured Data (simuklated with some structure)

#simulate first group
g1 <- rmvnorm(n = 40, rep(42,tsf),diag(21, tsf))
#second group
g2 <- rmvnorm(n = 40, rep(46,tsf), diag(21,tsf))
#third group
g3 <- rmvnorm(n = 40, rep(50,tsf), diag(21,tsf))

clear_structure <- rbind(g1, g2, g3) %>% as_tibble()
clear_structure$Group <- rep(c("1","2","3"), times = c(40,40,40)) %>% factor()

## Simulate Multivariate Data
new_dat <- cbind(clear_structure[,-ncol(clear_structure)], no_structure[,1:num_noise_features], clear_structure$Group)
names(new_dat) <- c(paste0("V",1:(tsf + num_noise_features)), "Group")



## Data visualization (parallel coordinate plots)
## some quick visual assessment
ggparcoord(new_dat, 
           columns = 1:(ncol(new_dat)-1),
           groupColumn = ncol(new_dat),
           scale = "globalminmax") +
  labs(x = "Dimension",
       y = "Simulated Value",
       title = "Structure With Noise Features") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_viridis_d(option = "D", end = 0.8)



#### Calculate Clusterings of the Data
final_dat <- new_dat[,-ncol(new_dat)]

# Properly Tune (need to assess Gap at several S, K values). 
gapvals <- list()
wvals <- seq(1.1,5, length = 4)
for(s in 1:length(wvals)){
 gapvals[[s]] <-  clusGap(new_dat[,-ncol(new_dat)], FUN = sparseclust1, K.max = 6, wb = wvals[s], B = 50, spaceH0 = "original", d.power = 2)
}
gapvals[[1]]



df1 <- gaps_mc$Tab  %>% as_tibble()

#df1$W <- rep(wvals, each = 1) %>% factor()
df1$K <- rep(1:6, times = 1)

ggplot(df1, aes(K, gap, color = W)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d(option = "D", end = 0.9) + theme_bw() + ggtitle(paste0("Gap Statistics with # Noise Features: ", sparsity_list[j])) + geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8)


### Generate Clusters


smc1 <- SparseMonoClust(new_dat[,-ncol(new_dat)], wbound = 3)
clusters <- smc1$clustob$membership

smc1$w %>% plot()
smc1$clustob$frame

mc1 <- MonoClust(new_dat[,-ncol(new_dat)], nclusters = 3)
mc1$membership
mc1$frame

## sparse clustering
hc1 <- HierarchicalSparseCluster(as.matrix(new_dat[,-ncol(new_dat)]),  wbound = 3)
hc1$ws %>% plot()
hc1$ws == smc1$w # these are equivalent
cutree(hc1$hc,3)

## Sparse vs. Non-Sparse
par(mfrow= c(1,2))
plot(smc1$clustob); title("Sparse Cluster Splits")
plot(mc1); title("Non-Sparse Splits")
## Classification Error Rate


## Need a way to assess cluster reality (right or wrong) 

# pick a "medoid" 
medoids = new_dat[c(1, 41, 81),]
ggparcoord(medoids, 
           columns = 1:(ncol(new_dat)-1),
           groupColumn = ncol(new_dat),
           scale = "globalminmax") +
  labs(x = "Dimension",
       y = "Simulated Value",
       title = "Medoids") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_viridis_d(option = "D", end = 0.8)

####
df_cluster <- tibble(TrueLabel = rep(1:3,each = 40) %>% factor, SparseMCLab = clusters %>% factor(levels = c(4,5,3), labels = c('1','2','3')), MCLab = mc1$membership %>% factor(levels = c(4,5,3), labels = c('1','2','3')), SparseHCLab = cutree(hc1$hc,3) %>% factor)
df_cluster

## These seem to perfectly cluster - potentially because of the clear structure
ggplot(df_cluster) + geom_point(aes(SparseMCLab, TrueLabel)) + ggtitle("Sparse MC: Labels vs. Truth")
ggplot(df_cluster) + geom_point(aes(MCLab, TrueLabel)) + ggtitle("MC: Labels vs. Truth")
ggplot(df_cluster) + geom_point(aes(SparseHCLab, TrueLabel)) + ggtitle("Sparse HC: Labels vs. Truth")













