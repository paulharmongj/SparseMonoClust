## Motivating Example plot
set.seed(42222)
x <-  tibble(a = rnorm(n = 500, 0, 1), b = rnorm(n = 500, 0, 1))
y <- tibble(a = rnorm(n = 500, 5, 1), b = rnorm(n = 500, 5, 3))
  
#(5,2),diag(1, 2))

z <- bind_rows(x,y)
names(z) <- c("V1","V2")
z$Group <- rep(c(1,2), each = 500) %>% as.factor()

O <- ggplot(z, aes(V1,V2, shape = Group, color = Group)) + geom_point() + theme_bw() + ggtitle("Actual Clusters") + scale_color_viridis_d(begin = 0.2, end = 0.8) + theme(plot.title = element_text(hjust = 0.5))





mc_out <- MonoClust(z[,1:2], nclusters = 2)
mc_out$frame
z$MC1 <- mc_out$membership %>% as.factor()

A <- ggplot(z, aes(V1,V2, shape = Group, color = MC1)) + geom_point() + theme_bw() + ggtitle("Motivating Example: 2-D Monothetic Clustering") + scale_color_viridis_d(begin = 0.2, end = 0.8) + theme(plot.title = element_text(hjust = 0.5))



source("SparseClusteringFunctions.R")
spmc_out <- SparseMonoClust(z[,1:2], wbound = 1.1, nclusters = 2)
z$SMC1 <- spmc_out$clustob$membership %>% as.factor()
ggplot(z, aes(V1,V2, shape = Group, color = SMC1)) + geom_point() + theme_bw() + ggtitle("Motivating Example: Sparse Monothetiic Clustering") + scale_color_viridis_d(begin = 0.2, end = 0.8) + theme(plot.title = element_text(hjust = 0.5))


oracle_out <- MonoClust(z[,1], nclusters = 2)
z$SMC2 <- oracle_out$membership %>% as.factor()
B <- ggplot(z, aes(V1,V2, shape = Group, color = SMC2)) + geom_point() + theme_bw() + ggtitle("Motivating Example: Sparse Monothetiic Clustering") + scale_color_viridis_d(begin = 0.2, end = 0.8) + theme(plot.title = element_text(hjust = 0.5))



z$KM1 <- kmeans(z[,1:2], centers = 2)$cluster %>% factor()
ggplot(z, aes(V1,V2, shape = Group, color = KM1)) + geom_point() + theme_bw() + ggtitle("Motivating Example")

z$KM1 <- kmeans(z[,1], centers = 2)$cluster %>% factor()
ggplot(z, aes(V1,V2, shape = Group, color = KM1)) + geom_point() + theme_bw() + ggtitle("Motivating Example")

z$GroupLab <- factor(z$Group, levels = c(1,2), labels = c(2,3))

change_num <- function(x){as.character(x) %>% as.numeric()}
results <- z %>%dplyr::select(c(GroupLab, SMC2, MC1))
results_n <- z %>%dplyr::select(c(GroupLab, SMC2, MC1)) %>% sapply(change_num) %>% as_tibble()

library(dplyr)
results %>% mutate(SMCtrue = SMC2 ==GroupLab, MCtrue = MC1 ==GroupLab) %>% dplyr::select(SMCtrue, MCtrue) %>% sapply(sum)
smrand <- rand.index(results_n$GroupLab, results_n$SMC2) %>% round(2)
mcrand <- rand.index(results_n$GroupLab, results_n$MC1)  %>% round(2)

#### Create Plot
library(ggpubr)

A <- ggplot(z, aes(V1,V2, shape = Group, color = MC1)) + geom_point() + theme_bw() + ggtitle("2-D Monothetic Clustering") + 
  scale_color_viridis_d(begin = 0.2, end = 0.8) + theme(plot.title = element_text(hjust = 0.5)) + annotate("label", x = -1, y = 12, fontface = "bold",color = "red", label = paste0("Rand Index: ",mcrand))
B <- ggplot(z, aes(V1,V2, shape = Group, color = SMC2)) + geom_point() + theme_bw() + ggtitle("Sparse Monothetiic Clustering") + 
  scale_color_viridis_d(begin = 0.2, end = 0.8) + theme(plot.title = element_text(hjust = 0.5)) + annotate("label", x = -1, y = 12, fontface = "bold", color = "red",label = paste0("Rand Index: ",smrand))




plot1 <- ggarrange(O,A,B, common.legend = TRUE, nrow = 1)
annotate_figure(plot1, top = text_grob("Motivating Example: Monothetic Clustering", 
                                      color = "black", face = "bold", size = 14))












### Appendix: $$$$$$
# load library MASS
library(MASS)

# set seed and create data vectors
set.seed(98989)
sample_size <- 1000									
sample_meanvector <- c(0, 10)								
sample_covariance_matrix <- matrix(c(1, 0, 2, 1),
                                   ncol = 2)

# create bivariate normal distribution
sample_distribution <- mvrnorm(n = sample_size,
                               mu = sample_meanvector,
                               Sigma = sample_covariance_matrix)

plot(sample_distribution)

