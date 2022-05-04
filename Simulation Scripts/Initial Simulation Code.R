## Paul Harmon

## Simulation Studies for Sparse Monothetic Clustering
## Scenarios to Consider
# + Independent 3 multivariate normal groups with varying signal to Noise Ratio
#     + 100% Signal
#     + 50% Signal
#     + 10% Signal


# Gap Statistic - Does it Recover the same number of groups?
# Does it filter out the noise features? 
# Cluster Comparison? 
library(ggplot2)
library(monoClust)
library(dplyr)
library(tidyr)
library(magrittr)
library(mvtnorm)
library(mnormt)
library(GGally)
library(cluster)
source("C:/Users/paulh/Documents/Doctoral Work/Sparcl Witten/sparcl/sparclproject/SparseClusteringFunctions.R")


#### INITIAL PARAMATERS ################################################################
total_noise_features <- 100
num_noise_features <- 10
tsf <- 10 #stands for total (clear) structure features
tss <- 10 #stands for total some structure features

data_plots_list <- list()
sparsity_list <- c(5,10,50)
gap_plot_list <- list()
gapvals_list <- list()


for(j in 1:length(sparsity_list)){
  
  num_noise_features <- sparsity_list[j]
  
  
  #### SIMULATE DATA #####################################################################
  #No structure: Simluates 20 rv's with mean 50 and variance 25 (independent)
  #set.seed(71919)
  no_structure <- rmvnorm(n = 90, mean = rep(50,total_noise_features), sigma = diag(25,total_noise_features)) %>% as_tibble()
  
  ###Clear Structure
  #set.seed(71919)
  #simulate first group
  g1 <- rmvnorm(n = 30, rep(40,tsf),diag(8, tsf))
  #second group
  g2 <- rmvnorm(n = 30, rep(50,tsf), diag(8,tsf))
  #third group
  g3 <- rmvnorm(n = 30, rep(60,tsf), diag(8,tsf))
  
  clear_structure <- rbind(g1, g2, g3) %>% as_tibble()
  clear_structure$Group <- rep(c("1","2","3"), times = c(30,30,30)) %>% factor()
  
  
  ## Some Structure: Based on the same method as clear structure
  # but this has slightly larger variances and closer means in each dimension
  #set.seed(72019)
  g1_s <- rmvnorm(n = 30, rep(42,tss),diag(21, tss))
  #second group
  g2_s <- rmvnorm(n = 30, rep(46,tss), diag(25,tss))
  #third group
  g3_s <- rmvnorm(n = 30, rep(50,tss), diag(21,tss))
  
  some_structure <- rbind(g1_s, g2_s, g3_s) %>% as_tibble()
  some_structure$Group <- rep(c("1","2","3"), times = c(30,30,30)) %>% factor()
  
  
  ## Simulate Multivariate Data
  new_dat <- cbind(clear_structure[,c(1:10)], no_structure[,1:num_noise_features], clear_structure$Group)
  names(new_dat) <- c(paste0("V",1:(10 + num_noise_features)), "Group")
  
  
  
  ## Data visualization (parallel coordinate plots)
  ## some quick visual assessment
  data_plots_list[[j]] <- ggparcoord(new_dat, 
                                     columns = 1:(ncol(new_dat)-1),
                                     groupColumn = ncol(new_dat),
                                     scale = "globalminmax") +
    labs(x = "Dimension",
         y = "Simulated Value",
         title = "Structure With Noise Features") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_color_viridis_d(option = "D", end = 0.8)
  
  
  #### Assess Optimal Gap and Sparsity Tuning
  
  #fill up a matrix
  wvals <- seq(1.1,3, length = 5)
  kvals <- 2:5
  gapvals <- list()
  
  
  # Single For Loop - max K at 6, and iterate through different w values
  for(i in 1:length(wvals)){
    gapvals[[i]]<- clusGap(new_dat[,-ncol(new_dat)], FUN = sparseclust1, K.max = 6, wb = wvals[i], B = 50, spaceH0 = "original", d.power = 2)
  }
  gapvals_list[[j]] <- gapvals #tracks across different runs
  
  ### Create Plot
  
  df1 <- rbind(gapvals[[1]]$Tab,
               gapvals[[2]]$Tab,
               gapvals[[3]]$Tab,
               gapvals[[4]]$Tab,
               gapvals[[5]]$Tab) %>% as_tibble()
  
  df1$W <- rep(wvals, each = 6) %>% factor()
  df1$K <- rep(1:6, times = 5)
  
  gap_plot_list[[j]] <- ggplot(df1, aes(K, gap, color = W)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d(option = "D", end = 0.9) + theme_bw() + ggtitle(paste0("Gap Statistics with # Noise Features: ", sparsity_list[j])) + geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8)
  
}


### Assess Gap Statistics and Other Outputs
library(ggpubr)
ggarrange(plotlist = gap_plot_list, ncol = 3, nrow = 1) 









