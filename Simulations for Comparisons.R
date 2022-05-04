## Paul Harmon
## May 1, 2022
### Simulations: 


library(devtools)
source("SparseClusteringFunctions.R")
library(rCOSA)
library(ggplot2)
library(monoClust)
library(dplyr)
library(tidyr)
library(magrittr)
library(mvtnorm)
library(mnormt)
library(GGally)
library(cluster)
library(sparcl)
library(fossil)

### Simulate Data for Assessment of Sparse Monothetic Clustering 
### Idea is to see how often we recover a better cluster structure in the presence of noise features 

## MonoClust, Sparse Monoclust, and Oracle Method (i.e. using just varying features)
## 

# Goal: To show that sparse monothetic clustering retains 4 clusters more often


checkRand <- function(vec){
  rand.index(as.numeric(as.character(vec)), as.numeric(as.character(df_lab$TrueCluster)))
}

checkCER <- function(vec){
  CER(as.numeric(as.character(vec)), as.numeric(as.character(df_lab$TrueCluster)))
}

#### For Computational Efficiency we reduce # Clusters to 4 ####
## 6 equally sized classes of 20 observations for n = 120
## p = 200
## 20 features differ between classes

### Global simulation parameters
n_sims = 10
wvals <- seq(1.1,8, by = 0.5)

#### INITIAL PARAMATERS ################################################################
total_noise_features <- 190
num_noise_features <- 190
tss <- 10 #stands for total some structure features

data_plots_list <- list()
gap_plot_list <- list()
new_dat_list <- list()

# list of gap values based on different SE rules
gapvals_list <- list()

#cluster labels
mc_out_labels <- mc_oracle_labels <- hc_out_labels <- matrix(0, nrow = 120, ncol = n_sims)
mc_out_list <- mc_oracle_list <- list()
mc_clusgap <- mc_oracle_clusgap <- hc_clusgap <- list()


start = Sys.time()
for(i in 1:n_sims){
  
  
  #### SIMULATE DATA #####################################################################
  #No structure: Simluates 20 rv's with mean 52 and variance 30 (independent)
  #set.seed(71919)
  no_structure <- rmvnorm(n = 120, mean = rep(52,total_noise_features), sigma = diag(30,total_noise_features)) %>% as_tibble()
  
  
  ## Some Structure: Based on the same method as clear structure
  # but this has slightly larger variances and closer means in each dimension
  #set.seed(72019)
  g1_s <- rmvnorm(n = 20, rep(42,tss),diag(20, tss))
  #second group
  g2_s <- rmvnorm(n = 20, rep(46,tss), diag(20,tss))
  #third group
  g3_s <- rmvnorm(n = 20, rep(50,tss), diag(20,tss))
  #4th group
  g4_s <- rmvnorm(n = 20, rep(54,tss),diag(20, tss))
  #5th group
  g5_s <- rmvnorm(n = 20, rep(58,tss),diag(20, tss))
  #6th group
  g6_s <- rmvnorm(n = 20, rep(62,tss),diag(20, tss))
  
  
  some_structure <- rbind(g1_s, g2_s, g3_s,g4_s, g5_s, g6_s) %>% as_tibble()
  some_structure$Group <- rep(c(1:6), each = 20) %>% factor()
  
  
  ## Simulate Multivariate Data
  new_dat <- cbind(some_structure[,c(1:tss)], no_structure[,1:num_noise_features], some_structure$Group)
  names(new_dat) <- c(paste0("V",1:(tss + num_noise_features)), "Group")
  
  #Save data for later
  new_dat_list[[i]] <- new_dat
  
  
  ## Data visualization (parallel coordinate plots)
  ## some quick visual assessment
  data_plots_list[[i]] <- 
    ggparcoord(new_dat, 
               columns = 1:(ncol(new_dat)-1),
               groupColumn = ncol(new_dat),
               scale = "globalminmax") +
    labs(x = "Dimension",
         y = "Simulated Value",
         title = "Structured Data With Additional Noise Features") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_color_viridis_d(option = "D", end = 0.8)
  
  
  #### Generate Monothetic Clusterings on All datasets
  ## Tuning: Set to Optimal Number of Clusters
  mc_clusgap[[i]] <- clusGap(new_dat[,-ncol(new_dat)], FUN = MonoClust1, K.max = 10, B = 50, spaceH0 = "original", d.power = 2)
  mc_out <- MonoClust(as.data.frame(new_dat[,-ncol(new_dat)]), nclusters = 4)
  mc_out_labels[,i] <- mc_out$membership
  mc_out_list[[i]] <- mc_out
  
  #### Generate Oracle Clusterings (only uses the first signal-bearing columns)
  ## Tuning: Set to Optimal Number of Clusters
  mc_oracle_clusgap[[i]] <- clusGap(new_dat[,1:tss], FUN = MonoClust1, K.max = 10, B = 50, spaceH0 = "original", d.power = 2)
  mc_oracle <- MonoClust(as.data.frame(new_dat[,1:tss]), nclusters = 4)
  mc_oracle_labels[,i] <- mc_oracle$membership
  mc_oracle_list[[i]] <- mc_oracle
  
  #### Hclust
  ## Tuning: For S, set to Optimal Number of Clusters
  hc_clusgap[[i]] <-  clusGap(new_dat[,1:tss], FUN = MonoClust1, K.max = 10, B = 50, spaceH0 = "original", d.power = 2)
  hc_out <- hclust(dist(new_dat[,-ncol(new_dat)]), method = "ward.D2")
  hc_out_labels[,i] <- cutree(hc_out, k = 4)
  
  #### Sparse MonoClust
  ## Tuning: Only For S, set to Optimal number of clusters
  

  
  
  #### Tune Simultaneously for S,K ####
#  gapvals_list[[i]] <- list() #initializes list of lists
  # 
  # for(j in 1:length(wvals)){
  #   # Inner For Loop controls the iteration through s, clusGap will iterate through K automatically
  #   gapvals_list[[i]][[j]]<- clusGap(new_dat[,-ncol(new_dat)], FUN = sparseclust1, K.max = 6, wb = wvals[j], B = 50, spaceH0 = "original", d.power = 2)
  #   
  #   
  #   # To Do in additional datasets
  #   #### Generate Sparse Clustering from Optimal K, S For Calculation of CER
  #   #### Save Weights to Dataframe/Tibble 
  #   
  #   
  # } #inner for loop iterating over grid of w-vals for s tuning
  
}#outer for loop iterating over number of simulations

end = Sys.time() #

## Save Outputs For Later Use
write.csv(mc_out_labels, "mc_out_5122.csv")
write.csv(mc_oracle_labels, "mc_oracle_5122.csv")
write.csv(hc_out_labels, "hc_out_5122.csv")
saveRDS()


saveRDS(mc_out_list, paste0(Sys.Date(),"mc_out.RDS"))
saveRDS(mc_oracle_list, paste0(Sys.Date(),"mc_oracle.RDS"))
saveRDS(new_dat_list, paste0(Sys.Date(), "newdata.RDS"))

saveRDS(mc_clusgap, "mc_clusgap.RDS")

#saveRDS(gapvals_list, paste0(Sys.Date(),"gapvals.RDS"))



#### Assessment of Cluster Fit

df_lab <- as_tibble(cbind(mc_oracle_labels, mc_out_labels, hc_out_labels))
names(df_lab) <- rep(c("Oracle", "MC", "HC"), each = n_sims)

df_lab$TrueCluster <- new_dat$Group


## Create Rand Index

#
df_rand <- sapply(df_lab, checkRand) %>% as_tibble()
df_cer <- sapply(df_lab, checkCER) %>% as_tibble()
df_rand$Method <- names(df_lab)
df_cer$Method <- names(df_lab)

A <- df_rand %>% group_by(Method) %>% summarize(Mean = mean(value)) #Rand Index
B <- df_cer %>% group_by(Method) %>% summarize(Mean = mean(value)) #Cluster Error Rate

library(knitr)
dfa <- cbind(A,B)[,c(1,2,4)]
names(dfa) <- c("Method", "Rand", "CER")

kable(dfa)













#Helper Functions
extractWs <- function(smc_list_ob){
  smc_list_ob$w
}

extractClassSparse <- function(smc_list_ob){
  smc_list_ob$clustob$membership
}

extractClass <- function(mc_list_ob){
  mc_list_ob$membership
}


sapply(mc_out_list, extractClass)


















