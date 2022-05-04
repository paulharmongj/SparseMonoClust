## Paul Harmon
## April 23, 2022
### Simulations: 

library(devtools)
#devtools::install_github('mkampert/rCOSA')
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

### Simulate Data for Assessment of Gap Statistic
### Idea is to see how often we recover the right number of clusters in the presence of noise features 

## want plots of optimal K's across runs for MonoClust, Sparse Monoclust, and Oracle Method (i.e. using just varying features)
## then show beanplots or boxplots of weights across the runs

# Goal: To show that sparse monothetic clustering retains 4 clusters more often


#### For Computational Efficiency we reduce # Clusters to 4 ####
## 4 equally sized classes of 20 observations for n = 80
## p = 100
## 20 features differ between classes

### Global simulation parameters
n_sims = 10
wvals <- seq(1.1,8, by = 0.5)

#### INITIAL PARAMATERS ################################################################
total_noise_features <- 80
num_noise_features <- 80
tsf <- 20 #stands for total (clear) structure features
tss <- 20 #stands for total some structure features

data_plots_list <- list()
sparsity_list <- c(5,10,50)
gap_plot_list <- list()
new_dat_list <- list()

# list of gap values based on different SE rules
gapvals_list <- list()

start = Sys.time()
for(i in 1:n_sims){
  
  
  #### SIMULATE DATA #####################################################################
  #No structure: Simluates 20 rv's with mean 50 and variance 30 (independent)
  #set.seed(71919)
  no_structure <- rmvnorm(n = 80, mean = rep(48,total_noise_features), sigma = diag(30,total_noise_features)) %>% as_tibble()
  
  
  ## Some Structure: Based on the same method as clear structure
  # but this has slightly larger variances and closer means in each dimension
  #set.seed(72019)
  g1_s <- rmvnorm(n = 20, rep(42,tss),diag(15, tss))
  #second group
  g2_s <- rmvnorm(n = 20, rep(46,tss), diag(15,tss))
  #third group
  g3_s <- rmvnorm(n = 20, rep(50,tss), diag(15,tss))
  #4th group
  g4_s <- rmvnorm(n = 20, rep(54,tss),diag(15, tss))
  
  
  some_structure <- rbind(g1_s, g2_s, g3_s,g4_s) %>% as_tibble()
  some_structure$Group <- rep(c(1:4), each = 20) %>% factor()
  
  
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
  
  
  
  
  
  #### Tune Simultaneously for S,K ####
  gapvals_list[[i]] <- list() #initializes list of lists
  
  for(j in 1:length(wvals)){
    # Inner For Loop controls the iteration through s, clusGap will iterate through K automatically
    gapvals_list[[i]][[j]]<- clusGap(new_dat[,-ncol(new_dat)], FUN = sparseclust1, K.max = 6, wb = wvals[j], B = 50, spaceH0 = "original", d.power = 2)
    
   
    # To Do in additional datasets
    #### Generate Sparse Clustering from Optimal K, S For Calculation of CER
    #### Save Weights to Dataframe/Tibble 
    
    
  } #inner for loop iterating over grid of w-vals for s tuning
  
}#outer for loop iterating over number of simulations

end = Sys.time() #takes about 3 hours to run in this iteration

## Save Outputs For Later Use
saveRDS(gapvals_list, paste0(Sys.Date(),"gapvals.RDS"))
saveRDS(newdat_list, paste0(Sys.Date(), "newdata.RDS"))



# Computation time
end - start


#### Read in Results from Simulations
new_dat_list <- readRDS("simulation425_data.rds")
#gapvals_list <- 

#CSV of gap vals
df_total <- read.csv("GapStatisticsData424.csv")



# Assess best GAP across each of the different methods across all the simulations
kmat_global <- matrix(0, ncol = length(gapvals_list[[1]]), nrow = length(gapvals_list))
kmat_first <- matrix(0, ncol = length(gapvals_list[[1]]), nrow = length(gapvals_list))
kmat_tibs <- matrix(0, ncol = length(gapvals_list[[1]]), nrow = length(gapvals_list))
kmat_firstSE <- matrix(0, ncol = length(gapvals_list[[1]]), nrow = length(gapvals_list))
kmat_globalSE <- matrix(0, ncol = length(gapvals_list[[1]]), nrow = length(gapvals_list))
colnames(kmat_global) <- colnames(kmat_first) <- colnames(kmat_tibs) <- colnames(kmat_firstSE) <- colnames(kmat_globalSE) <-wvals
rownames(kmat_global) <- rownames(kmat_first) <- rownames(kmat_tibs) <- rownames(kmat_firstSE) <- rownames(kmat_globalSE) <- 1:10

## double nested for loop to save a table with the info
for(i in 1:length(gapvals_list)){
  for(j in 1:length(gapvals_list[[i]])){
    df_gap <- gapvals_list[[i]][[j]]$Tab %>% as_tibble()
    kmat_global[i,j] <- maxSE(df_gap$gap, df_gap$SE.sim,  method = "globalmax")
    kmat_first[i,j] <- maxSE(df_gap$gap, df_gap$SE.sim,  method = "firstmax")
    kmat_tibs[i,j] <- maxSE(df_gap$gap, df_gap$SE.sim,  method = "Tibs2001SEmax")
    kmat_firstSE[i,j] <- maxSE(df_gap$gap, df_gap$SE.sim,  method = "firstSEmax")
    kmat_globalSE[i,j] <- maxSE(df_gap$gap, df_gap$SE.sim,  method = "globalSEmax")
  }
}



###Create a Faceted View
df_total <- rbind(kmat_global, kmat_first, kmat_tibs, kmat_firstSE, kmat_globalSE) %>% as_tibble()
df_total$Method <- rep(c("Global","First", "Tibshirani", "FirstSE", "GlobalSE"), each = 10)
#names(df_total) = c("S1","S1.575", "S2.05", "S2.525", "S3", "Method")

df_total %>% pivot_longer(1:5, names_to = "S_Value") %>% 
ggplot() + geom_bar(aes(factor(value), group = Method), stat = 'count') + facet_wrap(~Method) + 
  ggtitle("Optimal K for All Sparsity Values")  + theme_bw() + theme(plot.title = element_text(hjust = 0.5))


### Create Plots for Simultaenous Tuning

plot_ClusGap <- function(gapvals_list){
  #input takes a list of gapvals (inner list)
  extract_tab <- function(clusgap) {as_tibble(clusgap$Tab, .name_repair = 'unique')}
  df_temp <- lapply(gapvals_list, extract_tab) %>% bind_rows()
  # df <- rbind(gapvals_temp[[1]]$Tab,
  #              gapvals_temp[[2]]$Tab,
  #              gapvals_temp[[3]]$Tab,
  #              gapvals_temp[[4]]$Tab,
  #              gapvals_temp[[5]]$Tab) %>% as_tibble()
  
  df_temp$W <- rep(wvals, each = 6) %>% factor()
  df_temp$K <- rep(1:6, times = length(wvals))
  
  return(df_temp)
}

df_total <- lapply(gapvals_list,plot_ClusGap) %>% bind_rows()
df_total$Simulation <- rep(1:10, each = 84)
  
#Plot of K with lines colored by W
  ggplot(df_total, aes(K, gap, color = W)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d(option = "D", end = 0.9) + 
    theme_bw() + ggtitle("Gap Statistics Across 10 Runs") + ylab("Gap Statistic") + xlab("K: Number of Clusters") + 
    geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8) + 
    facet_wrap(~Simulation)
  
#Plot of W with lines colored by K (like Brodinova et al, 2019)
  ggplot(df_total, aes(W, gap, color = factor(K), group = K)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d("K", option = "D", end = 0.9) + 
    theme_bw() + ggtitle("Gap Statistics Across 10 Runs") + ylab("Gap Statistic") + xlab("S: Sparsity Parameter") + 
    geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8) + 
    facet_wrap(~Simulation)
  

#print only the first simulation
a<-   df_total %>% filter(Simulation == 1) %>% 
    ggplot(aes(K, gap, color = W, group = W)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d(option = "D", end = 0.9) + 
    theme_bw() + ggtitle("Gap Statistics Across K Values") + ylab("Gap Statistic") + xlab("K: Number of Clusters") + 
    geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8) + theme(plot.title = element_text(hjust = 0.5))
  
b <-   df_total %>% filter(Simulation == 1) %>%
    ggplot(aes(W, gap, color = factor(K), group = K)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d("K", option = "D", end = 0.9) + 
    theme_bw() + ggtitle("Gap Statistics across S Values") + ylab("Gap Statistic") + xlab("S: Sparsity Parameter") + 
    geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8) + theme(plot.title = element_text(hjust = 0.5))
  
  df_total %>% filter(Simulation == 1) %>%
    ggplot(aes(W, gap, color = factor(K), group = K)) + geom_point() + scale_color_viridis_d("K", option = "D", end = 0.9) + 
    theme_bw() + ggtitle("Gap Statistics Across 1 Run") + ylab("Gap Statistic") + xlab("S: Sparsity Parameter") #
  #+  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8) 
  
library(ggpubr)  
plot <- ggarrange(a,b)  
annotate_figure(plot, top = text_grob("Gap(s,k) Statistics Presented In Two Ways ", 
                                      color = "black", face = "bold", size = 14))
  
  # Ideally, we could see the optimal choice from each simulation
df_total %>% group_by(Simulation) %>% dplyr::max()
  

#save df_total as data frame
write.csv(df_total, "GapStatisticsData424.csv")


#function to determine optimal 
#df_gap needs to have gap, gapWE, W, and K defined
df_test <- filter(df_total, Simulation == 1)
multiGap <- function(df_gap){
  # order by W and K
  df_temp <- df_gap %>% arrange(W,K) %>% mutate(Lower = gap - SE.sim, Upper = gap + SE.sim)
  
  #Global Max:
  df_temp[which.max(df_temp$gap),]
  
  #First Max:
  df_temp %>% group_by(W) %>%
    mutate(maxGap = max(gap), maxGapName = paste0("NumClusters: ", which.max(gap))) %>% View()
  
  #Global SE Max (Dudoit)
  
  
  #Tibshirani 2001: Smallest K such that f(k) is >= f(k+1) - s(K+1)
  
  
}




### Plots the Features By level of Sparsity (and potentially based on K)
#plot the features weights For each of the optimal methods by simulation
#new_dat_list[[1]]
#consider K fixed at 4 for now
smc_out_list <- list()
wdf_list <- list()

weight_plot_list <- list()
weight_mat_list <- list() #list of plot matrices
rand_mat <- matrix(0, nrow = length(wvals), ncol = 10) #instantiate a matrix to store rand indexes by s val

#Helper Functions
extractWs <- function(smc_list_ob){
  smc_list_ob$w
}

extractClass <- function(smc_list_ob){
  smc_list_ob$clustob$membership
}

library(fossil)

for(i in 1:length(wvals)){
  
  #inner for-loop generates a list of sparse monoclust objects
  for(j in 1:length(new_dat_list)){
    smc_out_list[[j]] <- SparseMonoClust(new_dat_list[[j]][,-ncol(new_dat)], nclusters = 4, wbound = wvals[i])
  }
  

  #Extracts Wvalues and generates a plot that runs over the 10 simulated datasets
  Wvals <- sapply(smc_out_list, extractWs)
  
  wdf <- Wvals %>% as_tibble()
  wdf_list[[i]] <- wdf 
  #Plot
  weight_plot_list[[i]] <- wdf %>% pivot_longer(1:10, names_to = "Simulation") %>% mutate(Feature = rep(1:100, each = 10)) %>%
    ggplot(aes(Feature, value, group = Feature, fill = factor(Feature))) + geom_boxplot() + ggtitle("Feature Weights") + 
    geom_vline(xintercept = 20.5, linetype = 2, color = "tomato2") + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + guides(fill = "none")
  

  #Generate classes from each of the 10 sparse monoclusts
  df_class <- sapply(smc_out_list, extractClass) %>% as_tibble()
  df_class$TrueGroup <- new_dat$Group %>% as.character()
  
  
  rand.vec <- rep(0, 10)
  for(j in 1:10){
    col_name <- colnames(df_class)
    rand.vec[j] <- rand.index(as.numeric(unlist(df_class[,j])), as.numeric(df_class$TrueGroup))
  }
  
  rand_mat[i,] <- rand.vec 
  
  # Let's Compare
  #rand.vec1 <- rand.vec #s = 2.1
  #rand.vec2 <- rand.vec #s = 5
  #rand.vec3 <- rand.vec #s = 1.1
  #rand.vec4 <- rand.vec # s = 10
}  


### Create Plot of Output/Rand Values
  df_rv <- as_tibble(rand_mat)
  names(df_rv) <- paste0("Sim",1:10)
  df_rv$S = rep(wvals)
  
  df_rv %>% pivot_longer(1:10, names_to = "Sim", values_to = "Rand") %>% 
  ggplot() + geom_density(aes(Rand, color = factor(S), group = S), alpha = 0.6) + 
    ggtitle("Density of Rand Index Values") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) 
  





#geom_histogram(aes(Rand, y=..density.., fill = factor(S), group = S), alpha = 0.6, bins = 15, position = "identity") + 


