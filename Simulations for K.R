## Paul Harmon
## April 23, 2022
### Simulations: 

library(devtools)
devtools::install_github('mkampert/rCOSA')

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

### Simulate Data for motivating example
### Idea is to get a single realization with too spares penalty

## 6 equally sized classes of 120 observations
## p = 200
## 20 features differ between classes

#### INITIAL PARAMATERS ################################################################
total_noise_features <- 180
num_noise_features <- 180
tsf <- 20 #stands for total (clear) structure features
tss <- 20 #stands for total some structure features

data_plots_list <- list()
sparsity_list <- c(5,10,50)
gap_plot_list <- list()
gapvals_list <- list()


#for(j in 1:length(sparsity_list)){
  
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
  g2_s <- rmvnorm(n = 30, rep(46,tss), diag(21,tss))
  #third group
  g3_s <- rmvnorm(n = 30, rep(50,tss), diag(21,tss))
  #4th group
  g4_s <- rmvnorm(n = 30, rep(54,tss),diag(21, tss))
  #fifth group
  g5_s <- rmvnorm(n = 30, rep(58,tss), diag(21,tss))
  #sixth group
  g6_s <- rmvnorm(n = 30, rep(62,tss), diag(21,tss))
  
  some_structure <- rbind(g1_s, g2_s, g3_s,g4_s,g5_s,g6_s) %>% as_tibble()
  some_structure$Group <- rep(c(1:6), each = 30) %>% factor()
  
  
  ## Simulate Multivariate Data
  new_dat <- cbind(some_structure[,c(1:tss)], no_structure[,1:num_noise_features], some_structure$Group)
  names(new_dat) <- c(paste0("V",1:(tss + num_noise_features)), "Group")
  
  
  
  ## Data visualization (parallel coordinate plots)
  ## some quick visual assessment
  #data_plots_list[[j]] <- 
  ggparcoord(new_dat, 
             columns = 1:(ncol(new_dat)-1),
             groupColumn = ncol(new_dat),
             scale = "globalminmax") +
    labs(x = "Dimension",
         y = "Simulated Value",
         title = "Structure With Noise Features") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_color_viridis_d(option = "D", end = 0.8)
  

  
#### Simulation Study on Tuning:
  
  #fill up a matrix
  wvals <- seq(1.1,3, length = 5)
  kvals <- 2:4
  gapvals <- list()
  

  
  ### Tune for K First - clusgap on ALL features  
  # Single Run - max K at 6, and iterate through different w values
  start = Sys.time()
    gapvals_1<- clusGap(new_dat[,1:200], FUN = sparseclust1, K.max = 10, wb = wvals[1], B = 50, spaceH0 = "original", d.power = 2)
   end = Sys.time()
  end - start

  
df1 <- gapvals_1$Tab %>% as_tibble()
df1$K <- 1:10

gap_plot <-  ggplot(df1, aes(K, gap)) + geom_line(alpha = .8) + geom_point() + scale_color_viridis_d(option = "D", end = 0.9) + theme_bw() +
              ggtitle("Gap Statistics") + geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, alpha = 0.8) + 
              ylab("Gap Value") + scale_x_continuous(breaks = 1:10, labels = 1:10)

end - start

temp1 <- SparseMonoClust(rawdata = new_dat[,1:200], wbound = wvals[1], nclusters = 4)
new_dat2 <- new_dat

new_dat2$Label <- temp1$clustob$membership
rand.index(as.numeric(as.character(new_dat2$Label)), as.numeric(as.character(new_dat2$Group)))

ggparcoord(new_dat2, 
           columns = 1:20,
           groupColumn = ncol(new_dat),
           scale = "globalminmax") +
  labs(x = "Dimension",
       y = "Simulated Value",
       title = "Structure With Noise Features") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_viridis_d(option = "D", end = 0.8)

df_weight <- tibble(Feature = 1:200, Weight = temp1$w)
weight_plot <- ggplot(df_weight, aes(Feature, Weight)) + theme_bw() + ggtitle("Feature Weights: Too Much Sparsity") + geom_point() + 
  annotate("label", x = 5, y = 1, label = round(df_weight$Weight[5],3), hjust = -.5) + 
  annotate("label", x = 14, y = .02, label = round(df_weight$Weight[14],3), hjust = -.5, vjust = 0.1) + 
  annotate("label", x = 18, y = .09, label = round(df_weight$Weight[18],3), hjust = -1) +
  theme(plot.title = element_text(hjust = 0.5)) 


plot(temp1$w)  
  
### Create Combination Plot with Gap, Weights, and Monothetic Splits

#monothetic Splits
plot(temp1$clustob); title("Monothetic Splits")
frameplot <- recordPlot()

library(ggpubr)
plots2 <- ggarrange(gap_plot, weight_plot, nrow = 1)

library(cowplot)
plot_grid(frameplot, gap_plot, weight_plot , ncol =3, labels = c("","",""))


  
  
### Simultaneous Tuning: Loop assessed at a grid of w-values and iterating over K

  # Single For Loop - max K at 6, and iterate through different w values
  for(i in 1:length(wvals)){
    gapvals[[i]]<- clusGap(new_dat[,1:200], FUN = sparseclust1, K.max = 10, wb = wvals[i], B = 50, spaceH0 = "original", d.power = 2)
  }
  
    
  
  
  
  
  
  
  
  
 #####  Try Clustering Methods ####

  
cosa_out <- cosa2(new_dat[,-ncol(new_dat)])
  
rw_dat <- cosa_out$W * new_dat[,-ncol(new_dat)]
mc_out <- MonoClust(as.data.frame(rw_dat), nclusters = 6)
mc_out %>% plot()

table(new_dat$Group, mc_out$membership)

  
  
  

  
  
  
  
  
  
  
  
  
  
  
