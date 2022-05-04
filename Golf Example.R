# Read in data

library(readr)
pgolf <- read_csv("C:/Users/paulh/Downloads/Paul Riverside Tracked Rounds - Riverside.csv")
View(pgolf)



library(monoClust)





topar <- as.matrix(pgolf[-1,2:19]) - matrix(rep(unlist(pgolf[1,2:19]), nrow(pgolf)-1), ncol = 18) %>% as_tibble()
apply(pgolf[,2:19], 1, sum)
apply(topar, 1, sum)

topar$Year <- pgolf$Year[-1]


# Temporal Changes

tp20 <- dplyr::filter(topar, Year == 2020)
tp21 <- dplyr::filter(topar, Year == 2021)
mc_all <- MonoClust(topar[,1:18], nclusters = 3)
mc_20 <- MonoClust(tp20[,1:18], nclusters = 3)
mc_21 <- MonoClust(tp21[,1:18], nclusters = 3)

smc_all <- SparseMonoClust(topar[,1:18], nclusters = 3, wbound = 1.1)
plot(smc_all$clustob)
plot(smc_all$w)

par(mfrow = c(1,3))
plot(mc_all)
plot(mc_20)
plot(mc_21)

#create plots
topar$MC1 <- mc_all$membership %>% factor()
topar$ID <- 1:nrow(topar)
topar %>% pivot_longer(1:18, names_to = "Hole", values_to = "Score") %>% 
  ggplot(aes(Hole, Score, group = ID, color = MC1)) + geom_point() + geom_smooth(fill = NA) + ggtitle("All Scores")


# 2020
#create plots
tp20$MC1 <- mc_20$membership %>% factor()
tp20$ID <- 1:nrow(tp20)
tp20 %>% pivot_longer(1:18, names_to = "Hole", values_to = "Score") %>% 
  ggplot(aes(Hole, Score, group = ID, color = MC1)) + geom_point() + geom_smooth(fill = NA) + ggtitle("2020 Scores")




# 2021
#create plots
tp21$MC1 <- mc_21$membership %>% factor()
tp21$ID <- 1:nrow(tp21)
tp21 %>% pivot_longer(1:18, names_to = "Hole", values_to = "Score") %>% 
  ggplot(aes(Hole, Score, group = ID, color = MC1)) + geom_point() + geom_smooth(fill = NA) + ggtitle("2021 Scores")



















