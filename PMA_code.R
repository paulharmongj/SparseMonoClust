# ------------------------------- #
#### R script to play with PMA ####
# ------------------------------- #
install.packages('PMA')
library(PMA)

#bring in a dataset
library(readr)
library(tibble)
library(magrittr)
library(dplyr)
predicted_house <- read_csv("~/predicted_house.csv")

predicted_house_numeric <- predicted_house %>% select(YEAR_BUILT,BEDROOMS, TOT_BATH,ACRES)
predicted_house_numeric_scale <- predicted_house_numeric %>% scale(center = TRUE, scale = TRUE)
# --------------------- # 
#### PCA Comparisons ####
# --------------------- #

#standard pca - already standardized
pc1 <- prcomp(predicted_house_numeric_scale, center = FALSE, scale = FALSE)
summary(pc1)


#sparse pca using Witten's PMA package - no optimization/tuning
pc_sparse1 <- SPC(predicted_house_numeric_scale, K = 3, sumabsv = 1.5)
pc_sparse1$prop.var.explained


#play with the tuning a bit
pcv <- SPC.cv(predicted_house_numeric_scale, nfolds = 5, niter = 10, sumabsv = seq(1,1.9, by = 0.1))
plot(pcv)
pcv$bestsumabsv
pcv$bestsumabsv1se



