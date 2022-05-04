### Sparse Clustering Functions

SparseMonoClust <- function(rawdata, wbound = 1.1, dissimilarity = "squared.distance", nclusters = 3, scale1 = FALSE, warnings = FALSE, divclust = FALSE, PMD = TRUE, COSA = FALSE){
  ## Define a l2norm function
  l2n <- function(vec){
    return(sqrt(sum(vec^2)))
  }
  
  ## Check Packages
  if(!require(sparcl)){install.packages('sparcl');library(sparcl)}
  if(!require(monoClust)){install.packages('monoClust');library(monoClust)}
  
  ## Optional Scaling of the original data
  if(scale1 == TRUE){data <- apply(rawdata, 2,scale)} else{data <- rawdata}
  ## Generates a warning that the units are standardized so no longer on the same units
  if(warnings == TRUE){print("Hey! You're not going to be on the original units!")}
  
  ## First Step: Perform Sparse Clustering to Get W Vector with Sparcl
  
  if(COSA == TRUE){
    data_weight <- cosa2(as.matrix(data))$W * data  #obtains the sparse dissimilarity weights from COSA and multiplies weights times original dataset
    
  }
  else{
  hc1 <- HierarchicalSparseCluster(as.matrix(data), wbound = wbound, dissimilarity = dissimilarity, silent = TRUE)
  
  
  ## Second Step: Calculate the Distance Matrix and L2Norm (to match scaling of the 'sparcl' u matrix)
  
  dist2 <- function(vec){(dist(vec)^2)}
  ds <- apply(data, 2, dist2)
  weighted_feature_dists <- sweep(ds, 2, hc1$ws, '*')
  single_weighted_dist <- apply(weighted_feature_dists, 1, sum)
  
  u2 <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  u2[lower.tri(u2)] <- single_weighted_dist
  u2[upper.tri(u2)] <- single_weighted_dist
  normed_u2 <- u2/l2n(c(u2)) #generates the normed u vector (identical to the hcl$u)
  
  ## Now - re-weight the original data with both the weights and the l2norm
  ws_final <- sqrt(hc1$ws/l2n(c(u2)))
  data_weight <- sweep(data, 2, ws_final, '*')}
  
  
  
  
  ## Passes the weighted data into MonoClust or DivClust depending on divclust argument
if(divclust == TRUE){
  #sparsedc1$weighted_data[,c(which(sapply(sparsedc1$weighted_data, sum) >0))
  mc <- divclust(as.data.frame(data_weight[,which(sapply(data_weight,sum)>0)]), K = nclusters)} #uses divclust
else{
  mc <- MonoClust(as.data.frame(data_weight), nclusters = nclusters)
}
  ## Post-Processing: Either another step in this function, or a method applied afterward
  
return(list(clustob = mc, u = normed_u2, weighted_data = data_weight, w = hc1$ws, sparclob = hc1))}

## For input into ClusGap
####
sparseclust1 <- function(x, k, wb){
  # Calcluates the sparse clustering for use in clusgap
  temp <- SparseMonoClust(rawdata = x, wbound = wb, nclusters = k) 
  Membership <- temp$clustob$membership
  return(list(cluster = Membership))}

# sparsedivclust1 <- function(x, k, wb){
#   # Calcluates the sparse clustering for use in clusgap
#   temp <- SparseMonoClust(rawdata = x, wbound = wb, nclusters = k, divclust = TRUE) 
#   Membership <- temp$clustob$membership
#   return(list(cluster = Membership))}

sparseclust2 <- function(x, k, wb){
  temp <- HierarchicalSparseCluster(x, wbound = wb, silent = TRUE, method = "complete")
  Membership <- cutree(temp$hc, k = k)
  return(list(cluster = Membership))
}

MonoClust1 <- function(x, k){
  temp <- MonoClust(as.data.frame(x), nclusters = k)
  Membership <- temp$membership
  return(list(cluster = Membership))}
