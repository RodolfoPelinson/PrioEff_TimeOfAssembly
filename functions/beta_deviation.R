#' @title Computing Beta-Deviation
#'
#' @description This function computes observed beta-diversity (community variability), expected beta-diversity and beta-deviation for a given grouping factor. It strongly rely on functions from package `vegan`
#'
#' @param com The community data. A species by site matrix containing abundance or presence/absence data.
#' @param group a grouping factor to compute community variability
#' @param method The chosen null model algorithm. See help file from the `permatswap` function in package `vegan`. Defaut is set to "quasiswap"
#' @param dist The beta-diversity metric to be used. Default is set to bray-curtis beta dissimilarity.
#' @param fixedmar character, stating which of the row/column sums should be preserved ("none", "rows", "columns", "both"). Defaut is to "both".
#' @param shuffle Character, indicating whether individuals ("ind"), samples ("samp") or both ("both") should be shuffled, see details.
#' @param strata Numeric vector or factor with length same as nrow(m) for grouping rows within strata for restricted permutations. Unique values or levels are used. It is usefull to keep gamma-diversity within treatments constant. Defaut is set to "NULL".
#' @param mtype Matrix data type, either "count" for count data, or "prab" for presence-absence type incidence data. Defaut is set to "count".
#' @param times Number of permuted matrices. Defaut is set to 1000.
#' @param burning Number of null communities discarded before proper analysis in sequential ("swap", "tswap") methods. Defaut is set to 0.
#' @param thin Number of discarded permuted matrices between two evaluations in sequential ("swap", "tswap") methods. Defaut is set to 1.
#' @param transform A character. Should the community matrix values be transformed? Default is set to "NULL". Options here are the same as for the `decostand` function from `vegan`.
#' @param seed Do you wish set a seed? Defaut is to "NULL"
#' @export
#'

beta_deviation <- function(com, method = "quasiswap", dist = "bray", fixedmar="both", shuffle = "both",group,
                           strata = NULL,  mtype = "count", times = 1000, burnin = 0, thin = 1, transform = NULL, seed = NULL, type = "centroid", bias.adjust = TRUE, keep_fill = TRUE){
  
  
  require(vegan) #This function strongly rely on package vegan
  
  #It is recommended to set a seed if you want to reproduce exact results multiple times, it can be done without the function using "set.seed()" or inside it using argument "seed"
  if(is.null(seed)==F){
    set.seed(seed)
  }
  
  if(is.null(strata)==FALSE){strata <- as.integer(strata)}
  
  #Generating null communities
  #It uses permatswap function in vegan to generate null communities. See help file from permatswap() for more information on the arguments of this function
    if(isTRUE(keep_fill)){
    nm <- permatswap(com, method = method, fixedmar=fixedmar, shuffle = shuffle,
                   strata = c(strata), mtype = mtype, times = times, burnin = burnin, thin = thin)
  }else{
    nm <- permatfull(com, fixedmar=fixedmar, shuffle = shuffle,
                     strata = c(strata), mtype = mtype, times = times, burnin = burnin, thin = thin)
  }
  
  
  
  
  null_coms<-array(NA, dim = c(dim(as.matrix(vegdist(nm$perm[[1]], method = dist))),times))
  null_distances <- matrix(NA,ncol = nrow(com),nrow = times)
  for(i in 1:length(nm$perm)){
    
    #Here we check if we want to transform abundance data before performing the analysis (for example, to presence/absence or relative abundance).
    if(is.null(transform)==F){
      nm$perm[[i]] <- decostand(nm$perm[[i]], method = transform)
    }
    null_dissim <- vegdist(nm$perm[[i]], method = dist)
    
    #This uses betadisper from vegan to create distances to centroid for each null community generated
    beta_null <- betadisper(null_dissim, group = group, type = type, bias.adjust = bias.adjust)
    null_distances[i,] <- beta_null$distances
    null_coms[,,i] <- as.matrix(null_dissim)
  }
  
  #Computing mean expected distances to centroid and standard deviations
  null_distances_mean <- colMeans(null_distances)
  null_distances_sd <- apply(null_distances, 2, sd)
  
  
  #computing mean dissimilarities of the null communities
  null_exp <- matrix(NA,nrow = nrow(null_coms[,,i]), ncol = ncol(null_coms[,,i]))
  for(i in 1:dim(null_coms)[1]){
    for(j in 1:dim(null_coms)[2]){
      null_exp[i,j]<-mean(null_coms[i,j,])
    }
  }
  
  #computing standard deviation of dissimilarities of the null communities
  null_exp_sd <- matrix(NA,nrow = nrow(null_coms[,,i]), ncol = ncol(null_coms[,,i]))
  for(i in 1:dim(null_coms)[1]){
    for(j in 1:dim(null_coms)[2]){
      null_exp_sd[i,j]<-sd(null_coms[i,j,])
    }
  }
  
  
  #Computing distance to centroid and dissimilarities for the original communities
  if(is.null(transform)==F){
    com_t <- decostand(com, method = transform)
  }else{com_t <- com}
  
  obs_dist <- vegdist(com_t,method = dist)
  obs_beta <-  betadisper(obs_dist, group = group, type = type, bias.adjust = T)
  obs_distances <- obs_beta$distances
  obs_dist <- as.matrix(obs_dist)
  
  #Computing deviations from the null for each distance to centroid observation (one value for each pond)
  dev_distances <- (obs_distances -null_distances_mean)/null_distances_sd
  
  #Computing deviations from the null for each dissimilarity among each pair of communities
  dev_dist<-(obs_dist-null_exp)/null_exp_sd
  dev_dist[is.na(dev_dist)] <- 0
  
  
  
  #Creating and index for what rows and columns from the dissimilarity matrices we want for each factor level
  #We only want dissimilarities within group treatments
  id_dist <- list()
  for (i in 1:length(levels(group))){
    id_dist[[i]] <- which(group == levels(group)[i])
  }
  
  #Creating a vector with the mean expected dissimilarities
  dists_null <- list()
  for(i in 1:length(id_dist)){
    x <-null_exp[id_dist[[i]],id_dist[[i]]]
    dists_null[[i]] <- x[lower.tri(x,diag = F)]
  }
  dists_null_vector <- unlist(dists_null)
  
  #Creating a vector with the observed dissimilarities
  dists_obs <- list()
  for(i in 1:length(id_dist)){
    x <-obs_dist[id_dist[[i]],id_dist[[i]]]
    dists_obs[[i]] <- x[lower.tri(x,diag = F)]
  }
  dists_obs_vector <- unlist(dists_obs)
  
  #Creating a vector with the beta deviations
  deviation <- list()
  for(i in 1:length(id_dist)){
    x <-dev_dist[id_dist[[i]],id_dist[[i]]]
    deviation[[i]] <- x[lower.tri(x,diag = F)]
  }
  deviation_vector <- unlist(deviation)
  
  #Creating a vector with the treatments associated with each dissimilarity value
  dists_treatment <- list()
  for(i in 1:length(id_dist)){
    x <-dev_dist[id_dist[[i]],id_dist[[i]]]
    dists_treatment[[i]] <- rep(levels(group)[i],length(dists_obs[[i]]))
  }
  dists_treatment_vector <- unlist(dists_treatment)
  dists_treatment_vector <- as.factor(dists_treatment_vector)
  
  
  #Converting dissimilarity matrices to distance objects
  dev_dist <- as.dist(dev_dist)
  null_exp <- as.dist(null_exp)
  obs_dist <- as.dist(obs_dist)
  
  return(list(deviation_dissimilarity_matrix = dev_dist,
              expected_dissimilarity_matrix = null_exp,
              observed_dissimilarity_matrix = obs_dist,
              expected_dissimilarity = dists_null_vector,
              observed_dissimilarity = dists_obs_vector,
              deviation_dissimilarity = deviation_vector,
              observed_distances = obs_distances,
              expected_distances = null_distances_mean,
              deviation_distances = dev_distances,
              treatment_factor = dists_treatment_vector,
              null_cons = nm$perm,
              null_distances = null_distances
  ))
  
}
