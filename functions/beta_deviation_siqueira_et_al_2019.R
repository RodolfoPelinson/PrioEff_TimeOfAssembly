#' @title Computing Beta-Deviation
#'
#' @description This function computes observed beta-diversity (community variability), expected beta-diversity and beta-deviation for a given grouping factor. It is a changed version of the "beta_deviation" function to use the same null model as used in "Siqueira et al. 2020"
#'
#' @param com The community data. A species by site matrix containing abundance or presence/absence data.
#' @param group a grouping factor to compute community variability
#' @param keep.gamma Should permutations be restrict to each treatment?
#' @param dist The beta-diversity metric to be used. Default is set to bray-curtis beta dissimilarity.
#' @param times Number of permuted matrices. Default is set to 1000.
#' @param transform A character. Should the community matrix values be transformed? Default is set to "NULL". Options here are the same as for the `decostand` function from `vegan`.
#' @param seed Do you wish set a seed? Default is to "NULL"
#' @export
#' @references #'Siqueira, Tadeu, Victor S. Saito, Luis M. Bini, Adriano S. Melo, Danielle K. Petsch, Victor L. Landeiro, Kimmo T. Tolonen, Jenny Jyrkänkallio‐Mikkola, Janne Soininen, and Jani Heino. “Community Size Can Affect the Signals of Ecological Drift and Niche Selection on Biodiversity.” Ecology 101, no. 6 (2020): e03014. \url{https://doi.org/10.1002/ecy.3014}.
#'

beta_deviation_siqueira_et_al_2019 <- function(com, dist = "bray", keep.gamma = NULL, group,
                                               times = 1000, seed = NULL, transform = NULL, thin = NULL){
  
  require(vegan) #This function strongly rely on package vegan
  
  #It is recommended to set a seed if you want to reproduce exact results multiple times, it can be done without the function using "set.seed()" or inside it using argument "seed"
  if(is.null(seed)==F){
    set.seed(seed)
  }
  
  #Generating null communities
  #It uses permatswap function in vegan to generate null communities. See help file from permatswap() for more information on the arguments of this function
  ###############################################################################
  nulls_frac <- list()
  if(is.null(keep.gamma) == FALSE){
    levels <- levels(keep.gamma)
    nlevels <- length(levels(keep.gamma))
    for(i in 1:nlevels){
      frac_com <- com[which(keep.gamma == levels[i]),]
      new_com <- as.data.frame(matrix(nrow = nrow(frac_com) * ncol(frac_com), ncol = 3))
      
      colnames(new_com) <- c("species", "sites", "abundance")
      for(j in 1:ncol(frac_com)){
        new_com$species[(1+(nrow(frac_com))*(j-1)):(nrow(frac_com)*j)] <- rep(colnames(frac_com)[j],nrow(frac_com))
        new_com$sites[(1+(nrow(frac_com))*(j-1)):(nrow(frac_com)*j)] <- rownames(frac_com)
        new_com$abundance[(1+(nrow(frac_com))*(j-1)):(nrow(frac_com)*j)] <- frac_com[,j]
      }
      new_com$species <- as.numeric(as.factor(new_com$species))
      new_com$sites <- as.numeric(new_com$sites)
      new_com_expanded <- new_com[rep(seq.int(1, nrow(new_com)), new_com$abundance), ]
      
      nulls_frac[[i]] <- replicate(times, null_assembly(new_com_expanded, species.col = "species", site.col = "sites"), simplify = F)
    }
    
    
    cols <- as.numeric(as.factor(colnames(com)))
    for(i in 1:nlevels){
      for (j in 1:times){
        null_frac_cols <- as.numeric(colnames(nulls_frac[[i]][[j]]))
        missing <- setdiff(cols,null_frac_cols)
        add <- matrix(0,nrow = dim(nulls_frac[[i]][[j]])[1], ncol = length(missing))
        colnames(add) <- missing
        nulls_frac[[i]][[j]]<-cbind(nulls_frac[[i]][[j]], add)
        nulls_frac[[i]][[j]]<-nulls_frac[[i]][[j]][,order(as.numeric(colnames(nulls_frac[[i]][[j]])))]
      }
    }
    
    nulls <- list()
    for(i in 1:times){
      null_temporary<- rbind(nulls_frac[[1]][[i]], nulls_frac[[2]][[i]])
      for(j in 3:nlevels){
        null_temporary<- rbind(null_temporary, nulls_frac[[j]][[i]])
      }
      null_temporary<-null_temporary[,cols]
      null_temporary<-null_temporary[order(as.numeric(rownames(null_temporary))),]
      nulls[[i]] <- null_temporary
    }
    
    
  }else{
    new_com <- as.data.frame(matrix(nrow = nrow(com) * ncol(com), ncol = 3))
    colnames(new_com) <- c("species", "sites", "abundance")
    for(i in 1:ncol(com)){
      new_com$species[(1+(nrow(com))*(i-1)):(nrow(com)*i)] <- rep(colnames(com)[i],nrow(com))
      new_com$sites[(1+(nrow(com))*(i-1)):(nrow(com)*i)] <- rownames(com)
      new_com$abundance[(1+(nrow(com))*(i-1)):(nrow(com)*i)] <- com[,i]
    }
    new_com$species <- as.numeric(as.factor(new_com$species))
    new_com$sites <- as.numeric(new_com$sites)
    new_com_expanded <- new_com[rep(seq.int(1, nrow(new_com)), new_com$abundance), ]
    
    nulls <- replicate(times, null_assembly(new_com_expanded, species.col = "species", site.col = "sites"), simplify = F)
  }
  
  
  
  null_coms<-array(NA, dim = c(dim(as.matrix(vegdist(nulls[[1]], method = dist))),times))
  null_distances <- matrix(NA,ncol = nrow(com),nrow = times)
  for(i in 1:length(nulls)){
    
    #Here we check if we want to transform abundance data before performing the analysis (for example, to presence/absence or relative abundance).
    if(is.null(transform)==F){
      nulls[[i]] <- decostand(nulls[[i]], method = transform)
    }
    null_dissim <- vegdist(nulls[[i]], method = dist)
    
    #This uses betadisper from vegan to create distances to centroid for each null community generated
    beta_null <- betadisper(null_dissim, group = group, type = "centroid", bias.adjust = T)
    null_distances[i,] <- beta_null$distances
    null_coms[,,i] <- as.matrix(null_dissim)
  }
  
  #Computing mean expected distances to centroid and standard deviations
  null_distances_mean <- colMeans(null_distances)
  null_distances_sd <- apply(null_distances, 2, sd)
  
  
  #computing mean dissimilarities of the null communities
  null_exp <- matrix(NA,nrow = dim(null_coms)[1], ncol = dim(null_coms)[2])
  for(i in 1:dim(null_coms)[1]){
    for(j in 1:dim(null_coms)[2]){
      null_exp[i,j]<-mean(null_coms[i,j,])
    }
  }
  
  #computing standard deviation of dissimilarities of the null communities
  null_exp_sd <- matrix(NA,nrow = dim(null_coms)[1], ncol = dim(null_coms)[2])
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
  obs_beta <-  betadisper(obs_dist, group = group, type = "centroid", bias.adjust = T)
  obs_distances <- obs_beta$distances
  obs_dist <- as.matrix(obs_dist)
  
  #Computing deviations from the null for each distance to centroid observation (one value for each pond)
  dev_distances <- (obs_distances - null_distances_mean)/null_distances_sd
  
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
              null_cons = nulls,
              null_distances = null_distances
  ))
  
}
