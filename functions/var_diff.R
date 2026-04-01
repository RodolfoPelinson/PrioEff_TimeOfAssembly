com <- comm_AM4_drop_atrasado
groups <- Exp_design_drop_atrasado$treatments
type <- "median"
iter <- 100
method = "bray"


var_diff <- function(com = NULL, lvs = NULL, groups, method, type, iter, ...){
  
  #Convert groups to factor
  groups <- as.factor(groups)
  
  #Find group levels
  groups_levels <- levels(groups)
  
  
  
  if(is.null(com) == FALSE){
    #Compute distance matrix for all communities
    dist_com <- vegdist(com, method = method)
    
    #Compute centroid distances for all groups
    cent_dists <- betadisper(dist_com, group = groups, type = type, ...)
    
    #PCoAs
    vectors <- cent_dists$vectors
  }
  
  else{
    vectors <- lvs
  }

  
  
  ####### Em cada tratamento, tenho que pegar 4, fazer a média, guardar e fazer a média das médias.
  
  
  diff_cent <- rep(NA, iter)
  for(p in 1:iter){
    
    #### Compute the WITHIN groups variability ###########################################

    mean_within_cent_groups <- rep(NA, length(groups_levels))
    for(i in 1:length(groups_levels)){
      
      within_vectors <- matrix(NA, nrow = length(groups_levels), ncol = ncol(vectors))
      colnames(within_vectors) <- colnames(vectors)
      within_vectors <- data.frame(within_vectors)
      
      #take communities from specific group levels
      subsample_vectors <- vectors[groups == groups_levels[i],]
      
      within_vectors <- subsample_vectors[sample(1:nrow(subsample_vectors), length(groups_levels)),]
      
      if(type == "centroid"){
        centroid <- colMeans(within_vectors)
      }
      
      if(type == "median"){
        centroid <- apply(within_vectors, MARGIN = 2, median)
      }
      
      within_vectors_cent <-  rbind(within_vectors, centroid)
      
      dist_within_vectors <- as.matrix(dist(within_vectors_cent, method = "euclidean"))
      
      dist_cent <- dist_within_vectors[nrow(dist_within_vectors),1:(ncol(dist_within_vectors)-1)]
      
      if(type == "centroid"){
        mean_within_cent_groups[i] <- mean(dist_cent)
      }
      
      if(type == "median"){
        mean_within_cent_groups[i] <- median(dist_cent)
      }
      
    }
    
    
    if(type == "centroid"){
      mean_within_cent <- mean(mean_within_cent_groups)
    }
    
    if(type == "median"){
      mean_within_cent <- median(mean_within_cent_groups)
    }    
    
    
    
    
    #### Compute the AMONG groups variability ############################################
    among_vectors <- matrix(NA, nrow = length(groups_levels), ncol = ncol(vectors))
    colnames(among_vectors) <- colnames(vectors)
    among_vectors <- data.frame(among_vectors)
    
    for(i in 1:length(groups_levels)){
      
      #take communities from specific group levels
      subsample_vectors <- vectors[groups == groups_levels[i],]
      
      among_vectors[i,] <- subsample_vectors[sample(1:nrow(subsample_vectors), 1),]
      
    }
    if(type == "centroid"){
    centroid <- colMeans(among_vectors)
    }
    
    if(type == "median"){
      centroid <- apply(among_vectors, MARGIN = 2, median)
    }
    
    among_vectors_cent <-  rbind(among_vectors, centroid)
    
    dist_among_vectors <- as.matrix(dist(among_vectors_cent, method = "euclidean"))
    
    dist_cent <- dist_among_vectors[nrow(dist_among_vectors),1:(ncol(dist_among_vectors)-1)]
    
    if(type == "centroid"){
      mean_among_cent <- mean(dist_cent)
    }
    
    if(type == "median"){
      mean_among_cent <- median(dist_cent)
    }
    
    
    
    ##### Compute difference ######################
    
    diff_cent[p] <- mean_among_cent - mean_within_cent
    
  }
  
  
  mean_diff_cent <- mean(diff_cent)
  p_greater_cent <-  length(diff_cent[diff_cent<=0])/iter
  #quantile(diff_cent, c(0.025, 0.975))
  CI <- quantile(diff_cent, c(0.05, 1))
  results_cent <- c(mean = mean_diff_cent, LCL = CI[1], UCL = CI[2], p_value = p_greater_cent)
  
  return(results_cent)
  
  
  
  
  
  #Distance rout
  
  #mean_treat_dist <- rep(NA, length(groups_levels))
  #compute mean dist for each treatment
  #for(i in 1:length(mean_treat_dist)){
    
    #take communities from specific group levels
   # subsample_com <- com[groups == groups_levels[i],]
    
    #Compute distance matrix for all communities
    #dist_com <- vegdist(subsample_com, method = method)
    
    #if(type == "centroid"){
    #  mean_treat_dist[i] <- median(dist_com)
    #}
    
    #if(type == "median"){
    #  mean_treat_dist[i] <- median(dist_com)
    #}
    
  #}
  
  #if(type == "centroid"){
  #  mean_within_dist <- median(mean_treat_dist)
  #}
  
  #if(type == "median"){
  #  mean_within_dist <- median(mean_treat_dist)
  #}
  
  
  #compute mean distance among groups
  
  #diff_dist <- rep(NA, iter)
  
  #for(p in 1:iter){
    
  #  among_coms <- matrix(NA, nrow = length(groups_levels), ncol = ncol(com))
  #  colnames(among_coms) <- colnames(com)
  #  among_coms <- data.frame(among_coms)
    
  #  for(i in 1:length(groups_levels)){
      
      #take communities from specific group levels
  #    subsample_com <- com[groups == groups_levels[i],]
      
  #    among_coms[i,] <- subsample_com[sample(1:nrow(subsample_com), 1),]
      
  #  }
    
  #  dist_among_com <- vegdist(among_coms, method = method)
    
  #  if(type == "centroid"){
  #    mean_among_dist <- median(dist_among_com)
  #  }
    
  #  if(type == "median"){
  #    mean_among_dist <- median(dist_among_com)
  #  }
    
  #  diff_dist[p] <- mean_among_dist - mean_within_dist
    
    
    
  #}
  
  #mean_diff_dist <- mean(diff_dist)
  #p_greater_dist <- length(diff_dist[diff_dist<=0])/iter
  #quantile(diff_dist, c(0.025, 0.975))
  #CI <- quantile(diff_dist, c(0.05, 1))
  #results_dist <- c(mean = mean_diff_dist, LCL = CI[1], UCL = CI[2], p_value = p_greater_dist)
  
  #return(list(distance = results_dist, centroid = results_cent))
  
  
}


