remove_sp <- function(com, n_sp = 0, INDEX = NULL, any = FALSE){
  
  if(is.null(INDEX)){
    com_oc <- decostand(com, method = "pa")
    com <- com[,colSums(na.omit(com_oc)) > n_sp]
    return(com)
  }else{
    INDEX <- as.factor(INDEX)
    sps <- matrix(NA, nrow = length(levels(INDEX)), ncol = ncol(com))
    keep <- rep(NA, ncol(com))  
    
    for(i in 1:length(levels(INDEX))){
      splited <- com[INDEX == levels(INDEX)[i],]
      com_oc <- decostand(splited, method = "pa")
      sps[i,] <- colSums(na.omit(com_oc)) > n_sp
    }
    
    for(i in 1:ncol(sps)){
      if(any == TRUE){
        keep[i] <- any(sps[,i] == TRUE)
      }else{
        keep[i] <- any(sps[,i] == FALSE) == FALSE
      }
    }
    
    com <- com[,keep]
    return(com)
    
    
  }
  
  


}
