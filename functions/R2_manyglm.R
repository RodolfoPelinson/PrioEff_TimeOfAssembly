R2_manyglm <- function(model, model_null){
  n<-nrow(model$y)
  
  #R2 <- c(1 - ( exp(logLik(model_null)-logLik(model))^(2/n)) )
  
  R2_c <- c(1 - ( exp(logLik(model_null)-logLik(model))^(2/n)) ) / (1 - exp(logLik(model_null))^(2/n) )
  
  names(R2_c) <- colnames(model$y)
  
  p <- nrow(data.frame(model$coefficients)) - 1
  
  #adj_R2 <- 1 - (1 - R2) * (n-1)/(n-p-1)
  
  adj_R2 <- R2_c * ((n-p)/(n-1))
  
  ind_R2 <- list(R2 = R2_c, adj_R2 = adj_R2)
  
  com_R2 <- list(R2 = mean(R2_c), adj_R2 = mean(adj_R2))
  
  return(list(ind_R2 = ind_R2, com_R2 = com_R2))
}

