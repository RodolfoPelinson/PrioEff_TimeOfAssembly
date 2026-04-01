
extract_mles <- function(mod_gllvm){
  mles <- mod_gllvm$params$Xcoef[,1]
  se_mles <- mod_gllvm$sd$Xcoef[,1]
  
  z <- qnorm(0.975)
  
  upper_CI <- mles + (se_mles * z)
  lower_CI <- mles - (se_mles * z)
  
  mles_frame <- data.frame(mles, upper_CI, lower_CI)
  
  return(mles_frame)
}



