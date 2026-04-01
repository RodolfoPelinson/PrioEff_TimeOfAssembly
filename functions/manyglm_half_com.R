
manyglm_half_com <- function(comm,  comm_pred, covariables = NULL, n_rep = 10, n_cores = 1, multicore = TRUE){
  
require(mvabund)
require(doParallel)
require(pbapply)

  
colnames(comm) <- gsub(" ", "_", colnames(comm))  
colnames(comm_pred) <- gsub(" ", "_", colnames(comm_pred))  
cov <- covariables

list_com <- rep(list(comm_AM4_drop_atrasado), n_rep)

run_models <- function(y){
  
  
  n_species_used <- round(ncol(y)/2)
  
  first_half <- sample(1:ncol(y), n_species_used)
  first_half <- first_half[order(first_half)]
  #second_half <- c(1:ncol(y))[1:ncol(y) %in% first_half == FALSE]
  
  
  comm_first_half <- y[,first_half]
  #comm_second_half <- y[,second_half]
  
  comm_first_half_mv <- mvabund(comm_first_half)
  #comm_second_half_mv <- mvabund(comm_second_half)
  
  comm_pred_first_half <- comm_pred[,colnames(comm_pred) %in% colnames(comm_first_half_mv) == FALSE]
  #comm_pred_second_half <- comm_pred[,colnames(comm_pred) %in% colnames(comm_second_half_mv) == FALSE]
  
  if(is.null(cov) == FALSE){
    comm_pred_first_half <- data.frame(comm_pred_first_half, cov)
    #comm_pred_second_half <- data.frame(comm_pred_second_half, cov)
    formula_first_half_null <- formula( paste( "comm_first_half_mv ~", paste(colnames(cov), collapse = " + ") ) )
    #formula_second_half_null <- formula( paste( "comm_second_half_mv ~", paste(colnames(cov), collapse = " + ") ) )
  }else{
    formula_first_half_null <- formula( paste( "comm_first_half_mv ~ 1") )
    #formula_second_half_null <- formula( paste( "comm_second_half_mv ~ 1") )
  }
  

  
  
  formula_first_half <- formula( paste( "comm_first_half_mv ~", paste(colnames(comm_pred_first_half), collapse = " + ") ) )
  
  #formula_second_half <- formula( paste( "comm_second_half_mv ~", paste(colnames(comm_pred_second_half), collapse = " + ") ) )
  
  
  mod_first_half_null <- manyglm(formula_first_half_null,  family = "negative.binomial", data  = comm_pred_first_half)
  mod_first_half <- manyglm(formula_first_half,  family = "negative.binomial", data  = comm_pred_first_half)
  #mod_second_half_null <- manyglm(formula_second_half_null,  family = "negative.binomial", data  = comm_pred_second_half)
  #mod_second_half <- manyglm(formula_second_half,  family = "negative.binomial", data  = comm_pred_second_half)
  
  R2_first_half <- R2_manyglm(mod_first_half, mod_first_half_null)
  anova_first_falf <- anova(mod_first_half_null, mod_first_half, resamp = "montecarlo", show.time = "all", nBoot = 999)
  
  #R2_second_half <- R2_manyglm(mod_second_half, mod_second_half_null)
  #anova_second_falf <- anova(mod_second_half_null, mod_second_half, resamp = "montecarlo", show.time = "all", nBoot = 999)
  
  #R2 <- c(R2_first_half$com_R2$R2, R2_second_half$com_R2$R2)
  #R2_adj <- c(R2_first_half$com_R2$adj_R2, R2_second_half$com_R2$adj_R2)
  #p <- c(anova_first_falf$table$`Pr(>Dev)`[2], anova_second_falf$table$`Pr(>Dev)`[2])
  
  R2 <- c(R2_first_half$com_R2$R2)
  R2_adj <- c(R2_first_half$com_R2$adj_R2)
  p <- c(anova_first_falf$table$`Pr(>Dev)`[2])
  
  return(list(R2_adj = R2_adj,
              R2 = R2,
              p = p))
         
         
         
         
}

if(isTRUE(multicore)){
  message(paste("Using ", n_cores, " cores"))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  #.export = ls(globalenv())
  export_list <- c("mvabund", "manyglm", "anova.manyglm", "R2_manyglm")
  clusterExport(cl = cl, varlist = export_list)
  result_list <- pblapply(list_com, FUN = run_models, cl = cl)
  stopCluster(cl)
}else{
  result_list <- lapply(X = list_com, FUN = run_models)
}




R2_list <- list()
R2_adj_list <- list()
p_list <- list()

for(i in 1:length(result_list)){
  R2_list[[i]] <- result_list[[i]]$R2
  R2_adj_list[[i]] <- result_list[[i]]$R2_adj
  p_list[[i]] <- result_list[[i]]$p
}
  
R2_vec <- unlist(R2_list)
R2_adj_vec <- unlist(R2_adj_list)
p_vec <- unlist(p_list)


R2_adj_vec <- R2_adj_vec[R2_vec>=0]
p_vec <- p_vec[R2_vec>=0]
R2_vec <- R2_vec[R2_vec>=0]

p_value <- 1 - (length(p_vec[p_vec<=0.05]) / length(p_vec))
R2_mean <- mean(R2_vec)
R2_adj_mean <- mean(R2_adj_vec)

CI_R2 <- quantile(R2_vec, probs = c(0.975, 0.025))
CI_R2_adj <- quantile(R2_adj_vec, probs = c(0.975, 0.025))


result <- list(R2_mean = R2_mean,
               R2_adj_mean = R2_adj_mean,
               p_value = p_value,
               CI_R2 = CI_R2, 
               CI_R2_adj = CI_R2_adj,
               R2_vec = R2_vec,
               R2_adj_vec = R2_adj_vec,
               p_vec = p_vec)

return(result)

} 

