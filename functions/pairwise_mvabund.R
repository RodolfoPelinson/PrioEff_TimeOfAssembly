
pairwise_mvabund <- function(com, pairwise, covariable = NULL, by = NULL, family = "negative.binomial", p_method = "fdr", remove_species = NULL, seed = 1, composition = FALSE,...) {
  
  #############################################
  
  # Se quiser fazer paralelo, tem que programar especificamente
  ####################
  
  
  pairwise <- as.factor(as.character(pairwise))
  
  if(is.null(covariable) == FALSE){
    covariable <- as.factor(as.character(covariable))
  }
  

  #if(is.null(by)){
    levels <- levels(pairwise)
    n_levels <- length(levels)
    
    n_comparisons <- choose(length(levels), 2)
    
    comparisons <- matrix(NA, nrow = 2, ncol = n_comparisons)
    
    results <- matrix(NA, nrow = n_comparisons, ncol = 3)
    results_rows <- rep(NA, n_comparisons)
    
    models_null <- list()
    models_predictor <- list()
    com_y <- list()
    predictors <- list()
    
    ##############################SANTO GRAAAU
    k <- 0
    for (i in 1:(n_levels - 1) ) {
      for (j in (i + 1):n_levels) {
        k=k+1
        results_rows[k] <- paste(levels[i], 'vs', levels[j])
        comparisons[,k] <- c(levels[i], levels[j])
      }
    }
    
    rownames(results) <- results_rows
    colnames(results) <- c("Dev", "p", "adjusted p")
    
    
    #####################################
    
    for(i in 1:n_comparisons){
      com_subset <- com[pairwise == comparisons[1,i] | pairwise == comparisons[2,i],]
      
      if(is.null(remove_species)==FALSE){
        com_subset <- com_subset[,colSums(decostand(com_subset, method = "pa")) > remove_species]
      }else{
        com_subset <- com_subset[,colSums(decostand(com_subset, method = "pa")) > 0]
      }
      
      pairwise_subset <- as.factor(as.character(pairwise[pairwise == comparisons[1,i] | pairwise == comparisons[2,i]]))
      
      
      

        if(is.null(covariable)){
          
          com_subset <- mvabund(com_subset)
          data_pairwise <- data.frame(pred = pairwise_subset)
          fit0 <- manyglm(com_subset ~ 1, family="negative.binomial", composition = composition, data = data_pairwise)
          fit1 <- manyglm(com_subset ~ pred, family="negative.binomial", composition = composition, data = data_pairwise)
          
          #if(mean(fit1$sd$Xcoef[,1]) < 0.1){
          
          #  if(method == "VA"){method <- "EVA"}else{method <- "VA"}
          #  fit0 <- gllvm(y = com_subset, family = family, num.lv = num.lv, method = "EVA", control.start  = list(method = method, ...), seed = seed)
          #  fit1 <- gllvm(y = com_subset, formula = ~ pairwise_subset, X = data.frame(pairwise_subset), family = family, num.lv = num.lv, method = "EVA", control.start  = list(method = method, ...), seed = seed)
          
          #  if(method == "EVA"){method <- "VA"}else{method <- "EVA"}
          #}
          
          models_null[[i]] <- fit0
          models_predictor[[i]] <- fit1
          com_y[[i]] <- com_subset
          predictors[[i]] <- data.frame(treatments = pairwise_subset)
          
        }else{
          com_subset <- mvabund(com_subset)
          
          covariable_subset <- covariable[pairwise == comparisons[1,i] | pairwise == comparisons[2,i]]
          fit0 <- manyglm(com_subset ~ covariable_subset, family="negative.binomial", composition = composition, data = data.frame(covariable_subset, pairwise_subset))
          fit1 <- manyglm(com_subset ~ pairwise_subset + covariable_subset, family="negative.binomial", composition = composition, data = data.frame(covariable_subset, pairwise_subset))
          
          #if(mean(fit1$sd$Xcoef[,1]) < 0.1){
          #  
          #  if(method == "VA"){method <- "EVA"}else{method <- "VA"}
          #  fit0 <- gllvm(y = com_subset, formula = ~ covariable_subset, X = data.frame(covariable_subset), family = family, num.lv = num.lv, method = method, control.start  = list(...), seed = seed)
          #  fit1 <- gllvm(y = com_subset, formula = ~ pairwise_subset + covariable_subset, X = data.frame(covariable_subset, pairwise_subset), family = family, num.lv = num.lv, method = method, control.start  = list(...), seed = seed)
          #  if(method == "EVA"){method <- "VA"}else{method <- "EVA"}
          #}else
          
          models_null[[i]] <- fit0
          models_predictor[[i]] <- fit1
          com_y[[i]] <- com_subset
          predictors[[i]] <- data.frame(block2 = covariable_subset, treatments = pairwise_subset)
          
        }
      
      
      
      

      anova <- anova(fit0, fit1, ...)
      
      
      Dev <- as.numeric(anova$table$Dev[2])
      p <- as.numeric(anova$table$`Pr(>Dev)`[2])
      
      results[i,1] <- Dev
      results[i,2] <- p
      
      message(paste((i/n_comparisons)*100),"%")
    }
    
    p_ajusted <- as.numeric(p.adjust(results[,2], method = p_method))
    
    sig <- rep("", nrow(results))
    
    for(i in 1:length(sig)){
      if(p_ajusted[i] <= 0.1){sig[i] <- "."}
      if(p_ajusted[i] <= 0.05){sig[i] <- "*"}
      if(p_ajusted[i] <= 0.01){sig[i] <- "**"}
      if(p_ajusted[i] <= 0.005){sig[i] <- "***"}
    }
    
    results[,3] <- p_ajusted
    
    results[,1] <- round(results[,1], 3)
    results[,2] <- round(results[,2], 5)
    results[,3] <- round(results[,3], 5)
    
    results <- data.frame(results, sig=sig)
    
    names(models_null) <- results_rows
    names(models_predictor) <- results_rows
    names(com_y) <- results_rows
    names(predictors) <- results_rows
    
    
    results_list <- list(models_null = models_null,
                         models_predictor = models_predictor,
                         results = results,
                         com_y = com_y,
                         predictors = predictors)
    
 # }
  
  return(results_list)
  
  
  
}

