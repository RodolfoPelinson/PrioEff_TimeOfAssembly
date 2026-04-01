
pairwise_gllvm <- function(com, pairwise, num.lv = 0, covariable = NULL, by = NULL, family = "negative.binomial", p_method = "fdr", remove_species = 1, method = "EVA", seed = 1, starting.val  = "res", n_runs = 1, multicore = FALSE, n.init = 20, n.init.max = 10, row.eff = FALSE,...) {
  
  #############################################
  
  # Se quiser fazer paralelo, tem que programar especificamente
  ####################
  
  
  pairwise <- as.factor(as.character(pairwise))
  
  if(is.null(covariable) == FALSE){
    covariable <- as.factor(as.character(covariable))
  }
  
  ######################
  
  if(is.null(by)){
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
      
      
      
      if(isTRUE(multicore)){
        if(is.null(covariable)){
          
          X <-data.frame(pairwise_subset)
          
          fit00 <- gllvm_mc(y = com_subset, family = family, num.lv = num.lv, method = method, starting.val = starting.val, seed = seed, n_runs = n_runs, export_list = list("gllvm"), n_cores = detectCores()/2)
          fit01 <- gllvm_mc(y = com_subset, formula = ~ pairwise_subset, X = X, family = family, num.lv = num.lv, method = method, starting.val = starting.val, seed = seed, n_runs = n_runs, export_list = list("gllvm"), n_cores = detectCores()/2)
          
          #if(mean(fit1$sd$Xcoef[,1]) < 0.1){
          
          #  if(method == "VA"){method <- "EVA"}else{method <- "VA"}
          #  fit0 <- gllvm(y = com_subset, family = family, num.lv = num.lv, method = "EVA", control.start  = list(method = method, ...), seed = seed)
          #  fit1 <- gllvm(y = com_subset, formula = ~ pairwise_subset, X = data.frame(pairwise_subset), family = family, num.lv = num.lv, method = "EVA", control.start  = list(method = method, ...), seed = seed)
          
          #  if(method == "EVA"){method <- "VA"}else{method <- "EVA"}
          #}
          
          fit0<- fit00$best
          fit1<- fit01$best
          
          rm(fit00)
          rm(fit01)
          
          
          models_null[[i]] <- fit0
          models_predictor[[i]] <- fit1
          com_y[[i]] <- com_subset
          predictors[[i]] <- data.frame(treatments = pairwise_subset)
          
        }else{
          
          covariable_subset <- covariable[pairwise == comparisons[1,i] | pairwise == comparisons[2,i]]
          
          X <- data.frame(covariable_subset, pairwise_subset)
          
          fit00 <- gllvm_mc(y = com_subset, formula = ~ covariable_subset, X = X, family = family, num.lv = num.lv, method = method, starting.val = starting.val, seed = seed, n_runs = n_runs, export_list = list("gllvm"), n_cores = detectCores()/2)
          fit01 <- gllvm_mc(y = com_subset, formula = ~ pairwise_subset + covariable_subset, X = X, family = family, num.lv = num.lv, method = method, starting.val = starting.val, seed = seed, n_runs = n_runs, export_list = list("gllvm"), n_cores = detectCores()/2)
          
          #if(mean(fit1$sd$Xcoef[,1]) < 0.1){
          #  
          #  if(method == "VA"){method <- "EVA"}else{method <- "VA"}
          #  fit0 <- gllvm(y = com_subset, formula = ~ covariable_subset, X = data.frame(covariable_subset), family = family, num.lv = num.lv, method = method, control.start  = list(...), seed = seed)
          #  fit1 <- gllvm(y = com_subset, formula = ~ pairwise_subset + covariable_subset, X = data.frame(covariable_subset, pairwise_subset), family = family, num.lv = num.lv, method = method, control.start  = list(...), seed = seed)
          #  if(method == "EVA"){method <- "VA"}else{method <- "EVA"}
          #}else
          
          fit0 <- fit00$best
          fit1 <- fit01$best
          
          rm(fit00)
          rm(fit01)
          
          
          models_null[[i]] <- fit0
          models_predictor[[i]] <- fit1
          com_y[[i]] <- com_subset
          predictors[[i]] <- data.frame(block2 = covariable_subset, treatments = pairwise_subset)
          
        }
      }else{
        if(is.null(covariable)){
          fit0 <- gllvm(y = com_subset, family = family, num.lv = num.lv, method = method, control.start  = list(starting.val = starting.val,n.init = n.init, n.init.max = n.init.max,  ...), seed = seed, row.eff = row.eff)
          fit1 <- gllvm(y = com_subset, formula = ~ pairwise_subset, X = data.frame(pairwise_subset), family = family, num.lv = num.lv, method = method, control.start  = list(starting.val = starting.val, n.init = n.init, n.init.max = n.init.max, ...), seed = seed, row.eff = row.eff)
          
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
          covariable_subset <- covariable[pairwise == comparisons[1,i] | pairwise == comparisons[2,i]]
          fit0 <- gllvm(y = com_subset, formula = ~ covariable_subset, X = data.frame(covariable_subset), family = family, num.lv = num.lv, method = method, control.start  = list(starting.val = starting.val,n.init = n.init, n.init.max = n.init.max,...), seed = seed, row.eff = row.eff)
          fit1 <- gllvm(y = com_subset, formula = ~ pairwise_subset + covariable_subset, X = data.frame(covariable_subset, pairwise_subset), family = family, num.lv = num.lv, method = method, control.start  = list(starting.val = starting.val,n.init = n.init, n.init.max = n.init.max, ...), seed = seed, row.eff = row.eff)
          
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
      }
      
      

      gc()
      
      anova <- anova(fit0, fit1)
      

      Dev <- as.numeric(anova$D[2])
      p <- as.numeric(anova$P.value[2])
      
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
    
    return(results_list)
  }
  
  
  
  
  
  
  
  
 # if(is.null(by) == FALSE){
    
#    if(is.null(covariable) == FALSE){
 #     if(by == covariable){
  #      stop("covariable can be the same as the by argument")
  #    }
  #  }
    
 #   levels_by <- levels(by)
    
#    for(i in 1:length(levels_by)){
#      com_subset0 <- com[by == levels_by[i],]
#      pairwise_subset0 <- as.factor(as.character(pairwise[by == levels_by[i]]))
#      
#      if(is.null(covariable) == FALSE){
#        covariable_subset0 <- as.factor(as.character(covariable[by == levels_by[i]]))
#      }
#      
#     levels <- levels(pairwise_subset0)
#      n_levels <- length(levels)
      
#      n_comparisons <- choose(length(levels), 2)
      
#      comparisons <- matrix(NA, nrow = 2, ncol = n_comparisons)
      
#      results <- matrix(NA, nrow = n_comparisons, ncol = 3)
#      results_rows <- rep(NA, n_comparisons)
      
      ##############################SANTO GRAAAU
#      k <- 0
#      for (t in 1:(n_levels - 1) ) {
#        for (j in (t + 1):n_levels) {
#          k=k+1
#          results_rows[k] <- paste(levels_by[i], "-" ,levels[t], 'vs', levels[j])
#          comparisons[,k] <- c(levels[t], levels[j])
#        }
#      }
      
#      rownames(results) <- results_rows
#      colnames(results) <- c("Dev", "p", "adjusted p")
      
      
      #####################################
      
#      for(j in 1:n_comparisons){
#        com_subset <- com_subset0[pairwise_subset0 == comparisons[1,j] | pairwise_subset0 == comparisons[2,j],]
#        pairwise_subset <- as.factor(as.character(pairwise_subset0[pairwise_subset0 == comparisons[1,j] | pairwise_subset0 == comparisons[2,j]]))
        
#        if(is.null(covariable) == FALSE){
#          covariable_subset <- covariable_subset0[pairwise_subset0 == comparisons[1,j] | pairwise_subset0 == comparisons[2,j]]
#          fi0 <- manyglm(com_subset ~  covariable_subset, family = family)
#          fi1 <- manyglm(com_subset ~  covariable_subset + pairwise_subset, family = family)
#        }
#        if(is.null(covariable)){
#          fi0 <- manyglm(com_subset ~  1, family = family)
#          fi1 <- manyglm(com_subset ~  pairwise_subset, family = family)
#        }
        
#        anova <- anova.manyglm(fi0, fi1, test = test, resamp = resamp, nBoot = nBoot)
        
#        Dev <- anova$table$Dev[2]
#        p <- anova$table$`Pr(>Dev)`[2]
        
#        results[j,1] <- Dev
#        results[j,2] <- p
        
        
        
#      }
      
#      if(i==1){
#        results_previous <- results
#      }else{results_previous <- rbind(results_previous, results)}
      
#    }
#    p_adjusted <- p.adjust(results_previous[,2], method = p_method)
#    results_previous[,3] <- p_adjusted
    
#    return(results_previous)
#  }
  
}

