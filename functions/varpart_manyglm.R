
varpart_manyglm<- function(resp, pred, DF_adj_r2 = FALSE, ...){
  
  length_preds <- lapply(pred, FUN = ncol)
  
  pred <- pred[length_preds > 0]
  
  
  resp_mv <- mvabund(resp)
  
  if(length(pred) > 5 | length(pred) < 2){
    stop("variation partitioning can only be done for 2, 3, 4 or 5 predictor sets")
  }
  
  
  combinations <- list()
  for(i in 1:length(pred)){
    combinations[[i]] <- combn(1:length(pred),(1:length(pred))[i])
  }
  
  
  #preparing fractions 
  
  model_null <- manyglm(formula = resp_mv ~ 1, data = pred, ...)
  
  R2s_comb <- list()
  models_comb <- list()
  R2s_comb_sp <- list()
  
  #i se refere ao número de elementos na lista combinations
  for(i in 1:length(combinations)){
    
    R2s <- list()
    R2s_sp <- list()
    models <- list()
    names_z <- rep(NA, ncol(combinations[[i]]))
    
    #z se refere ao número de colunas em cada elemento da lista combinations
    for(z in 1:ncol(combinations[[i]])){
      
      
      #new_pred recebe os preditores presentes na primeira linha, coluna z
      new_pred <- pred[[   combinations[[i]][1,z]   ]]
      
      #Se o numero de linhas de cada coluna de i for maior que 1, pula essa etapa, se não adiciona mais preditores da linha j
      if(nrow(combinations[[i]])> 1){
        for(j in 2:nrow(combinations[[i]])){
          new_pred <- data.frame(new_pred, pred[[ combinations[[i]][j,z] ]] )
        }
      }
      
      
      right_fm <- paste(colnames(new_pred), collapse = " + ")
      formula <- formula(paste("resp_mv ~", right_fm))
      
      model <- manyglm(formula = formula, data = new_pred, family = "negative.binomial")
      
      R2 <- R2_manyglm(model, model_null)
      
      if(DF_adj_r2 == TRUE){
        R2s[[z]] <- R2$com_R2$adj_R2
        R2s_sp[[z]] <- R2$ind_R2$adj_R2
      }else{
        R2s[[z]] <- R2$com_R2$R2
        R2s_sp[[z]] <- R2$ind_R2$R2
      }
      models[[z]] <- model
      
      
      names_z[z] <- paste(combinations[[i]][,z], collapse = "-")
      
      
      
      
    }
    
    names(R2s) <- names_z
    names(models) <- names_z
    names(R2s_sp) <- names_z
    
    
    R2s_comb[[i]]<- R2s
    R2s_comb_sp[[i]]<- R2s_sp
    models_comb[[i]] <- models
    
    
  }
  
  
  R2 <- as.list(unlist(R2s_comb))
  
  R2_sp <- R2s_comb_sp[[1]]
  for(i in 2:length(R2s_comb_sp)){
    R2_sp <- c(R2_sp, R2s_comb_sp[[i]])
  }
  
  models <- models_comb[[1]]
  for(i in 2:length(models_comb)){
    models <- c(models, models_comb[[i]])
  }

  names_orig <- names(R2)
  names <- names_orig
  for(i in 1:length(pred)){
    names <- gsub(paste(i), names(pred)[i], names)
  }
  
  names(R2) <- names
  names(R2_sp) <- names
  names(models) <- names
  

  #Valores de R2 para todas as combinações de grupos de preditores.
  
  
  
  # Identificação das frações "puras"
  names_orig_new <- gsub("-", "",names_orig)
  only_mods <- nchar(names_orig_new) == length(pred)-1
  mods <- nchar(names_orig_new) == 1
  full_mods <- nchar(names_orig_new) == length(pred)
  
  R2_only <- unlist(R2)[full_mods] - unlist(R2)[only_mods]
  
  sp_names <- names(R2_sp[[1]])
  
  R2_sp_only <- list()
  for(i in 1:length(R2_sp[only_mods])){
    R2_sp_only[[i]] <- unlist(R2_sp[full_mods]) - unlist(R2_sp[only_mods][[i]])
    names(R2_sp_only[[i]]) <- sp_names
  }
  names(R2_sp_only) <- names(R2_sp[only_mods])
  

  names_only_mods <- rep(NA, length(pred))
  for(i in 1:length(pred)){
    names_only_mods[i] <- which(grepl(i, names_orig[only_mods]) == FALSE)
  }
  
  for(i in 1:length(pred)){
    names_only_mods <- gsub(paste(i), names(pred)[i], names_only_mods)
  }
  
  names(R2_only) <- names_only_mods
  names(R2_sp_only) <- names_only_mods
  
  
  R2_only <- R2_only[match(names(R2_only), names(R2)[c(1:length(pred))])]
  R2_sp_only <- R2_sp_only[match(names(R2_sp_only), names(R2_sp)[c(1:length(pred))])]
  
  
  R2_fractions_com <- data.frame(R2_full_fraction = unlist(R2)[c(1:length(pred))], R2_pure_fraction = R2_only)
  
  R2_fractions_sp <- list(R2_full_fraction = as.data.frame(R2_sp[c(1:length(pred))]), R2_pure_fraction = as.data.frame(R2_sp_only))
  
  R2_fractions_sp$R2_pure_fraction
  
  result <- list( R2_fractions_com = R2_fractions_com,
                  R2_fractions_sp = R2_fractions_sp,
                  #R2_models = unlist(R2),
                  R2_models = R2,
                  R2_models_sp = as.data.frame(R2_sp),
                  models = models,
                  model_null = model_null)
  
  return(result)
  
  
}
