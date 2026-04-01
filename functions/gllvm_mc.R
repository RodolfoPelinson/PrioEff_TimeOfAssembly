gllvm_mc <- function(y, X = NULL, seed = NULL, formula = NULL, family,
                     num.lv = 1, method = "VA", starting.val = "res",
                     n_runs = 50, n_cores = 1, export_list, row.eff = FALSE, ...){
  require(gllvm)
  require(doParallel)
  require(pbapply)
  
  if(is.null(seed) == FALSE){  set.seed(seed)}
  
  
  coms <- list()

  for(i in 1:n_runs){
    coms[[i]] <- y
  }
  
  
  my_gllvm <- function(y){
    mod <- gllvm(y,X = X, formula = formula,  family = family,
          num.lv = num.lv, method = method, control.start  = list(n.init = 1, n.init.max = 1, starting.val = starting.val, jitter.var = 10), row.eff = row.eff, ...)
  
  return(mod)  
  }

  message(paste("Using ", n_cores, " cores"))
  #save.image()
    
  cl <- makeCluster(n_cores, type="PSOCK")
  #registerDoParallel(cl)
  
  #clusterExport(cl = cl, varlist = list("my_gllvm","gllvm", "method" ,"num.lv", "family", "starting.val", "X", "coms", "Exp_design"))
  
  
  clusterExport(cl = cl, varlist = export_list)
  
  
  mods <- pblapply(coms, FUN = my_gllvm, cl = cl)
  
  #mods <- lapply(coms, FUN = my_gllvm)
  

  stopCluster(cl)

  AICcS <- unlist(lapply(mods, AICc.gllvm))
  
  ord_AICc<-order(AICcS)
  
  AICcS_ordered <- AICcS[ord_AICc]
  
  q0 <- 1
  second_best <- 2
  third_best <- 3
  
  q50 <- round(n_runs / 2)
  q5 <- round(n_runs / 20)
  q1 <- round(n_runs / 100)
  q10 <- round(n_runs / 10)
  
  if(q1== 0){q1 <- 1}
  if(q5== 0){q5 <- 1}
  if(q10== 0){q10 <- 1}
  if(q50== 0){q50 <- 1}
  
  
  second_best_2<-match(AICcS_ordered[second_best], AICcS)
  third_best_2<-match(AICcS_ordered[third_best], AICcS)
  
  q0_2<-match(AICcS_ordered[q0], AICcS)
  q1_2<-match(AICcS_ordered[q1], AICcS)
  q5_2<-match(AICcS_ordered[q5], AICcS)
  q10_2<-match(AICcS_ordered[q10], AICcS)
  q50_2<-match(AICcS_ordered[q50], AICcS)
  
  
  second_best <- mods[[second_best_2]]
  third_best <- mods[[third_best_2]]
  quantile_0 <- mods[[q0_2]]
  quantile_1 <- mods[[q1_2]]
  quantile_5 <- mods[[q5_2]]
  quantile_10 <- mods[[q10_2]]
  quantile_50 <- mods[[q50_2]]
  
  #quantile_50 <- mods[[q50[1]]]
  
  #plot(AICcS_ordered ~ 1:length(n_runs), type = "lines", ylab = "AICc", "Model")
  #abline(v = AICcS_ordered[1], col = "red")
  #abline(v = AICcS_ordered[q1], col = "red")
  #abline(v = AICcS_ordered[q5], col = "red")
  #abline(v = AICcS_ordered[q10], col = "red")
  #abline(v = AICcS_ordered[q50], col = "red")
  
  
  result <- list(best = quantile_0,
                 second_best = second_best,
                 third_best = third_best,
                 quantile_1 = quantile_1,
                 quantile_5 = quantile_5,
                 quantile_10 = quantile_10,
                 quantile_50 = quantile_50,
                 AICc = AICcS_ordered)
  #result <- AICcS
  rm(mods)
  

return(result)  
  
}


