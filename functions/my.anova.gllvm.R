my.anova.gllvm <- function(object, ...) {
  objects <- list(object, ...)
  
  if (length(objects) < 2)
    stop("At least two objects are needed for tests.")
  if (any(!(sapply(objects, function(x)inherits(x,"gllvm")))))
    stop("The function 'anova.gllvm' can only be used for a gllvm object.")
  
  tt <- sapply(objects, function(x)
    x$method)
  if (!(all(tt == "VA") | all(tt == "LA") | all(tt == "EVA")))
    stop("The objects are not comparable when they are fitted using different methods.")
  
  y <- object$y
  n <- NROW(y)
  p <- NCOL(y)
  diff <- sapply(objects, function(x) sum(x$y - y))
  if (any(!(diff == 0)))
    stop("The objects can not be compared")
  
  df.list <- sapply(objects, function(x)
    attr(logLik.gllvm(x), "df"))
  objects_order <- objects[order(df.list)]
  formulas <- sapply(objects_order, function(x)
    formula(x$terms))
  
  #if(test=="LR"){
  df.list <- sapply(objects_order, function(x) attr(logLik.gllvm(x), "df"))
  ll.list <- sapply(objects_order, logLik.gllvm)
  
  D <- 2 * (ll.list[-1] - ll.list[1:(length(df.list) - 1)])
  df.chisq <- (df.list[-1] - df.list[1:(length(df.list) - 1)])
  Pval <- 1 - pchisq(D, df.chisq)
  paste("Model", 1:length(objects_order))
  if(Pval <= 0.01){sig <- "*"}else{sig<-""}
  result <- data.frame( Resid.Df = n * p - df.list, D = c(0, D), Df.diff = c(0, df.chisq), P.value = c("", signif(Pval)),`sig?`= c("", sig))
  #}
  if (any(result$Df > 20))
    warning( "This test was not designed for tests with a df.diff larger than 20 so the P-value should be treated as approximate.\n")
  for (i in 1:length(objects_order)) {
    formchar <- as.character(formulas[[i]])
    if (length(formchar) == 3)
      formchar <- formchar[c(2, 1, 3)]
    cat("Model ", i, ": ", formchar, "\n")
  }
  return(result)
  
}
