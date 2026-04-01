get_scaled_lvs <- function(object, alpha = 0.5, ...){

  choose.lvs <- getLV(object, ...)
  choose.lv.coefs <- object$params$theta

  #Compute the singular-value decomposition of a rectangular matrix for ratation
  do_svd <- svd(object$lvs)

  svd_rotmat_sites <- do_svd$v
  svd_rotmat_species <- do_svd$v


  p <- NCOL(object$y)
  n <- NROW(object$y)
  num.lv <- object$num.lv

  idx <- matrix(TRUE, ncol = num.lv,  nrow = p)

  bothnorms <- vector("numeric", ncol(choose.lv.coefs))
  for (i in 1:ncol(choose.lv.coefs)) {
    bothnorms[i] <- sqrt(sum(choose.lvs[, i]^2)) * sqrt(sum(choose.lv.coefs[idx[, i], i]^2))
  }

  #Scale sites
  scaled_cw_sites <- t(t(choose.lvs)/sqrt(colSums(choose.lvs^2)) *       (bothnorms^alpha))

  #Scale species
  scaled_cw_species <- choose.lv.coefs
  for (i in 1:ncol(scaled_cw_species)) {
    scaled_cw_species[, i] <- choose.lv.coefs[, i]/sqrt(sum(choose.lv.coefs[idx[,
                                                                                i], i]^2)) * (bothnorms[i]^(1 - alpha))
  }


  choose.lvs <- scaled_cw_sites %*% svd_rotmat_sites
  choose.lv.coefs <- scaled_cw_species %*% svd_rotmat_species

  result <- list(sites = choose.lvs, species = choose.lv.coefs)

  return(result)
}
