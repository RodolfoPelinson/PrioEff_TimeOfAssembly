#' @title Null assembly
#'
#' @description This function is to run the null model used by "Siqueira et al. 2020"
#'
#' @export
#' @references #'Siqueira, Tadeu, Victor S. Saito, Luis M. Bini, Adriano S. Melo, Danielle K. Petsch, Victor L. Landeiro, Kimmo T. Tolonen, Jenny Jyrkänkallio‐Mikkola, Janne Soininen, and Jani Heino. “Community Size Can Affect the Signals of Ecological Drift and Niche Selection on Biodiversity.” Ecology 101, no. 6 (2020): e03014. \url{https://doi.org/10.1002/ecy.3014}.
#'


null_assembly <- function(test, species.col = "species", site.col = "site")
{
  a <- sample(test[ , which(colnames(test) == site.col)],
              length(test[ , which(colnames(test) == site.col)]))
  
  sim_table <- tapply(rep(1, length(a)),
                      list(a, as.vector(test[ , which(colnames(test) == species.col)])),
                      sum)
  
  sim_table[is.na(sim_table)] <- 0
  return(sim_table)
}
