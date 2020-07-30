#######################################################################
# Simulations
#  1) specify/read in tree structure
#  2) simulate Gaussian diffusion along the trees
#  3) simulate the data (multivariate binary at the tips)
#######################################################################


#' Simulate from latent class models with subject-specific indicators
#'
#'
#' @param n sample size
#' @param itemprob item probabilities
#' @param classprob class probabilities
#' @param fit fitted object
#'
#'
#' @return a list; \code{x}: data; \code{z} a vector of integer indicators
#' of class membership.
#'
#' @seealso \code{\link[BayesLCA]{blca}}
#'
#' @export
lotR_blca <- function (n, itemprob = 0.5, classprob = 1, fit = NULL)
{
  if (is.null(fit)) {
    itemprob <- as.matrix(itemprob)
    G <- nrow(itemprob)
    M <- ncol(itemprob)
  }
  else {
    itemprob <- fit$itemprob
    classprob <- fit$classprob
    G <- nrow(itemprob)
    M <- ncol(itemprob)
  }
  x <- matrix(runif(n * M), nrow = n)
  classvec <- as.vector(rmultinom(1, n, prob = classprob))
  ind <- c(0, cumsum(classvec))
  z <- rep(NA,n)
  for (g in 1:G) {
    x[(ind[g] + 1):ind[g + 1], ] <- t(t(x[(ind[g] +
                                                          1):ind[g + 1], ]) < itemprob[g, ]) * 1
    z[(ind[g] + 1):ind[g + 1]] <- g
    }
  make_list(x,z)
}
