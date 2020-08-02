#' Simulated multivariate binary data.
#'
#' Simulated data for testing the \code{\link{lcm_tree}} function.
#' Used in the example.
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{Y}{A matrix with n=1000 rows and J=18 columns.}
#'   \code{Z_obs}{A two column integer matrix. The first column is from
#'   1 to n. The second column is a mix of \code{NA} and integers between 1 and
#'   K (number of latent classess). NA represents unobserved latent class
#'   membership; an integer represent observed latent class.}
#'   \item{curr_leaves}{A character vector of length n=1000.
#'   \code{curr_leaves[i]} indicates which leaf does observation \code{i}
#'   belong to.}
#'   \item{truth}{a list: \code{tau}: the class probabilities (a single LCM);
#'   \code{theta}: the class-specific response probabilities;
#'   \code{prop_obs}: proportion of subjects with known classes;
#'   \code{Z} an integer vector of length n=1000, that indicates the true
#'   class memberships.}
#' }
"lotR_example_data"
