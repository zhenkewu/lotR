#' Simulated multivariate binary data.
#'
#' Simulated data for testing the [lcm_tree()] function.
#' Used in the example.
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{Y}{A matrix with n=1000 rows and J=18 columns.}
#'   `Z_obs`{A two column integer matrix. The first column is from
#'   1 to n. The second column is a mix of `NA` and integers between 1 and
#'   K (number of latent classess). NA represents unobserved latent class
#'   membership; an integer represent observed latent class. `NULL` if
#'   no observation has observed class membership.}
#'   \item{curr_leaves}{A character vector of length n=1000.
#'   `curr_leaves[i]` indicates which leaf does observation `i`
#'   belong to.}
#'   \item{truth}{a list: `tau`: the class probabilities (a single LCM);
#'   `theta`: the class-specific response probabilities;
#'   `prop_obs`: proportion of subjects with known classes;
#'   `Z` an integer vector of length n=1000, that indicates the true
#'   class memberships.}
#' }
"lotR_example_data"
