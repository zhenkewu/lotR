## This file documents all example data.

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
#'   K (number of latent classes). `NA` represents unobserved latent class
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


#' Edgelist for constructing an example tree.
#'
#' A matrix representing edges
#' from the outcome tree formed by cardiovascular diseases
#' in category 7.4
#' (diseases of arteries, arterioles, and capillaries)
#' of the multilevel Clinical Classifications Software (CCS)
#' hierarchical disease classification system.
#' See vignette("moretrees") for details of how to construct
#' a tree from this edgelist.
#'
#' @format The edge list is a matrix with 15 rows with 2 columns:
#' the left column
#' represents parent nodes or categories, and the right column
#' represents children or subcategories of the parents.
#' There is one row for every parent-child pair.
"lotR_example_edges"


#' Simulated multivariate binary data. (with tree structure)
#'
#' Simulated data for testing the [lcm_tree()] function.
#' Used in the example. (K=3)
#'
#' @format A named list with the following elements:
#' see the returned values of [simulate_lcm_tree()]
#' @seealso [simulate_lcm_tree()]
"lotR_example_data_tree"


#' Simulated multivariate binary data. (with tree structure)
#'
#' Simulated data for testing the [lcm_tree()] function.
#' Used in the example. (K=2)
#'
#' @format A named list with the following elements:
#' see the returned values of [simulate_lcm_tree()]
#' @seealso [simulate_lcm_tree()]
"lotR_example_data_tree_K2"

