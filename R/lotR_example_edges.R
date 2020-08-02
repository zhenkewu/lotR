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
