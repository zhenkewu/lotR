#######################################################################
# Simulations
#  1) specify/read in tree structure
#  2) simulate Gaussian diffusion along the trees
#  3) simulate the data (multivariate binary at the tips)
#######################################################################

#' Simulate data and subject-specific indicators from latent class models
#'
#'
#' @param n sample size
#' @param itemprob item probabilities
#' @param classprob class probabilities
#' @param fit fitted object
#'
#'
#' @return a list; `x`: data; `z` a vector of integer indicators
#' of class membership.
#'
#' @seealso [BayesLCA::blca()]
#' @importFrom stats rmultinom runif
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


#' Simulate data and subject-specific indicators from tree-structured latent class models
#'
#' The observations belong to leaves that may further belong to a few groups, each
#' with its own K-class probabilities. We assume all the leaves share the same
#' set of K class-specific response probability profiles.
#'
#' @param n sample size
#' @param itemprob item probabilities; this is shared across leaf nodes; K by J
#' @param mytree see [design_tree()]
#' @param pi_mat class probabilities for  pL leaf nodes; it is pL by K.
#' @param h_pau a p-dim vector of positive values that indicates the branch lengths
#' between a node `u` and its parent `pa(u)`
#'
#' @return a list
#' \describe{
#' \item{Y}{observations leaf by leaf}
#' \item{curr_leaves}{leaf names, need to be for each row of Y}
#' \item{truth}{a list that contains the simulation truth:
#' \itemize{
#' \item `Z` true class indicators for all observations
#' \item `itemprob` a K by J matrix of response probability profiles
#' \item `pi_mat` the eta_v transformed to pi_v; pL by K
#' \item `h_pau` a vector of p values, each representing the
#' branch length between the node `u` and its parent node `pa(u)`
#' }
#' }
#' }
#'
#' @examples
#'
#' library(igraph)
#' n = 1000
#' tau   <- c(0.6,0.3,0.1)
#' itemprob <- rbind(rep(rep(c(0.9, 0.9), each = 1),9),
#'                   rep(rep(c(0.5, 0.5), each = 1),9),
#'                   rep(rep(c(0.1, 0.1), each = 1),9))
#'
#' data("lotR_example_edges")
#' mytree <- igraph::graph_from_edgelist(lotR_example_edges, directed = TRUE)
#' # Plot tree
#'
#' nodes  <- names(igraph::V(mytree))
#' leaves <- names(igraph::V(mytree)[degree(mytree, mode = "out") == 0])
#' pL = length(leaves)
#' p  = length(igraph::V(mytree))
#'
#' ######
#' K = nrow(itemprob)
#' # specify the nodes that have non-trivial alpha_u, this was called
#' # xi_u, because xi_u = s_u*alpha_u and s_u = 1 if we set it in the simulation.
#' alpha_mat = rbind(logit(prob2stick(tau)[-K]),
#'                   c(-1,-0.5),
#'                   c(1,0.5),
#'                   matrix(0,nrow=p-3,ncol=K-1)
#' )
#'
#' # get lists of ancestors for each leaf_ids:
#' d <- igraph::diameter(mytree,weights=NA)
#' # need to set weight=NA to prevent the use of edge lengths in determining the diameter.
#' ancestors <- igraph::ego(mytree, order = d + 1, nodes = leaves, mode = "in")
#' ancestors <- sapply(ancestors, names, simplify = FALSE)
#' ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a),
#'                     nodes = nodes, simplify = FALSE)
#' names(ancestors) <- leaves
#'
#' # calculate the class probabilities for all leaf nodes; each leaf node
#' # should have a K-dim vector that sums to one; Some nodes may share
#' # the same set of K-dim probability vector, others may differ. There are
#' # one or more groups of leaf nodes with distinct K-dim probability vectors.
#' # Note the branch lengths may also be used here.
#' pi_mat <- matrix(NA,nrow=pL,ncol=K)
#' for (v in seq_along(leaves)){
#'   pi_mat[v,-K] <- colSums(alpha_mat[ancestors[[v]],,drop=FALSE])
#'   pi_mat[v,] <- tsb(c(expit(pi_mat[v,-K]),1))
#' }
#'
#' # s = c(1, 1,1,0,0, rep(0,pL)) # effective nodes
#' h_pau = rep(1,p)
#'
#' lotR_example_data_tree <- simulate_lcm_tree(n,itemprob,mytree,pi_mat,h_pau)
#' #save the simulated data to the R package for illustration:
#' # save(lotR_example_data_tree, file = "data/lotR_example_data_tree.rda", compress = "xz")
#' @seealso [BayesLCA::blca()]
#' @importFrom stats rmultinom runif
#' @importFrom igraph V
#' @export
simulate_lcm_tree <- function (n, itemprob, mytree, pi_mat, h_pau)
{
  # a few quick calculations:
  K <- nrow(itemprob)
  nodes  <- names(V(mytree))
  leaves <- names(V(mytree)[degree(mytree, mode = "out") == 0])
  pL <- length(leaves)
  p  <- length(V(mytree))

  # simulate number of observations for each leaf node:
  N_sim <- 1:pL
  N_sim <- c(N_sim,sample(1:pL,size=n-pL,prob=rep(1/pL,pL),replace=TRUE)) # even leafs.
  N_sim <- as.integer(table(sort(N_sim)))

  # simulate the observations leaf by leaf: Y, curr_leaves, truth (
  # - true set of nodes with non-trivial increment xi_u, i.e.,
  # alpha_u, the actual values of these increments, the eta_v which
  # sums over all non-trivial ancestral nodes' value of alpha_u;
  # - the eta_v transformed to pi_v;
  # - the class-specific response probabilities)

  # simulate the multivariate responses
  Y_sim <- list()
  # simulate the class memberships
  Z_sim <- list()
  curr_leaves_sim <- list()
  for (v in 1:pL){
    simu <- lotR_blca(N_sim[v], itemprob = itemprob, classprob = pi_mat[v,])
    Y_sim[[v]] <- simu$x
    curr_leaves_sim[[v]] <- rep(v,N_sim[v])
    Z_sim[[v]] <- simu$z
  }

  Y <- do.call("rbind",Y_sim)
  curr_leaves <- leaves[do.call("c",curr_leaves_sim)]
  Z <- do.call("c",Z_sim)

  truth <- make_list(Z,itemprob, mytree, pi_mat,h_pau)
  make_list(Y,curr_leaves,truth)
}








