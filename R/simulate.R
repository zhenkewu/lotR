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
#' @examples
#'
#' lotR_blca(10,matrix(0.5,nrow=3,ncol=20),c(0.98,0.01,0.01))
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
    if (ind[g]<ind[g+1]){
    x[(ind[g] + 1):ind[g + 1], ] <- t(t(x[(ind[g] +
                                             1):ind[g + 1], ]) < itemprob[g, ]) * 1
    z[(ind[g] + 1):ind[g + 1]] <- g
    }
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
#' @param balanced by default is `TRUE` to uniformly assign observations to the leaf nodes;
#' otherwise set this to `FALSE`.
#' @param ratio for a pair of leaves; the sample ratios (larger vs smaller ones);
#' in the event of an odd number of leaves; the smaller leaf in the pair is kept.
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
#' # save(lotR_example_data_tree, file = "data/lotR_example_data_tree2.rda", compress = "xz")
#' @seealso [BayesLCA::blca()]
#' @importFrom stats rmultinom runif
#' @importFrom igraph V
#' @export
simulate_lcm_tree <- function (n, itemprob, mytree, pi_mat, h_pau,balanced=TRUE,ratio=4)
{
  # a few quick calculations:
  K <- nrow(itemprob)
  nodes  <- names(V(mytree))
  leaves <- names(V(mytree)[degree(mytree, mode = "out") == 0])
  pL <- length(leaves)
  p  <- length(V(mytree))

  # simulate number of observations for each leaf node:
  N_sim <- rep(1:pL,each=2) #at least two observations per leaf node.
  prob_vec <- rep(1/pL,pL)
  if (!balanced){
      prob_vec <- rep(c(1,ratio), c(floor(pL/2),pL-floor(pL/2))) # currently 1:4 ratio.
      prob_vec <- sample(prob_vec/sum(prob_vec))
    }
  N_sim <- c(N_sim,sample(1:pL,size=n-2*pL,prob=prob_vec,replace=TRUE)) # even leafs.
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



#' generate all permutations of 1 to n
#'
#' Copied from Museful's solution from here
#' <https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r>
#' This also mimics `RcppAlgos::permuteGeneral(n,n)`
#'
#' @param n an integer
#'
#' @examples
#' permutations(3)
#'
#' @returns a matrix of n columns and `factorial(n)` rows
#' @export
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}


#' Get optimal column permutation of a matrix B to match the matrix A
#'
#' Useful in simulations where the posterior sampling relabels the classes,
#' which is equivalent and cannot be told apart by the likelihood only.
#'
#' @param A a matrix (e.g., truth)
#' @param B another matrix (e.g., estimated); A and B must have the same dimensions;
#' B's columns may be permutated to best match those of A.
#'
#' @return a permutation of columns of B (represented by a permuted `1:ncol(A)`)
#'
#' @examples
#'
#' A <- matrix(c(1,2,3,4,5,6),nrow=2,ncol=3)
#' B <- A[,c(2,1,3)]
#' opt_colpermB(A,B) # should expect c(2,1,3)
#' @export
opt_colpermB <- function(A,B){
  if (dim(A)[1]!=dim(B)[1] && dim(A)[2]!=dim(B)[2]){
    stop("[lotR] matrix of different dimensions")}
  n <- ncol(A)
  col_perms  <- permutations(n)
  dist_mat   <- apply(col_perms, 1, function(col) sum(abs(A- B[,col,drop=FALSE])^2))
  optim_cols <- which.min(dist_mat)
  col_perms[optim_cols, ]
}

#' calculate root mean squared errors and biases
#'
#' used for comparing truth and estimates
#'
#' @param A a matrix (e.g., truth)
#' @param B another matrix (e.g., estimated); A and B must have the same dimensions;
#' B's columns must have already been permuted to best match those of A.
#'
#' @return a list
#' \describe{
#' \item{rmse_total}{root mean squared error for all entries}
#' \item{rmse_marg}{root mean squared error for all entries by column of A}
#' \item{frac_bias}{percent bias averaged over all entries}
#' \item{frac_bias_marg}{ percent bias by column of A}
#' }
#' @examples
#'
#' A <- matrix(c(1,2,3,4,5,6),nrow=2,ncol=3)
#' B <- A
#' rmse_bias_AB(A,B)
#' @export
rmse_bias_AB <- function(A,B){
  if (dim(A)[1]!=dim(B)[1] && dim(A)[2]!=dim(B)[2]){
    stop("[lotR] matrix of different dimensions")}
  rmse_total <- sqrt(mean((A-B)^2))
  rmse_marg <- sqrt(colMeans((A-B)^2))
  frac_bias  <- mean(abs(A-B)/abs(A))
  frac_bias_marg <- colMeans(abs(A-B)/abs(A))
  make_list(rmse_total,rmse_marg,frac_bias,frac_bias_marg)
}

#' calculate the misclassification rates
#'
#' the first is the truth, the second is the estimated ones
#'
#' @param a truth class labels
#' @param b estimated class labels; a and b are of the same length
#'
#' @return the overall and the class-specific misclassification errors
#' @examples
#'
#' compute_misrate(c(1,1,1,2,2,2,3,3,3),c(1,1,1,1,2,2,2,2,2))
#' @export
compute_misrate <- function(a,b){
  if (length(a)!=length(b)){stop("[lotR] the two vectors are of distinct lengths")}
  K <- length(unique(a))
  res <- matrix(NA,nrow=K,ncol=K)
  for (k1 in 1:K){
    for (k2 in 1:K){
      res[k1,k2] <- sum(a==k1 & b==k2)
    }
  }
  n <- length(a)
  c(1-sum(res[cbind(1:K,1:K)])/n,1-(res/rowSums(res))[cbind(1:K,1:K)])
}








