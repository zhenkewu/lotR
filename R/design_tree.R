#######################################################################
# functions to organize data around the tree structure
#######################################################################

#' Organize the data around the rooted binary weighted tree
#'
#' NB: currently this minimally built; need some checking functions
#'
#' @param Y `N` by `J` binary data matrix; rows for subjects; columns
#' for binary measurements/features
#' @param leaf_ids Character string for leaf nodes for each observation
#' @param mytree a tree (an `igraph` object) that contains the node,
#'              edge, edge-length ("weight") information.
#' @param weighted_edge logical: `TRUE` for using the branch lengths then
#' the `mytree` must have this info; if `FALSE`, every edge, including
#' an imaginary edge leading to the root node, is set to have length `1`.
#' @param Z_obs Default is `NULL`; A two-column matrix of (id, class indicator);
#' the number of rows equals the number of observations;
#' an entry of 2nd column is `NA` if the subject
#' in that row has an unknown class indicator.
#' Importantly, the rows will be reordered according to the reordered `Y`.
#'
#' @return A list of data and tree information for model fitting
#' \itemize{
#' \item `Y` A matrix of `N` by `J`; binary measurements with rows ordered by
#'       leaf groups  (`leaf_ids`).
#' \item `A` A matrix of `p` by `p`; each column contains some `1`s, indicating
#' the node in that column is an ancestor of the node represented in the row.
#' Ancestor of a node include that node itself.
#' \item `A_leaves` A matrix of `pL` by `p`; A submatrix of `A` that represents
#' the ancestor matrix but only for leaves
#' \item `leaf_ids` A vector of `N`integers; ordered by the leaves as
#'              specified by `mytree`.
#' \item `leaf_ids_units` A list of length `pL`, each element
#'            is a vector of subject ids belonging to each leaf node
#' \item `leaf_ids_nodes` a list of length `p`, each element
#'            is a vector of integers (between `1` and `pL`; id
#'            is only for leaf nodes) indicating the leaf nodes.
#' \item `ancestors` a list of length `pL`,
#'      each element is the vector of ancestors (between `1` and `p`; id is among
#'      all nodes)
#' \item `edge_lengths` a list of length `pL`,
#'       each element is a numeric vector of edge lengths from the root node
#'       to the leaf. It is computed based on `E(mytree)$weight`. It is `NULL`
#'       if `E(mytree)$weight` is `NULL`
#' \item `h_pau` a numeric vector of length `p`; each value is
#' the edge length from u to its parent (if u is a root node, then the value is 1).
#' This vector by default is all `1`s. If `weighted_edge=TRUE`, `h_pau`
#' is set to `E(mytree)$weight`, the input edge weights.
#' \item `v_units` a vector of length equal to the total number of rows in `Y`;
#' each element is an integer between `1 and `pL`,
#' indicating which leaf does the observation belong to.
#' \item `subject_id_list` a list of length `p`; each element is a vector of
#' subject ids
#' that are in the leaf descendants of node `u` (internal or leaf node)
#' \item `ord` the permutation to order the original rows to produce the final
#' ordering of the rows of `Y`.
#' }
#' @export
#' @import igraph
#'
design_tree <- function(Y,leaf_ids,mytree,weighted_edge=FALSE,Z_obs = NULL){ # by default, not weighted tree.
  # warning("using hard-coded info\n")
  # Y        <- dat_mge[!is.na(match_ind),ind_EL] # these are not ordered yet.
  # leaf_ids <- dat_mge[!is.na(match_ind),"ct_MLST"]
  # mytree   <- thetree_igraph
  # E(mytree)$weight <- thetree$edge.length
  # rootnode="Node1"
  # weighted_edge <- FALSE

  if (!is.character(leaf_ids)) stop("[lotR] `leaf_ids` is not a character object.")
  if (!igraph::is.igraph(mytree)) stop("[lotR] 'mytree' is not a graph object.")
  if (!igraph::is.directed(mytree)) stop("[lotR] 'mytree' is not directed.")

  cat("\n\n [lotR] ", c("unweighted", "weighted")[weighted_edge+1], "tree...\n\n")

  # extract nodes:
  nodes  <- names(igraph::V(mytree))
  leaves <- names(igraph::V(mytree)[igraph::degree(mytree, mode = "out") == 0])
  rootnode <- names(igraph::V(mytree)[igraph::degree(mytree, mode = "in") == 0])
  if(!setequal(unique(leaf_ids), leaves))
  {stop("[lotR] Not all `leaf_ids` are leaves of tree")}

  nodes0 <- nodes
  # Re-order nodes to have internal nodes first, then leaves: <------- ordered nodes.
  nodes <- c(nodes[!(nodes %in% leaves)], leaves)




  # Get levels for specifying groups of hyperparameters
  if (is.null(igraph::V(mytree)$levels)) {
    # The default is to have two levels: one for leaf nodes, one for internal nodes
    levels <- rep(1, length(nodes))
    levels[nodes %in% leaves] <- 2
  } else {
    # Otherwise, use the levels supplied; can be a single level too:
    levels <- igraph::V(mytree)[nodes]$levels    # <--  ordered.
    levels <- levels[match(names(igraph::V(mytree)), nodes)] # <--- ordered. same as above???
  }
  if (sum(table(levels) < 5) > 0) {
    warning("[lotR] Some levels contain fewer than five nodes: this may lead to problems
            with estimating variance parameters. Recommend increasing the number
            of nodes per level.")
  }

  # sizes:
  p  <- length(nodes)
  pL <- length(leaves)
  n  <- nrow(Y)

  # ancestor matrix:
  A <- igraph::as_adjacency_matrix(mytree, sparse = TRUE) # directed graph; upper diag mat
  A <- A[nodes, nodes] # re-order rows/columns to mirror nodes
  A <- Matrix::expm(Matrix::t(A)) # get ancestors: lower diag mat (1: column j is anc for the row)
  A[A > 0 ] <- 1
  A <- Matrix::Matrix(A, sparse = TRUE)


  # 16 x 16 sparse Matrix of class "dtCMatrix"
  # [[ suppressing 16 column names ‘7.4’, ‘7.4.1’, ‘7.4.2’ ... ]]
  #
  # 7.4     1 . . . . . . . . . . . . . . .
  # 7.4.1   1 1 . . . . . . . . . . . . . .
  # 7.4.2   1 . 1 . . . . . . . . . . . . .
  # 7.4.3   1 . . 1 . . . . . . . . . . . .
  # 7.4.4   1 . . . 1 . . . . . . . . . . .
  # 7.4.1.1 1 1 . . . 1 . . . . . . . . . .
  # 7.4.1.2 1 1 . . . . 1 . . . . . . . . .
  # 7.4.1.3 1 1 . . . . . 1 . . . . . . . .
  # 7.4.2.0 1 . 1 . . . . . 1 . . . . . . .
  # 7.4.2.1 1 . 1 . . . . . . 1 . . . . . .
  # 7.4.2.2 1 . 1 . . . . . . . 1 . . . . .
  # 7.4.3.0 1 . . 1 . . . . . . . 1 . . . .
  # 7.4.3.1 1 . . 1 . . . . . . . . 1 . . .
  # 7.4.3.2 1 . . 1 . . . . . . . . . 1 . .
  # 7.4.4.1 1 . . . 1 . . . . . . . . . 1 .
  # 7.4.4.2 1 . . . 1 . . . . . . . . . . 1

  # Sort by leaf_ids, where order is specified by ordering in 'mytree':
  ord <- order(ordered(leaf_ids, levels = leaves)) # <--- leaves.
  Y   <- Y[ord,,drop=FALSE]
  leaf_ids <- leaf_ids[ord] # length = n.
  if (!is.null(Z_obs)){
    Z_obs <- Z_obs[ord,,drop=FALSE]
  }
  # get lists of ancestors for each leaf_ids:
  d <- igraph::diameter(mytree,weights=NA)
  # need to set weight=NA to prevent the use of edge lengths in determining the diameter.
  ancestors <- igraph::ego(mytree, order = d + 1, nodes = leaves, mode = "in")
  ancestors <- sapply(ancestors, names, simplify = F)
  ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a),
                      nodes = nodes, simplify = F)
  names(ancestors) <- leaves

  # get lists of which units correspond to each leaf_ids
  leaf_ids_units <- sapply(leaves, function(v) which(leaf_ids == v), simplify = FALSE)
  names(leaf_ids_units) <- leaves

  # get lists of leaf_ids are descendants of each node
  descendants <- igraph::ego(mytree, order = d + 1, nodes = nodes, mode = "out")
  descendants <- sapply(descendants, names)

  # leaf.descendants:
  leaf_ids_nodes <- sapply(descendants, function(d, leaves) which(leaves %in% d),
                           leaves = leaves,simplify = FALSE)

  # Change leaf_ids to integers
  leaf_ids <- sapply(leaf_ids, function(o) which(leaves == o))

  # The branch lengths (NB: need to check this is correct.)
  edge_lengths        <- vector("list",length=pL)
  names(edge_lengths) <- names(ancestors)
  for (j in 1:pL){
    res <- shortest_paths(mytree, # here we are using the original mytree, not "nodes" which is reordered.
                          which(vertex_attr(mytree)$name==rootnode), # root.
                          which(vertex_attr(mytree)$name==leaves[j]), # leaf.
                          output="both")
    edge_lengths[[j]] <- E(mytree)$weight[sapply(res$epath[[1]],
                                                 function(e) which(E(mytree) == e))]
  }

  parent_nm <- c(rootnode,
                 names(unlist(lapply(igraph::ego(mytree, order = 1,
                                                 nodes = nodes, mode = "in")[-1],"[",-1))))
  h_pau       <- rep(1,p) # equally weighted with default lengths 1.
  if (weighted_edge){
    for (u in 2:p){
      res <- shortest_paths(mytree, # here we are using the original mytree, not "nodes" which is reordered.
                            which(vertex_attr(mytree)$name==parent_nm[u]), #parent
                            which(vertex_attr(mytree)$name==nodes[u]), # leaf.
                            output="both")

      h_pau[u] <- E(mytree)$weight[sapply(res$epath[[1]],
                                          function(e) which(E(mytree) == e))]
    }
    h_pau <- pmax(h_pau,1e-8) # just to make the weights non-zero. <-----------------------------------------IS THIS A MUST?
  }
  #h_pau[1] <- mean(h_pau[-1])

  A_leaves <- A[nodes%in%leaves,]
  v_units <- NULL
  for (v in 1:pL){
    leaf_list_tmp <- leaf_ids_units[v]
    units         <- unlist(leaf_list_tmp)
    v_units       <- c(v_units,unlist(mapply(rep,v,unlist(lapply(leaf_list_tmp,length)))))
  }

  # leaf_ids_units: pL length, each is the ids of individuals.
  # leaf_ids_nodes: p length, each is the indicator (between 1 and pL) the leaves.
  subject_id_list <- list() # for each internal or leaf node
  for (u in 1:p){
    leaf_desc            <- leaf_ids_nodes[[u]]
    subject_id_list[[u]] <- unlist(leaf_ids_units[leaf_desc])
    # if (u==120){
    #   for (v in leaf_desc){
    #     print(leaves[v])
    #     print(c(unlist(leaf_ids_units[v])))
    #   }
    # }
  }

  # print(subject_id_list)
  mytree <- permute(mytree,match(nodes,nodes0))
  make_list(Y,A,A_leaves,rootnode,
            leaf_ids,leaf_ids_units,leaf_ids_nodes,
            ancestors,edge_lengths,h_pau, # h_pau is used in fitting, so it can be that h_pau are 1s but the edge lengths are weighted.
            levels,v_units,subject_id_list,Z_obs,ord,mytree)
  #NB: by now, the node ids have already been reordered; it is only that the tree
  # itself has not ordered the ids in the tree.
}






