# This code is not designed as an example for a particular function,
# but it simulates the example data set: lotR_example_data_tree

library(igraph)
 set.seed(2020)
 n = 3000
 tau   <- c(0.6,0.3,0.1)
 itemprob <- rbind(rep(rep(c(0.9, 0.9), each = 1),9),
                   rep(rep(c(0.5, 0.5), each = 1),9),
                   rep(rep(c(0.1, 0.1), each = 1),9))

 data("lotR_example_edges")
 mytree <- igraph::graph_from_edgelist(lotR_example_edges, directed = TRUE)
 # Plot tree

 nodes  <- names(igraph::V(mytree))
 leaves <- names(igraph::V(mytree)[degree(mytree, mode = "out") == 0])
 pL = length(leaves)
 p  = length(igraph::V(mytree))

 ######
 K = nrow(itemprob)
 # specify the nodes that have non-trivial alpha_u, this was called
 # xi_u, because xi_u = s_u*alpha_u and s_u = 1 if we set it in the simulation.
 alpha_mat = rbind(logit(prob2stick(tau)[-K]),
                   c(-0.5,-1),
                   c(0.5,1),
                   matrix(0,nrow=p-3,ncol=K-1)
 )

 # get lists of ancestors for each leaf_ids:
 d <- igraph::diameter(mytree,weights=NA)
 # need to set weight=NA to prevent the use of edge lengths in determining the diameter.
 ancestors <- igraph::ego(mytree, order = d + 1, nodes = leaves, mode = "in")
 ancestors <- sapply(ancestors, names, simplify = FALSE)
 ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a),
                     nodes = nodes, simplify = FALSE)
 names(ancestors) <- leaves

 # calculate the class probabilities for all leaf nodes; each leaf node
 # should have a K-dim vector that sums to one; Some nodes may share
 # the same set of K-dim probability vector, others may differ. There are
 # one or more groups of leaf nodes with distinct K-dim probability vectors.
 # Note the branch lengths may also be used here.
 pi_mat <- matrix(NA,nrow=pL,ncol=K)
 for (v in seq_along(leaves)){
   pi_mat[v,-K] <- colSums(alpha_mat[ancestors[[v]],,drop=FALSE])
   pi_mat[v,] <- tsb(c(expit(pi_mat[v,-K]),1))
 }

 # s = c(1, 1,1,0,0, rep(0,pL)) # effective nodes
 h_pau = rep(1,p)

 lotR_example_data_tree <- simulate_lcm_tree(n,itemprob,mytree,pi_mat,h_pau)
 table(lotR_example_data_tree$truth$Z)

 #save the simulated data to the R package for illustration:
 save(lotR_example_data_tree, file = "data/lotR_example_data_tree.rda", compress = "xz")
