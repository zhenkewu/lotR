### Example workflow of `lotR`
### Zhenke Wu | zhenkewu@gmail.com
### August 04, 2020
### This example has true parameters following the tree structured priors.
rm(list=ls())
library(lotR)

data("lotR_example_edges")
library(igraph)
tr <- graph_from_edgelist(lotR_example_edges, directed = TRUE) # Plot tree
# Here in this example, we use equal weighted edges.

# If needed, first install ggtree
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)
ggtree(tr, ladderize = FALSE, layout = "slanted") +
  geom_tiplab(geom = "label") + geom_nodelab(geom = "label") +
  theme(plot.margin = unit(c(0, 1.5, 0, 0.2), "cm")) +
  coord_cartesian(clip = "off") + scale_y_reverse()

leaves <- names(igraph::V(tr)[degree(tr, mode = "out") == 0])

# set levels: leaves 2, non-leaves 1
igraph::V(tr)$levels <- rep(1,length(igraph::V(tr)))
igraph::V(tr)$levels[match(names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0]),
                   names(igraph::V(tr)))] <- 2

###############################################################################
## illustration of lotR on simulated data:
###############################################################################
data("lotR_example_data_tree") # assumes a single pi for all observations.
Y     <- lotR_example_data_tree$Y
Z_obs <- lotR_example_data_tree$Z_obs
curr_leaves <- lotR_example_data_tree$curr_leaves
theta <- lotR_example_data_tree$truth$itemprob
pi_mat <- lotR_example_data_tree$truth$pi_mat
# in the simulation, the 1st class has the highest level of probabilities,
# so we will reorder obtained estimates to match the meaning of classes.
Z     <- lotR_example_data_tree$truth$Z

## check if all leave names are present in the data.
setequal(leaves, unique(curr_leaves))
#> [1] TRUE
knitr::kable(table(curr_leaves))

K     <- nrow(theta)
p     <- length(V(tr))
J     <- ncol(Y)

dsgn0  <- design_tree(Y,curr_leaves,tr,weighted_edge = FALSE,Z_obs)

nrestarts <- 1
# doParallel::registerDoParallel(cores = nrestarts)
# log_dir <- tempdir()
# dir.create(log_dir)
set.seed(345083)
par(mfrow=c(3,3));plot(0,0)
mod0     <- lcm_tree(Y,curr_leaves,tr,
                     weighted_edge = !TRUE,
                     Z_obs         = Z_obs,
                     parallel      = TRUE,
                     hyper_fixed   = list(K=K,a=c(1,1),b=c(1,1),
                                          tau_update_levels = c(1,2),
                                          # should this be set as default?
                                          s_u_zeroset = NULL,
                                          s_u_oneset = c(1)),
                     # hyperparams_init = list(tau_1=matrix(9/4,nrow=2,ncol=K-1),
                     #                         tau_2=array(9/4,c(2,J,K))),
                     vi_params_init = list(prob=rep(0.95,p)),
                     random_init = !TRUE,
                     nrestarts     = nrestarts,
                     quiet         = !TRUE,
                     plot_fig      = TRUE,
                     shared_tau    = FALSE,
                     get_lcm_by_group = TRUE, #<--- check what are the prereq.
                     print_freq    = 10,update_hyper_freq = 50, max_iter = 5000,
                     tol           = 1e-6,
                     tol_hyper     = 1e-3,
                     allow_continue = FALSE)#,
#log_restarts =!TRUE,
#log_dir = log_dir)

###############################################################################
# print out summaries of group specific estimates:
###############################################################################
plot(mod0,layout = "slanted", horizontal = FALSE)

print(mod0)

which(mod0$mod$vi_params$prob>0.5)
unique(lotR_example_data_tree$truth$pi_mat)

