### Example workflow of `lotR`
### Zhenke Wu | zhenkewu@gmail.com
### August 02, 2020
rm(list=ls())
library(lotR)

# library(moretrees)
# data("moretreesExampleEdgelist")
# lotR_example_edges = moretreesExampleEdgelist
# save(lotR_example_edges, file = "data/lotR_example_edges.rda", compress = "xz")

data("lotR_example_edges")
library(igraph)
tr <- graph_from_edgelist(lotR_example_edges, directed = TRUE) # Plot tree
# here in this example, we use equal weighted edges.

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

# data("moretreesExampleData")
# Xcase <- moretreesExampleData$Xcase
# Xcontrol <- moretreesExampleData$Xcontrol
# Wcase <- moretreesExampleData$Wcase
# Wcontrol <- moretreesExampleData$Wcontrol
# outcomes <- moretreesExampleData$outcomes
leaves <- names(V(tr)[degree(tr, mode = "out") == 0])

# set levels: leaves 2, non-leaves 1
V(tr)$levels <- rep(1,length(V(tr)))
V(tr)$levels[match(names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0]),
                               names(igraph::V(tr)))] <- 2

###############################################################################
## illustration of lotR on simulated data:
###############################################################################

# # <------ DETAILS OF SIMULATED DATA:
# # simulation truth:
# tau   <- c(0.6,0.3,0.1)
# theta <- rbind(rep(rep(c(0.9, 0.9), each = 1),9),
#                rep(rep(c(0.5, 0.5), each = 1),9),
#                rep(rep(c(0.1, 0.1), each = 1),9))
#
# # tau   <- c(0.5,0.3,0.15,0.05)
# #
# # theta <- rbind(rep(rep(c(0.9, 0.9), each = 1),9),
# #                rep(rep(c(0.01, 0.01), each = 1),9),
# #                rep(rep(c(0.5, 0.5), each = 1),9),
# #                rep(rep(c(0.1, 0.7), each = 1),9)
# # )
#
# image(t(theta))
#
# fold     <- 1
# n        <- 1000
# prop_obs <- 0.0
# simu <- lotR_blca(n*fold, itemprob = theta, classprob = tau)
# Y <- simu$x
# Z <- simu$z
# Z_obs  <- as.matrix(cbind(1:nrow(Y),Z))
# if (prop_obs*nrow(Y)<1){
#   Z_obs <- NULL
# }else{
#   Z_obs[-sort(sample(1:nrow(Y),floor(prop_obs*nrow(Y)),replace=FALSE)),2] <- NA
# }
#
# curr_leaves <- sample(leaves,nrow(Y),replace = TRUE)
#
# lotR_example_data <- list(
#   Y = Y,
#   Z_obs = Z_obs,
#   curr_leaves = curr_leaves,
#   truth = make_list(tau,theta,prop_obs,Z)
# )
# # save the simulated data to the R package for illustration:
# save(lotR_example_data, file = "data/lotR_example_data.rda", compress = "xz")
# <------ DETAILS OF SIMULATED DATA.

data("lotR_example_data") # assumes a single pi for all observations.

Y     <- lotR_example_data$Y
Z_obs <- lotR_example_data$Z_obs
curr_leaves <- lotR_example_data$curr_leaves
tau <- lotR_example_data$truth$tau
theta <- lotR_example_data$truth$theta
# in the simulation, the 1st class has the highest level of probabilities,
# so we will reorder obtained estimates to match the meaning of classes.
Z     <- lotR_example_data$truth$Z

## check if all leave names are present in the data.
setequal(leaves, unique(curr_leaves))
#> [1] TRUE
knitr::kable(table(curr_leaves))

K     <- length(tau)
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
                                          #s_u_oneset = NULL),
                                          #s_u_zeroset = (1:p)[-c(1)],
                                          s_u_oneset = c(1)),
                     # hyperparams_init = list(tau_1=rep(9/4,2),
                     #                         tau_2=rep(9/4,2)),
                     # hyperparams_init = list(tau_1=matrix(9/4,nrow=2,ncol=K-1),
                     #                          tau_2=array(9/4,c(2,J,K))),
                     vi_params_init = list(prob=rep(0.95,p)),
                     random_init = !TRUE,
                     nrestarts     = nrestarts,
                     quiet         = !TRUE,
                     plot_fig      = !TRUE,
                     shared_tau    = FALSE,
                     get_lcm_by_group = TRUE, #<--- check what are the prereq.
                     print_freq    = 10,update_hyper_freq = 50, max_iter = 5000,
                     tol           = 1e-8,
                     tol_hyper     = 1e-4,
                     allow_continue = FALSE)#,
                     #log_restarts =!TRUE,
                     #log_dir = log_dir)

###############################################################################
# print out summaries of group specific estimates:
###############################################################################
plot(mod0,layout = "slanted", horizontal = FALSE)

mod0$prob_est$pi_collapsed
which(mod0$mod$vi_params$prob>0.5)


###############################################################################
# Compare with other existing methods (Variational Bayes and EM)
# We can do this comparison because the simulation truth is a single LCM.
# The lotR package can start from an unknown grouping, and shrink towards a
# single group.
###############################################################################

# after we check the grouping is identical to truth, we
# compare the proposed method against other methods (VB, EM):


# permute the estimated classes to best match the truth in terms of the response
# probability profiles.
itemprob_truth   <- t(theta)
itemprob_est     <- mod0$prob_est$theta_collapsed
# we set prob_gamma to be all zero except for the root node;
# so theta_collapsed is just the shared response probability profiles.
ord_rearrange      <- opt_colpermB(itemprob_truth,itemprob_est)
# permute estimated classes to best match truth.

proposed <- itemprob_est[,ord_rearrange]

# blca:
# VB:
fullLCM_vb <- BayesLCA::blca(dsgn0$Y,K,method="vb",verbose=FALSE)
# EM:
fullLCM_em <- BayesLCA::blca(dsgn0$Y,K,method="em",verbose=FALSE)

#
# compare via plots:
#
# png("inst/example_figure/comparison_with_std.png",width=10,height=9,units = "in",res=300)
par(mfcol=c(3,3))
image(t(theta),main="truth")

# differences
image((proposed-t(theta)),main="diff")
hist((proposed-t(theta)),xlim=c(-1,1))

image(proposed,main="proposed")
image(t(fullLCM_vb$itemprob),main="BayesLCA package: vb")
image(t(fullLCM_em$itemprob),main="BayesLCA package: em")

# pi estimate:
barplot(rbind(c(tau,rep(0,K-length(tau))),
              # mod0$prob_est$pi_collapsed,
              sort(mod0$prob_est$pi_collapsed,decreasing = TRUE),
              fullLCM_vb$classprob,
              fullLCM_em$classprob),
        main="class probabilities", beside=TRUE)#,
        #legend.text = c("truth","proposed","vb","em"))


## individual-specific variational probabilities of belonging to each class:
image(mod0$mod$vi_param$rmat[order(dsgn0$ord),ord_rearrange],
      main="inferred individual class prob",
      xlab="subjects",xaxt="n",yaxt="n")
image(BayesLCA:::unMAP(Z),main="true individual class memberships",
      xlab="subjects",xaxt="n",yaxt="n")
# dev.off()

plot(mod0$mod$ELBO_track,xlab="iteration",ylab="ELBO",
     main="Evidence Lower Bound",type="l")

# ## show the ELBO trajectory:
# library(plotly)
# tmp_df <- data.frame(iteration = 1:length(mod0$mod$ELBO_track),
#                      ELBO=mod0$mod$ELBO_track)
# plot_ly(tmp_df,x=~iteration,y = ~ELBO,type = 'scatter', mode = 'lines')




