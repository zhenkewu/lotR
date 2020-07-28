# still haven't fixed monotinicity :(/
source("inst/example/lotR_example.R")

V(thetree_igraph)$levels       <- rep(2,length(V(thetree_igraph)))
V(thetree_igraph)$levels[c(1)] <- 1

# V(thetree_igraph)$levels <- rep(1,length(V(thetree_igraph)))
# V(thetree_igraph)$levels[match(names(igraph::V(thetree_igraph)[igraph::degree(thetree_igraph, mode = "out") == 0]),
#                                names(igraph::V(thetree_igraph)))] <- 2
# V(thetree_igraph)$levels[1] <- 3 # separate rootnode.
tau   <- c(0.5,0.3,0.15,0.05)

theta <- rbind(rep(rep(c(0.01, 0.01), each = 1),9),
               rep(rep(c(0.5, 0.5), each = 1),9),
               rep(rep(c(0.1, 0.7), each = 1),9),
               rep(rep(c(0.9, 0.9), each = 1),9))
image(t(theta))

#
# simulated data
#
# Y <- BayesLCA::rlca(2663, itemprob = theta, classprob = tau)
# curr_leaves <- c(dat_mge[!is.na(match_ind),"ct_MLST"])

## five times more: () <---------- is it possible that the order of the outcomes matter in
## terms of estimation. Clean the code and push to github.
Y <- BayesLCA::rlca(2663*5, itemprob = theta, classprob = tau)
curr_leaves <- c(dat_mge[!is.na(match_ind),"ct_MLST"],
                 dat_mge[!is.na(match_ind),"ct_MLST"],
                 dat_mge[!is.na(match_ind),"ct_MLST"],
                 dat_mge[!is.na(match_ind),"ct_MLST"],
                 dat_mge[!is.na(match_ind),"ct_MLST"])


## two times more:
# Y <- BayesLCA::rlca(2663*2, itemprob = theta, classprob = tau)
# curr_leaves <- c(dat_mge[!is.na(match_ind),"ct_MLST"],
#                  dat_mge[!is.na(match_ind),"ct_MLST"]
#                  )


#
# real data:
#
# Y <- rbind(dat_mge[!is.na(match_ind),ind_EL]
# )
# curr_leaves <- c(dat_mge[!is.na(match_ind),"ct_MLST"]
# )

K     <- length(tau)
p     <- length(V(thetree_igraph))
J     <- ncol(Y)

dsgn0  <- design_tree(Y,curr_leaves,thetree_igraph,root_node = "Node1",weighted_edge = FALSE)

nrestarts <- 1
# doParallel::registerDoParallel(cores = nrestarts)
log_dir <- "restart_logs"
dir.create(log_dir)
set.seed(345083)
par(mfrow=c(3,3));plot(0,0)
mod0     <- lcm_tree(Y,curr_leaves,thetree_igraph,
                     rootnode      = "Node1", # <-- may be redundant?
                     weighted_edge = !TRUE,
                     parallel = TRUE,
                     hyper_fixed   = list(K=K,a=c(1,1),b=c(1,1),tau_update_levels = c(99), # two levels, slow; one level faster?
                                          #s_u_zeroset = NULL,s_u_oneset = NULL),
                                          #s_u_zeroset = (1:265)[-c(1,2,3,4,76,120)],s_u_oneset = c(1,2,3,4,76,120)),
                     s_u_zeroset = (1:265)[-c(1)],s_u_oneset = c(1)),
                     # hyperparams_init = list(tau_1=rep(9/4,2),
                     #                        tau_2=rep(9/4,2)),
                     #random_init = TRUE,
                     nrestarts     = nrestarts,
                     print_freq = 1,update_hyper_freq = 50, max_iter = 500,
                     allow_continue = FALSE,log_restarts =!TRUE, log_dir = log_dir)

# unlink(log_dir, recursive = T)
# closeAllConnections()

# summarize posterior results:
plot(mod0)

#image(mod0$prob_est$pi_collapsed)
mod0$prob_est$pi_collapsed
which(mod0$mod$vi_params$prob>0.5)




#colnames(dsgn0$Y) <- c("EL2","EL3", "EL40", "EL35", "EL44", "EL45", "EL46", "EL12", "EL18","EL19", "EL36", "EL37", "EL38", "EL39" ,"EL41", "EL42", "EL43", "EL50")
# blca:
# VB:

fullLCM_vb <- BayesLCA::blca(dsgn0$Y,K,method="vb",verbose=FALSE)
fullLCM_vb
image(t(fullLCM_vb$itemprob),main="BayesLCA package: vb")

# EM:
fullLCM_em <- BayesLCA::blca(dsgn0$Y,K,method="em",verbose=FALSE)
fullLCM_em
image(t(fullLCM_em$itemprob),main="BayesLCA package: em")


barplot(rbind(c(tau,rep(0,K-length(tau))),
              sort(mod0$prob_est$pi_collapsed,decreasing = TRUE),
              fullLCM_vb$classprob,#,
              fullLCM_em$classprob
),
main="class probabilities", beside=TRUE)#,legend.text = c("truth","proposed","vb","em"))


## poLCA:
# form = as.formula("cbind(EL2,EL3, EL40, EL35, EL44, EL45, EL46, EL12, EL18,
#                   EL19, EL36, EL37, EL38, EL39 ,EL41, EL42, EL43, EL50)~1")
# for (j in 1:ncol(dsgn0$Y)){dsgn0$Y[,j] <- as.integer(as.factor(dsgn0$Y[,j]))}
# res2 <- poLCA::poLCA(form, as.data.frame(dsgn0$Y),nclass=K)

proposed <- mod0$prob_est$theta_collapsed


png("inst/example_figure/comparison_with_std.png",width=10,height=9,units = "inches")
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
              mod0$prob_est$pi_collapsed,
              fullLCM_vb$classprob,
              fullLCM_em$classprob),
        main="class probabilities", beside=TRUE,legend.text = c("truth","proposed","vb","em"))
dev.off()


barplot(rbind(c(tau,rep(0,K-length(tau))),
              sort(mod0$prob_est$pi_collapsed,decreasing = TRUE),
              fullLCM_vb$classprob#,
              #fullLCM_em$classprob
),
main="class probabilities", beside=TRUE,legend.text = c("truth","proposed","vb","em"))
