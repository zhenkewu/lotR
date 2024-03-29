#######################################################################
# functions for initialization:
# 1. variational parameters (intermediate moments):
#     prob (p_u), (a_t,b_t) - for update rho, r_i
#     (tau_1_t, tau_2_t) - the variance for the gamma or alpha parameters
#     mu_gamma, mu_alpha, sigma_gamma, Sigma_alpha
# 2. hyperparameters:
#    - psi, phi, g_psi, g_phi - these are actually local variational parameters
#      but they are updated in the same function as tau_1 and tau_2, but more frequently.
#    - (tau_1, tau_2) - pre-specified hyperparameters - not fixed.
# 3. hyperfixed: constants, e.g., number of classes K; for
#        rho_l - beta(a_l,b_l), the (a_l,b_l) are fixed.
#
# NB: 1) need good initialization and random initialization
#     2) need to use fast code to fit LCM using particular groupings (BayesLCA)
#     3) theoretically, only a subset of parameter/hyperparameters
#        needs to be initialized, because some updates in the first iteration
#       can produce a valid value. But for completeness, it makes
#       sense to have a full initialization overall parameters. This
#       helps if the updating order is changed.
#     4) here should have a few functions
#        a) initialization function
#        b) update variational parameters
#        c) update the hyperparameters every d steps
#        d) a function that wraps a and b with initial values
#        e) a main wrapper function that reads in data, organize
#           data around the tree, check inputs for compatibility,
#           initialize, update the model parameters and hyperparameters,
#           summarize the posterior, compute separate LCMs based on the discovered
#           grouping of the observations.
#######################################################################


#' Initialize the variational inferential algorithm for latent
#' class model with tree structured shrinkage
#'
#' @param Y matrix of binary observations for LCM - rows ordered by leaves in the tree
#' @param A a p by p binary matrix: each row is the ancestors including the node itself;
#'          ordered by leaves in the tree
#' @param Z_obs a two-column matrix of integers, first column is subject ids,
#' second column is a mix of NA and integers, NA means unknown class indicators,
#' an integer indicates the known class.
#' @param leaf_ids_units The subject ids in each leaf node (a list)
#' @param leaf_ids_nodes the leaf descendants for each internal or leaf nodes (a list)
#' @param ancestors a numeric vector of ancestor nodes for each leaf node (a list of length equal to
#' the number of leaves)
#' @param h_pau a numeric vector (length = # nodes); the edge length between a node
#' and its parent node. The root node has no parent, because we suggest a separate
#' prior variance for the root node's gamma and alpha, we set the "edge-length" toward root
#' node as 1.
#' @param levels a numeric vector of integers from 1 to L, indicating for each node
#' (leaf or internal node) which set of hyperparameters to use. For example,
#' if we want the root node to have a separate `tau_1` and `tau_2`, can specify it to
#' its own level. Another example would be to have distinct sets of hyperparameters
#' for leaf and non-leaf nodes. The levels are pre-specified, not estimated.
#' @param vi_params the list of variational parameters. `mu_gamma`,
#' `mu_alpha`, `prob` (for s_u), `a_t`, `b_t`,`sigma_gamma`, `Sigma_alpha`
#' `prob`,`prob_gamma=c(1,rep(0,p-1))` (this is not updated).
#' @param hyperparams the list of hyperparameters, `tau_1` and `tau_2` -
#' these are initial specifications of the hyperparameters - they are updated by
#'  `tau_1_t`, `tau_2_t`; `psi`, `g_psi`, `phi`, `g_phi` (these
#' are not hyperparameters, but are updated in the same function as `tau_1_t` and `tau_2_t`),
#' but more frequently.
#' @param hyper_fixed a list of fixed hyperparameters, such as those
#' in the Beta priors for `rho`, e.g., `list(a=c(1,1,99),b=c(1,1,1))`,
#' `K`, number of classes.
#' for three levels, where the first and second levels are uniform Beta, the third level
#' has `Beta(99,1)` prior - this nearly ensures setting `s_u=1`.
#' @param random_init logical; `TRUE` for adding additional variability to the
#' initial values. This is needed if the algorithm needs multiple random starts
#' to choose the best converged values.
#' @param random_init_vals NB: fill out specific elements
#' @param subject_id_list,v_units,shared_tau see [update_vi_params()]
#'
#' @importFrom BayesLCA blca unMAP
#' @importFrom stats runif rnorm
#' @family internal function
#' @return a list `vi_params,hyperparams` containing the initial values.
#'
#' @export
initialize_tree_lcm <- function(Y,A,Z_obs,
                                leaf_ids_units,
                                leaf_ids_nodes,
                                ancestors,
                                h_pau,
                                levels,
                                vi_params,
                                hyperparams,
                                hyper_fixed,
                                random_init,
                                random_init_vals,
                                subject_id_list,v_units,shared_tau
){
  X <- 2*Y-1
  n <- nrow(Y)
  J <- ncol(Y)
  p <- length(unique(unlist(ancestors)))
  pL <- length(ancestors)
  Fg <- max(levels)
  cat("\n number of levels of hyperparameters: L = ", Fg," \n")
  if (!is.null(hyper_fixed$tau_update_levels)){
    cat("\n update tau_1/2 for levels: ", hyper_fixed$tau_update_levels, " \n")
  }
  if (!shared_tau){
    cat("\n For each node level, tau_1 of K-1 dimension; tau_2 is J by K \n")
  }else{
    cat("\n For each node level, tau_1 and tau_2 are both numeric \n")
  }
  K  <- hyper_fixed$K

  ## initialize:
  if (is.null(vi_params[["mu_gamma"]]) | is.null(vi_params[["mu_alpha"]])){
    # beta_map <- array(rnorm(p*J*K,0,1),c(p,J,K))         # transformed parameters for response probs.
    # eta_map  <- matrix(rnorm(p*(K-1)),nrow=p,ncol=K-1) # transformed parameters for class probs.
    # # fullLCM <- blca(Y,K,method="vb",verbose=FALSE)
    # # eta_map[1,]   <-  as.numeric(alr(colMeans(fullLCM$Z)))
    # # beta_map[1,,] <-  t(logit(fullLCM$itemprob))

    # regular
    beta_map <- array(NA,c(p,J,K))         # transformed parameters for response probs.
    eta_map  <- matrix(NA,nrow=p,ncol=K-1) # transformed parameters for class probs.
    start_class_member_prob <- matrix(runif(n*K),nrow=n,ncol=K)
    start_class_member_prob <- start_class_member_prob/rowSums(start_class_member_prob)
    if (!is.null(Z_obs)){
      for (i in which(!is.na(Z_obs[,2]))){
        start_class_member_prob[i,Z_obs[i,2]] <- 1
        start_class_member_prob[i,-Z_obs[i,2]] <- 0
      }
    }
    for (v in 1:p){
      u   <- leaf_ids_nodes[[v]] # this might be few for some leaves
      # Amazingly this can be fitted by VB method in BayesLCA.
      units <- unlist(leaf_ids_units[u])
      #invisible(capture.output(
      sub_id <- match(Z_obs[!is.na(Z_obs[,2]),1],units) # which in units are in Z_obs.
      if (!is.null(Z_obs) && sum(!is.na(sub_id))>0){
        goodstart <- start_class_member_prob[units,,drop=FALSE]
        mod <- blca(Y[units,,drop=FALSE],K,method="vb",verbose=FALSE,start.vals = goodstart)#)) # <---- This has random intialization, causing different initialization.
      } else{
        mod <- blca(Y[units,,drop=FALSE],K,method="vb",verbose=FALSE,start.vals="single")
      }
      # but this should not matter if hyperparams and vi_params are provided.
      tmp          <- sapply(mod$classprob,function(s) min(max(s,0.01),0.99))
      tau          <- (tmp/sum(tmp))
      eta_map[v,]  <- logit(prob2stick(tau)[-length(tau)])
      beta_map[v,,] <- logit(pmin(pmax(t(mod$itemprob),0.01),0.99))
      # cat("node ",v," initialized as:\n")
      # barplot(tau,main=v)
      # image(t(mod$itemprob),main=v)
      rm("tau")
    }
    # replace infinite values:
    # beta_marg <- aperm(array(logit(pmin(pmax(colMeans(Y),0.001),0.9999)),
    #                          c(J,K,p)),c(3,1,2))
    # eta_all <- t(replicate(p,c(eta_map[1,])))
    #
    # beta_map[is.infinite(beta_map)] <- beta_marg[is.infinite(beta_map)]
    # #beta_map[is.infinite(beta_map) & beta_map>0] <-   5
    # eta_map[is.infinite(eta_map)] <- eta_all[is.infinite(eta_map)]
    # #eta_map[is.infinite(eta_map) & eta_map>0] <-   5
    # # replace NaN values:
    # beta_map[is.na(beta_map)] <- beta_marg[is.na(beta_map)] #runif(length(is.na(beta_map))) #0?
    # eta_map[is.na(eta_map)]   <-  eta_all[is.na(eta_map)]#runif(length(is.na(eta_map)))
    # #print(eta_map)

    # replace infinite values:
    beta_map[is.infinite(beta_map) & beta_map<0] <- - 5
    beta_map[is.infinite(beta_map) & beta_map>0] <-   5
    eta_map[is.infinite(eta_map) & eta_map<0] <- - 5
    eta_map[is.infinite(eta_map) & eta_map>0] <-   5
    # replace NaN values:
    beta_map[is.na(beta_map)] <- 0
    beta_map[is.na(beta_map)] <- 0
    eta_map[is.na(eta_map)]   <- 0
    eta_map[is.na(eta_map)]   <- 0

    # transform from estimates to increments (eta to xi; beta to zeta):
    A_inv    <- solve(A)
    mu_gamma <- apply(beta_map,c(2,3),function(mat) as.matrix(A_inv%*%mat))
    mu_alpha <- apply(eta_map,2,function(v) as.matrix(A_inv%*%matrix(v,ncol=1)))
  }
  if (is.null(vi_params[["mu_gamma"]])){ # this means we have run the above block of code.
    # convert to list:
    vi_params$mu_gamma <- split_along_dim(mu_gamma,1)
  } else{
    check <- is.list(vi_params$mu_gamma) &&
      sum(sapply(vi_params$mu_gamma,function(x) sum(dim(x)==c(J,K))==2)) ==p
    if (!check) stop("[] incompatible initial value for 'mu_gamma'.")
  }
  if (random_init){
    vi_params$mu_gamma <- lapply(vi_params$mu_gamma,
                                 function(mu) matrix(mu+rnorm(J*K,sd = c(abs(mu))*random_init_vals$mu_gamma_sd_frac),nrow=J,ncol=K))
  }

  if (is.null(vi_params[["mu_alpha"]])){
    vi_params$mu_alpha <- lapply(1:p, function(v, mu) c(mu[v, ]),mu=mu_alpha)
  } else{
    check <- is.list(vi_params$mu_alpha) &&
      sum(sapply(vi_params$mu_alpha,function(x) length(x)==K-1)) == p
    if (!check) stop("[] incompatible initial value for 'mu_alpha'.")
  }
  if (random_init){
    vi_params$mu_alpha <- lapply(vi_params$mu_alpha,
                                 function(mu) mu+rnorm(length(mu),sd = abs(mu)*random_init_vals$mu_alpha_sd_frac))
  }

  # Initialize hyper-parameters:
  if (shared_tau){  # <--------------------------------------------------------- IF SHARED TAU_1's in a node.
    if (is.null(hyperparams[["tau_1"]])){ # tau for alphas.
      hyperparams$tau_1 <- sapply(1:Fg,
                                  function(l) mean((unlist(vi_params$mu_alpha[levels==l])^2/
                                                      rep(unlist(h_pau[levels==l]),each=K-1))))
    }else{
      check <- is.numeric(hyperparams$tau_1) &&
        length(hyperparams$tau_1) == Fg
      if (!check) stop("[] Incompatible intial values for 'tau_1' - mlogit class probabilities (alpha).")
    }
    if (random_init){
      hyperparams$tau_1 <- sapply(hyperparams$tau_1,
                                  function(tau) runif(1,min = tau*random_init_vals$tau_lims[1],
                                                      max = tau*random_init_vals$tau_lims[2]))
    }

    if (is.null(hyperparams[["tau_2"]])){ # tau for gammas
      mu_gamma_sq_over_h <- mapply(FUN=function(mat,h){mat^2/h}, mat=vi_params$mu_gamma, h = h_pau,SIMPLIFY=F)
      hyperparams$tau_2 <- sapply(1:Fg,function(l) mean(unlist(mu_gamma_sq_over_h[levels==l]))) # currently very large upon initialization
    }else{
      check <- is.numeric(hyperparams$tau_2) &&
        length(hyperparams$tau_2) == Fg
      if (!check) stop("[] Incompatible intial values for 'tau_2' - logit response probabilities (gamma).")
    }
    if (random_init){
      hyperparams$tau_2 <- sapply(hyperparams$tau_2,
                                  function(tau) runif(1,min = tau*random_init_vals$tau_lims[1],
                                                      max = tau*random_init_vals$tau_lims[2]))
    }
  } else{        # <--------------------------------------------------------- IF DISTINCT TAU_1's in a node.
    # Initialize hyper-parameters:
    if (is.null(hyperparams[["tau_1"]])){ # tau for alphas.
      hyperparams$tau_1 <- matrix(NA,nrow=Fg,ncol=K-1)

      meanfun<- rowMeans
      if (K==2){meanfun<-mean}
      for (l in 1:Fg){
        hyperparams$tau_1[l,]  <- meanfun(mapply(function(v1,v2){v1^2/v2},
                                                                v1 = vi_params$mu_alpha[levels==l],
                                                                v2 = h_pau[levels==l]))
        } # <--- this has problem if K=2; so need to fix.

    }else{
      check <- sum(abs(dim(hyperparams[["tau_1"]]) - c(Fg,K-1))<0.1)==2
      if (!check) stop("[] Incompatible intial values for 'tau_1' - mlogit class probabilities (alpha).")
    }
    if (random_init){
      hyperparams$tau_1 <- hyperparams$tau_1*matrix(runif(Fg*(K-1),
                                                          min = random_init_vals$tau_lims[1],
                                                          max = random_init_vals$tau_lims[2]),nrow=Fg,ncol=K-1)
    }

    if (is.null(hyperparams[["tau_2"]])){ # tau for gammas
      mu_gamma_sq_over_h <- mapply(FUN=function(mat,h){mat^2/h},
                                   mat=vi_params$mu_gamma, h = h_pau,SIMPLIFY=F)
      hyperparams$tau_2 <- array(0,c(Fg,J,K))
      for (l in 1:Fg){
          hyperparams$tau_2[l,,] <- Reduce("+",mu_gamma_sq_over_h[levels==l])/sum(levels==l)
      }
    }else{
      check <- is.numeric(hyperparams$tau_2) &&
        sum(dim(hyperparams$tau_2) == c(Fg,J,K))==3
      if (!check) stop("[] Incompatible intial values for 'tau_2' - logit response probabilities (gamma).")
    }
    if (random_init){
      hyperparams$tau_2 <- hyperparams$tau_2*array(runif(Fg*J*K,min = random_init_vals$tau_lims[1],
                                                         max = random_init_vals$tau_lims[2]),c(Fg,J,K))
    }
  }


  ## Initialize variational prob for s_u:
  if (is.null(vi_params[["prob"]])){
    vi_params$prob <- rep(0.5,p)
  }else{
    check <- is.numeric(vi_params$prob) &&
      length(vi_params$prob) == p &&
      sum(vi_params$prob >=0) == p &&
      sum(vi_params$prob <=1) == p
    if (!check) stop("[] Incompatible initial values for 'p_u' (for variational probability of s_u=1).")
  }
  if (random_init){
    prob <- vi_params$prob
    prob[prob>0.99] <- 0.99
    prob[prob<0.01] <- 0.01
    u <- log(prob/(1-prob))
    u <- u + rnorm(p)*random_init_vals$u_sd_frac*abs(u)
    vi_params$prob <- expit(u)
  }

  ## forcing identical gammas, i.e., response probability profiles across the nodes.
  if (is.null(vi_params[["prob_gamma"]])){ # this is used in update_vi_param and update_hyperparam
    vi_params$prob_gamma <- c(1,rep(0,p-1))
  }

  # the following code segment for a_t, b_t is identical to moretrees:
  if (is.null(vi_params[["a_t"]])) {
    vi_params$a_t <- numeric(Fg)
    for (f in 1:Fg) {
      # initialise these parameters using VI updates
      vi_params$a_t[f] <- hyper_fixed$a[f] + sum(vi_params$prob[levels == f])
    }
  } else {
    check <- is.numeric(vi_params$a_t) &&
      length(vi_params$a_t) == Fg &&
      sum(vi_params$a_t > 0) == Fg
    if (!check) stop("[] Incompatible initial value supplied for a_t")
  }
  if (is.null(vi_params[["b_t"]])) {
    vi_params$b_t <- numeric(Fg)
    for (f in 1:Fg) {
      # initialise these parameters using VI updates
      vi_params$b_t[f] <- hyper_fixed$b[f] + sum(1 - vi_params$prob[levels == f])
    }
  } else {
    check <- is.numeric(vi_params$b_t) &&
      length(vi_params$b_t) == Fg &&
      sum(vi_params$b_t > 0) == Fg
    if (!check) stop("[] Incompatible initial value supplied for b_t")
  }

  ## initialize for the local variational parameters psi and phi; they are
  # termed hyperparameter here, not because they are hyperparameters, but because
  # they are updated in the same function as as the hypeparameteters, tau_1_t and tau_2_t.
  # but more frequently.
  # NB: currently this is following moretrees, with the local variational parameters updated
  #     by something close to real update.

  # psi: v ijk: expected linear predictor squared.
  if (is.null(hyperparams[["psi"]])){
    zeta <- mapply(FUN = function(prob,mu) prob*mu,
                   prob = vi_params$prob,mu=vi_params$mu_gamma,SIMPLIFY=F)
    beta_v <- array(NA,c(pL,J,K)) #beta_map is directly estimated from the data as initialization; while zeta has prob.
    for (v in 1:pL){beta_v[v,,] <- Reduce('+',zeta[ancestors[[v]]])}
    hyperparams$psi <- abs(beta_v) # not exactly the update VI, but supposedly close: 1st moments.
  }  else {
    check <- is.numeric(hyperparams$psi) && sum(dim(hyperparams$psi)==c(pL,J,K))==3
    if(!check) stop("[] Incompatible initial values for 'psi' - local variational parameters")
  }
  if (random_init){
    hyperparams$psi <- abs(hyperparams$psi*
                             (1+rnorm(pL*J*K)*random_init_vals$psi_sd_frac))
  }
  hyperparams$g_psi <- g_fun.vec(hyperparams$psi)

  # phi: v km
  if (is.null(hyperparams[["phi"]])){ # pL by K-1
    xi <- mapply(FUN = function(prob,mu) c(prob*mu),
                 prob = vi_params$prob,mu=vi_params$mu_alpha,SIMPLIFY=F)
    hyperparams$phi <- matrix(0,nrow=pL,ncol=K-1)
    for (v in 1:pL){
      eta_v <- Reduce('+',xi[ancestors[[v]]])
      #print(eta_v)
      hyperparams$phi[v,] <- abs(eta_v)
    }
  }  else {
    check <- is.numeric(hyperparams$phi) &&
      sum(dim(hyperparams$phi)==c(pL,K-1))==2
    if(!check) stop("[] Incompatible initial values for 'phi' - local variational parameters")
  }
  if (random_init){
    hyperparams$phi <- abs(hyperparams$phi*
                             (1+rnorm(pL*(K-1))*random_init_vals$phi_sd_frac))
  }
  hyperparams$g_phi <- g_fun.vec(hyperparams$phi)

  # Initialize the multinomial variational parameters: n by K.
  if (is.null(vi_params[["rmat"]])){

    # fullLCM <- blca(Y,K,method="vb",verbose=FALSE)
    # vi_params$rmat <- fullLCM$Z

    Z <- unMAP(sample(1:K, size = n, replace = TRUE))
    if (ncol(Z) < K){
      Z <- cbind(Z, matrix(0, nrow = n, ncol = (K -ncol(Z))))
    }
    if (!is.null(Z_obs)){
      for (i in which(!is.na(Z_obs[,2]))){
        Z[i,Z_obs[i,2]] <- 1
        Z[i,-Z_obs[i,2]] <- 0
      }
    }
    vi_params$rmat <- Z
    # print(vi_params$rmat[!is.na(Z_obs[,2]),])

    # zeta <- mapply(FUN = function(prob,mu) prob*mu,
    #                prob = vi_params$prob,mu=vi_params$mu_gamma,SIMPLIFY=F)
    # lp <- array(0,c(n,J,K))
    # for (v in 1:pL){
    #   beta_v <- Reduce('+',zeta[ancestors[[v]]])
    #   for (i in leaf_ids_units[[v]]){
    #     lp[i,,] <- diag(X[i,,drop=FALSE])%*%beta_v
    #   }
    # }
    #
    # # calculate the eta_vk - eta_vm matrix:
    # xi <- mapply(FUN = function(prob,mu) c(prob*mu,0),
    #              prob = vi_params$prob,mu=vi_params$mu_alpha,SIMPLIFY=F)
    # lp_diff <- array(0,c(pL,K,K))
    # for (v in 1:pL){
    #   eta_v <- Reduce('+',xi[ancestors[[v]]])
    #   for (k in 1:K){
    #     for (m in (1:K)[-k]){
    #       lp_diff[v,k,m] <- eta_v[k] - eta_v[m]
    #     }
    #   }
    # }
    #
    # # calculate the multinomial variational probabilities,
    # # with the quadratic term set to zero as approximation; in actual VI update
    # # we will account for the quadratic term.
    # for (v in 1:pL){
    #   units <- leaf_ids_units[[v]]
    #   tmp   <- sum_m_neq_k(log(expit(hyperparams$phi[v,,]))+
    #                          (lp_diff[v,,]-hyperparams$phi[v,,])/2)
    #   tmp2 <- log(expit(hyperparams$psi[rep(v,length(units)),,,drop=FALSE]))+
    #     (lp[units,,,drop=FALSE]-hyperparams$psi[rep(v,length(units)),,,drop=FALSE])/2
    #   vi_params$rmat[units,] <- exp(apply(tmp2,c(1,3),sum)+t(replicate(length(units),tmp)))
    # }
    # vi_params$rmat <- sweep(vi_params$rmat,MARGIN=1,FUN="/",STATS=rowSums(vi_params$rmat))
  }


  # Initialize the variance parameters, sigma_gamma (for each u,j,k), Sigma_alpha(for each u):
  if (shared_tau){ # SHARED TAU ------------------->
    if (is.null(vi_params[["tau_1_t"]])){ # just read in initial tau_1's:
      vi_params$tau_1_t <- hyperparams$tau_1[levels]
    } else{
      check <- is.numeric(vi_params$tau_1_t) &&
        length(vi_params$tau_1_t) == p &&
        sum(vi_params$tau_1_t>0) ==p
      if (!check) stop("Incompatible intial values for 'tau_1_t'; for alpha")
    }

    if (is.null(vi_params[["tau_2_t"]])){
      vi_params$tau_2_t <- hyperparams$tau_2[levels]
    } else{
      check <- is.numeric(vi_params$tau_2_t) &&
        length(vi_params$tau_2_t) == p &&
        sum(vi_params$tau_2_t>0) ==p
      if (!check) stop("Incompatible intial values for 'tau_2_t'; for gamma")
    }

  } else{ # DISTINCT TAU ----------------------->
    if (is.null(vi_params[["tau_1_t"]])){ # just read in initial tau_1's:
      vi_params$tau_1_t <- split_along_dim(hyperparams$tau_1[levels,,drop=FALSE],1)
    } else{
      check <- length(vi_params$tau_1_t) == p &&
        sum(unlist(lapply(vi_params$tau_1_t,function(v) sum(v>0)==K-1))) == p# <--- NB: need fixing!
      if (!check) stop("Incompatible intial values for 'tau_1_t'; for alpha")
    }

    if (is.null(vi_params[["tau_2_t"]])){
      vi_params$tau_2_t <- split_along_dim(hyperparams$tau_2[levels,,,drop=FALSE],1)
    } else{
      check <- length(vi_params$tau_2_t) == p &&
        sum(unlist(lapply(vi_params$tau_2_t,function(v) sum(v>0)==K*J))) == p
      if (!check) stop("Incompatible intial values for 'tau_2_t'; for gamma")
    }
  }

  if (is.null(vi_params[["sigma_gamma"]])){ # depends on multinomial parameters.
    vi_params$sigma_gamma<- array(NA,c(p,J,K))
    for (u in 1:p){
      leaf_list_tmp <- leaf_ids_units[leaf_ids_nodes[[u]]]
      units         <- unlist(leaf_list_tmp)
      v_units_curr       <- unlist(mapply(rep,leaf_ids_nodes[[u]],unlist(lapply(leaf_list_tmp,length))))
      if (shared_tau){ # SHARED TAU:
        for (j in 1:J){
            vi_params$sigma_gamma[u,j,] <- 1/(2*colSums(vi_params$rmat[units,,drop=FALSE] *
                                                          hyperparams$g_psi[v_units_curr,j,])+1/(vi_params$tau_2_t[u]*h_pau[u])) # <-- zero h_pau?
        }
      } else{ # DISTINCT TAU:
        for (j in 1:J){
            vi_params$sigma_gamma[u,j,] <- 1/(2*colSums(vi_params$rmat[units,,drop=FALSE] *
                                                          hyperparams$g_psi[v_units_curr,j,])+1/(vi_params$tau_2_t[[u]][j,]*h_pau[u])) # <-- zero h_pau?
        }
      }
    }
  } else{
    check <- sum(dim(vi_params$sigma_gamma)==c(p,J,K))==3 &&
      sum(vi_params$sigma_gamma>0) == p*J*K
    if (!check) stop("[] Incompatible initial value for sigma_gamma.")
  }

  if (is.null(vi_params[["Sigma_alpha"]])){
    vi_params$Sigma_alpha <- list()
    # for (u in 1:p){
    #   leaf_list_tmp <- leaf_ids_units[leaf_ids_nodes[[u]]]
    #   units <- unlist(leaf_list_tmp)
    #   v_units <- c(unlist(mapply(rep,leaf_ids_nodes[[u]],unlist(lapply(leaf_list_tmp,length)))))
    #   tmp <- matrix(0,nrow=K,ncol=K)
    #   for (k in 1:K){
    #     for (i in seq_along(units)){
    #       Gvk <- diag(hyperparams$g_phi[v_units[i],k,-k])
    #       tmp <- tmp + 2*vi_params$rmat[i,k]*t(D_k(k,K))%*%Gvk%*%D_k(k,K)
    #     }
    #   }
    #   if (h_pau[u]>0){
    #     vi_params$Sigma_alpha[[u]] <- solve(tmp[-K,-K]+diag(1/(vi_params$tau_1_t[[u]]*h_pau[u]))) # <-- zero h_pau?
    #   } else{
    #     vi_params$Sigma_alpha[[u]] <- diag(rep(0,K-1))
    #   }
    # }

    if (shared_tau){ # SHARED_TAU:
      for (u in 1:p){
        vi_params$Sigma_alpha[[u]] <- c(1/getC(u,hyperparams$g_phi,vi_params$rmat,
                                               vi_params$tau_1_t,h_pau,subject_id_list[[u]],v_units))
      }
    } else{ # DISTINCT TAU:
      for (u in 1:p){
        leaf_list_tmp <- leaf_ids_units[leaf_ids_nodes[[u]]]
        #units         <- unlist(leaf_list_tmp)
        #v_units       <- unlist(mapply(rep,leaf_ids_nodes[[u]],unlist(lapply(leaf_list_tmp,length))))
        vi_params$Sigma_alpha[[u]] <- c(1/getC_separate_tau(u,hyperparams$g_phi,vi_params$rmat,
                                                            do.call("rbind",vi_params$tau_1_t),h_pau,subject_id_list[[u]],v_units))
      }
    }

    #print(vi_params$Sigma_alpha)
  }else{
    check <- is.list(vi_params$Sigma_alpha) &&
      #sum(sapply(vi_params$Sigma_alpha,is.vector)) == p &&
      sum(sapply(vi_params$Sigma_alpha,function(x) {length(x)==K-1}))==p
    if (!check) stop("[] Incompatible initial value for Sigma_alpha. (variational variance for alpha_u)")
  }

  # ELBO:
  if (is.null(hyperparams$ELBO)){
    hyperparams$ELBO  <- 1E-16
  } else{
    hyperparams$ELBO <- hyperparams$ELBO[length(hyperparams$ELBO)]
  }

  make_list(vi_params,hyperparams)
}





