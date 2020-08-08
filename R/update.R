#' update variational parameters in the latent class models with observations
#' organized in a tree
#'
#' Used by [fit_lcm_tree()], which also invoke [update_hyperparams()]
#' to update hyperparameters and calculate evidence lower bound (ELBO)
#'
#' @param Y,A,leaf_ids_units,leaf_ids_nodes,ancestors,cardanc,v_units,h_pau,levels,subject_id_list outputs from
#' [design_tree()], which reorders the nodes by internal and leaf nodes; the observations are also
#' ordered from low- to high-indexed leaf nodes.
#' @param X,n,J,p,pL,Fg computed from the data by [lcm_tree()] at the beginning before the VI updates
#' @param  prob,prob_gamma,mu_gamma,mu_alpha,rmat,sigma_gamma,Sigma_alpha,tau_1_t,tau_2_t,a_t,b_t,
#' variational parameters updated by [update_vi_params()]
#' @param psi,g_psi,phi,g_phi,tau_1,tau_2,shared_tau parameters updated by `update_hyperparams()`
#' @param a,b,K,s_u_zeroset,s_u_oneset,Z_obs fixed hyperparameters not to be updated.
#'
#' @importFrom matrixStats logSumExp
#' @export
#' @family Internal VI functions
update_vi_params <- function(Y,A,Z_obs,
                             leaf_ids_units,
                             leaf_ids_nodes,
                             ancestors,cardanc,
                             v_units,
                             h_pau,
                             levels,
                             subject_id_list,
                             X,n,J,p,pL,Fg,# data and design
                             prob,prob_gamma,mu_gamma,mu_alpha,rmat,
                             sigma_gamma,Sigma_alpha,# Sigma_alpha is list of length p; in Rcpp this is matrix p by K-1.
                             tau_1_t,
                             tau_2_t,
                             a_t,b_t,
                             psi,g_psi,phi,g_phi,
                             tau_1,tau_2,shared_tau,
                             a,b,K,s_u_zeroset,s_u_oneset#fixed hyper-parameters not to be updated.
){
  X <- as.matrix(X)
  if(sum(abs(2*Y-1-X))!=0){stop("[lotR] data issue X and Y not matched.")}
  # iterate over internal + leaf nodes' (tau_1_t, tau_2_t, gamma_u, alpha_u, s_u):
  # Note: only after the entire block is updated we have a monotonic increase in ELBO!
  seq_update <- 1:p
  if (!is.null(s_u_zeroset)){ # not updating nodes that are set to zeros.
    if (length(setdiff(s_u_zeroset,1:p))!=0){stop("[lotR] 's_u_zeroset has elements not between 1 and p.'")}
    seq_update <- (1:p)[-s_u_zeroset]
    prob[s_u_zeroset] <- 0
  }
  if (!is.null(s_u_oneset)){
    prob[s_u_oneset] <- 1
    if (length(setdiff(s_u_oneset,1:p))!=0){stop("[lotR] 's_u_oneset has elements not between 1 and p.'")}
  }

  if (!exists("E_beta_sq") || !exists("E_eta_sq") || !exists("E_beta") || !exists("E_eta")){
    # calculate initial moments that are required in the VI updates (do so only when not available):
    moments_cpp <- get_moments_cpp(prob,prob_gamma,array(unlist(mu_gamma),c(J,K,p)),
                                   aperm(sigma_gamma,c(2,3,1)),
                                   as.matrix(do.call("rbind",mu_alpha)),
                                   as.matrix(do.call("rbind",Sigma_alpha)),
                                   ancestors,cardanc)
    # needed in updating rmat, and for update_hyperparams:
    E_beta_sq      <- moments_cpp$E_beta_sq
    E_eta_sq       <- moments_cpp$E_eta_sq
    E_beta         <- moments_cpp$E_beta
    E_eta          <- moments_cpp$E_eta
  }
  # update for all subjects their variational probabilities of belonging to each of K classes: q_t(Z^v_i)
  if (!is.null(Z_obs)){
    ind_unknown <- which(is.na(Z_obs[,2]))
    tmp  <- update_rmat_partial(ind_unknown,psi,g_psi,phi,g_phi,as.matrix(X),E_beta,E_eta,E_beta_sq,E_eta_sq,v_units)
    rmat[ind_unknown,] <- tmp[ind_unknown,] # don't update observed Z; it is initialized.
  } else{
    rmat  <- update_rmat(psi,g_psi,phi,g_phi,as.matrix(X),E_beta,E_eta,E_beta_sq,E_eta_sq,v_units)
  }

  for (u in seq_update){
    ## --- | update q_t(gamma_u,alpha_u,s_u):
    ## because it is a two-component mixture - two Gaussian components, with distinct
    #mean (e.g., 0 and mu_alpha[[u]])and variance parameters (e.g., tau_1_t[u]*h_pau[u] and Sigma_alpha[u,,])
    ##with component indicator s_u.
    # variational variance parameters for the s_u=0 component; they do not include h_pau[u].
    if (shared_tau){
      tau_1_t[u] <- tau_1[levels[u]] # tau_1 is of legnth Fg
      tau_2_t[u] <- tau_2[levels[u]] # tau_2 is of length Fg

      gamma_alpha_update <- update_gamma_alpha_subid(u,g_psi,g_phi,
                                                     tau_2_t[u],tau_1_t[u],
                                                     E_beta,as.matrix(prob_gamma[u]*mu_gamma[[u]],nrow=J,ncol=K),as.matrix(X),
                                                     E_eta, c(prob[u]*mu_alpha[[u]]),
                                                     rmat,h_pau,levels,subject_id_list[[u]],v_units)
      mu_gamma[[u]]    <-   gamma_alpha_update$resB*gamma_alpha_update$resA #  J by K
      sigma_gamma[u,,] <-   gamma_alpha_update$resA
      mu_alpha[[u]]    <-  c(gamma_alpha_update$resD)*c(gamma_alpha_update$resC)
      Sigma_alpha[[u]] <-  c(gamma_alpha_update$resC)

      # update w_u, q_t(s_u)
      w_u   <- digamma(a_t[levels[u]])-digamma(b_t[levels[u]])-
        # 0.5*J*K*log(tau_2_t[u]*h_pau[u])+0.5*sum(log(sigma_gamma[u,,]))+
        # 0.5*exp(logSumExp(c(gamma_alpha_update$logresBsq_o_A)))-
        0.5*(K-1)*log(tau_1_t[u]*h_pau[u])+0.5*sum(log(Sigma_alpha[[u]]))+
        0.5*exp(logSumExp(c(gamma_alpha_update$logresDsq_o_C)))
    } else{ # <--- if separate tau's.
      tau_1_t[[u]] <- tau_1[levels[u],] # tau_1 is Fg by K-1
      tau_2_t[[u]] <- tau_2[levels[u],,] # tau_2 is Fg by J by K.

      gamma_alpha_update <- update_gamma_alpha_subid_separate_tau(
        u,g_psi,g_phi,
        tau_2_t[[u]],tau_1_t[[u]],
        E_beta,as.matrix(prob_gamma[u]*mu_gamma[[u]],nrow=J,ncol=K),as.matrix(X),
        E_eta, c(prob[u]*mu_alpha[[u]]),
        rmat,h_pau,levels,subject_id_list[[u]],v_units)
      mu_gamma[[u]]    <-   gamma_alpha_update$resB*gamma_alpha_update$resA #  J by K
      sigma_gamma[u,,] <-   gamma_alpha_update$resA
      mu_alpha[[u]]    <-  c(gamma_alpha_update$resD)*c(gamma_alpha_update$resC)
      Sigma_alpha[[u]] <-  c(gamma_alpha_update$resC)

      # update w_u, q_t(s_u)
      w_u   <- digamma(a_t[levels[u]])-digamma(b_t[levels[u]])-
        # 0.5*sum(log(c(tau_2_t[[u]])*h_pau[u]))+0.5*sum(log(sigma_gamma[u,,]))+
        # 0.5*exp(logSumExp(c(gamma_alpha_update$logresBsq_o_A)))-
        0.5*sum(log(tau_1_t[[u]]*h_pau[u]))+0.5*sum(log(Sigma_alpha[[u]]))+
        0.5*exp(logSumExp(c(gamma_alpha_update$logresDsq_o_C)))
    }
    prob[u] <- expit(w_u)

    # updating ancestral node will impact the variational moments of the descendants:
    moments_cpp <- get_moments_cpp_eco(leaf_ids_nodes[[u]],
                                       E_beta,E_beta_sq,E_eta,E_eta_sq,
                                       prob,prob_gamma,array(unlist(mu_gamma),c(J,K,p)),
                                       aperm(sigma_gamma,c(2,3,1)),
                                       as.matrix(do.call("rbind",mu_alpha)),
                                       as.matrix(do.call("rbind",Sigma_alpha)),
                                       ancestors,cardanc)
    E_beta_sq      <- moments_cpp$E_beta_sq
    E_eta_sq       <- moments_cpp$E_eta_sq
    E_beta         <- moments_cpp$E_beta
    E_eta          <- moments_cpp$E_eta
  }

  # update variational parameters for q_t(rho):
  for (l in 1:Fg){
    a_t[l] <- a[l] + sum(prob[levels == l])
    b_t[l] <- b[l] + sum(1-prob[levels == l])
  }

  make_list(a_t,b_t,rmat,tau_1_t,tau_2_t,sigma_gamma,mu_gamma,Sigma_alpha,mu_alpha,
            prob,prob_gamma,
            E_beta_sq,E_eta_sq,E_beta,E_eta)
}

#' update hyperparameters in the latent class model with observations organized in
#' a tree
#'
#' @inheritParams update_vi_params
#' @param update_hyper Logical, `TRUE` or `FALSE` to indicate
#' whether to update `tau_1` and `tau_2`. This is computed at every iteration
#' in [fit_lcm_tree()]
#' @param E_beta_sq,E_eta_sq,E_beta,E_eta moments computed by [update_vi_params()]
#' @param tau_update_levels a numeric vector, specifies which levels of hyperparameters to update
#' @param quiet default to `FALSE`, which prints intermediate updates of hyperparameters
#'
#' @importFrom matrixStats logSumExp
#'
#' @return a list of updated hyperparameters: tau_1,tau_2,psi,g_psi,phi,g_phi,
#' along with a new ELBO value.
#' @export
update_hyperparams <- function(Y,A,
                               leaf_ids_units,
                               ancestors,
                               h_pau,
                               levels,
                               v_units,
                               X,n,J,p,pL,Fg, # redundant but streamlined in lcm_tree.
                               # data and design
                               prob,prob_gamma,mu_gamma,mu_alpha,rmat,
                               sigma_gamma,Sigma_alpha,
                               tau_1_t,tau_2_t, # <-- update tau_1, tau_2.(this is unique values); tau_1_t is a full vector.
                               a_t,b_t,
                               E_beta_sq,E_eta_sq,E_beta,E_eta,# vi parameters; !!! although this can be calculated from mu's sigma's
                               psi,g_psi,phi,g_phi,
                               tau_1,tau_2,shared_tau, # hyper-params
                               a,b, K,tau_update_levels,#fixed hyper-parameters not to be update.
                               update_hyper,quiet # called in 'fit_lcm' to update tau_1 and tau_2 or not.
){
  # calculate some useful quantities:
  if(sum(abs(2*Y-1-X))!=0){stop("[lotR] data issue X and Y not matched.")}
  # update local vi parameters, these can be moved to updating VI. In the model fitting, this is run at every iteration.
  psi   <- aperm(sqrt(E_beta_sq),c(3,1,2))
  phi   <- sqrt(E_eta_sq)

  g_psi <- g_fun.vec(psi)
  g_phi <- g_fun.vec(phi)

  if (shared_tau){
    # update hyper-parameters: tau_1, tau_2:
    expected_ss_alpha <- expected_ss_gamma <- numeric(p) # marginal variational posterior expectation of squared (alpha or gamma):
    for (u in 1:p){
      expected_ss_alpha[u] <- 1/h_pau[u]*(prob[u]*(
        exp(logSumExp(c(log(Sigma_alpha[[u]]),2*log(abs(mu_alpha[[u]])))))      )+
          (1-prob[u])*(K-1)*tau_1_t[u]*h_pau[u])
      expected_ss_gamma[u] <- 1/h_pau[u]*(prob_gamma[u]*(
        exp(logSumExp(c(c(log(sigma_gamma[u,,])),c(2*log(abs(mu_gamma[[u]]))))))        )+
          (1-prob_gamma[u])*J*K*tau_2_t[u]*h_pau[u])
    }

    if (update_hyper){ # computed at each iteration - TRUE to update the hyperparameters.
      for (l in 1:Fg){
        if (l %in% tau_update_levels){
          tau_1[l]  <- sum(expected_ss_alpha[levels==l])/((K-1)*sum(levels==l))
          tau_2[l]  <- sum(expected_ss_gamma[levels==l])/(J*K*sum(levels==l))
          # cat("> Updated tau_1 and tau_2; level ",l,":",  tau_1[l],", ",tau_2[l],". \n")
        }
      }
    }

    # update ELBO:
    expected_l_rho    <- digamma(a_t) - digamma(a_t+b_t)
    expected_l_1m_rho <- digamma(b_t) - digamma(a_t+b_t)

    # part 1: E_q(lower bound of joint distribution (all data and unknowns)):
    res1_2_13  <- get_line1_2_13_subid(psi,g_psi,phi,g_phi,rmat,E_beta,E_beta_sq,E_eta,E_eta_sq,as.matrix(X),v_units)
    line_1     <- res1_2_13$res1
    line_2     <- res1_2_13$res2
    line_3     <- sum(expected_l_rho[levels]*prob+expected_l_1m_rho[levels]*(1-prob))
    line_4     <- sum((a-1)*expected_l_rho+(b-1)*expected_l_1m_rho-mapply(lbeta,a,b))
    line_5     <- - sum(expected_ss_alpha/2/tau_1[levels])-sum((K-1)*log(2*pi*tau_1[levels]*h_pau))/2
    line_6     <- - sum(expected_ss_gamma/2/tau_2[levels])-sum(J*K*log(2*pi*tau_2[levels]*h_pau))/2 # this is prior, use hyperparams.

    # part 2: -E_{q_t}(log(q_t)):
    line_7     <- (J*K*sum(prob_gamma)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob_gamma[u]*log(sigma_gamma[u,,])))))/2
    line_8     <- J*K*sum(1-prob_gamma)/2 +J*K*sum(log(2*pi*tau_2_t*h_pau)*(1-prob_gamma))/2
    line_9     <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
    line_10    <- ((K-1)*sum(prob)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob[u]*log(Sigma_alpha[[u]])))))/2
    line_11    <-((K-1) / 2) * sum(1 - prob) + ((K-1) / 2) * sum(log(2 * pi * tau_1_t*h_pau) * (1 - prob))
    line_12    <- -1* sum((a_t - 1) * expected_l_rho + (b_t - 1) * expected_l_1m_rho - mapply(lbeta, a_t, b_t)) # a_t, b_t is of Fg dimension.
    line_13    <- res1_2_13$res3

  } else{ #<--------------- if distinct tau_1's in a node.
    # update hyper-parameters: tau_1, tau_2:
    expected_ss_alpha <- expected_ss_gamma <- vector("list",p) # marginal variational posterior expectation of squared (alpha or gamma):
    for (u in 1:p){
      expected_ss_alpha[[u]] <- (prob[u]*(c(Sigma_alpha[[u]])+c(mu_alpha[[u]]^2))+(1-prob[u])*c(tau_1_t[[u]])*h_pau[u])/h_pau[u]
      expected_ss_gamma[[u]] <- (prob_gamma[u]*(sigma_gamma[u,,]+mu_gamma[[u]]^2)+(1-prob_gamma[u])*tau_2_t[[u]]*h_pau[u])/h_pau[u]
    }

    if (update_hyper){ # computed at each iteration - TRUE to update the hyperparameters.
      for (l in 1:Fg){
        if (l %in% tau_update_levels){
          tau_1[l,]  <- colMeans(do.call("rbind",expected_ss_alpha[levels==l]))
          tau_2[l,,] <- apply(array(unlist(expected_ss_gamma[levels==l]),c(J,K,sum(levels==l))),c(1,2),mean)
          # if (!quiet){cat("> Updated tau_1 and tau_2; level ",l,":",  tau_1[l],", ",tau_2[l],". \n")}
        }
      }
    }

    # update ELBO:
    expected_l_rho    <- digamma(a_t) - digamma(a_t+b_t)
    expected_l_1m_rho <- digamma(b_t) - digamma(a_t+b_t)

    # part 1: E_q(lower bound of joint distribution (all data and unknowns)):
    res1_2_13  <- get_line1_2_13_subid(psi,g_psi,phi,g_phi,rmat,E_beta,E_beta_sq,E_eta,E_eta_sq,as.matrix(X),v_units)
    line_1     <- res1_2_13$res1
    line_2     <- res1_2_13$res2
    line_3     <- sum(expected_l_rho[levels]*prob+expected_l_1m_rho[levels]*(1-prob))
    line_4     <- sum((a-1)*expected_l_rho+(b-1)*expected_l_1m_rho-mapply(lbeta,a,b))
    line_5     <- - sum(sapply(1:p,function(u){sum(expected_ss_alpha[[u]]/2/tau_1[levels[u],])}))-0.5*(sum(sapply(1:p,function(u){sum(log(2*pi*tau_1[levels[u],]*h_pau[u]))})))
    line_6     <- - sum(sapply(1:p,function(u){sum(expected_ss_gamma[[u]]/2/tau_2[levels[u],,])}))-0.5*(sum(sapply(1:p,function(u){sum(log(2*pi*tau_2[levels[u],,]*h_pau[u]))}))) # this is prior, use hyerprior params.

    # part 2: -E_q(log(q)):
    line_7  <- (J*K*sum(prob_gamma)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob_gamma[u]*log(sigma_gamma[u,,])))))/2
    line_8  <- J*K*sum(1-prob_gamma)/2 + sum(mapply(FUN=function(pp,mat,hh){sum(pp*log(2*pi*mat*hh))},pp=1-prob_gamma,mat=tau_2_t,hh=h_pau))/2
    line_9  <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
    line_10 <- ((K-1)*sum(prob)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob[u]*log(Sigma_alpha[[u]])))))/2
    line_11 <- (K-1) / 2 * sum(1 - prob) + sum(mapply(FUN=function(pp,v,hh){sum(pp*log(2 * pi * v*hh) )},pp=1-prob,v=tau_1_t,hh=h_pau))/2
    line_12 <- -1* sum((a_t - 1) * expected_l_rho + (b_t - 1) * expected_l_1m_rho - mapply(lbeta, a_t, b_t)) # a_t, b_t is of L dimension.
    line_13 <- res1_2_13$res3
  }
  ELBO       <- line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 +
    line_8 + line_9 +line_10 + line_11 + line_12+ line_13

  # return results:
  make_list(ELBO,psi,g_psi,phi,g_phi,tau_1,tau_2)
}
