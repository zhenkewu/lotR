#' update variational parameters (core internal function)
#'
#' NB: needs speed up for multiplications over multiple arrays/lists of matrices.
#'
#' @param Y,A,outcomes_units,outcomes_nodes,h_pau,ancestors,levels outputs from
#' \code{design_tree()}
#' @param X,n,J,p,pL,Fg,K, computed from the data by \code{lcm_tree()}
#' @param  prob,mu_gamma,mu_alpha,rmat,sigma_gamma,Sigma_alpha,tau_1_t,tau_2_t,a_t,b_t,
#' variational parameters updated by \code{update_vi_params()}
#' @param psi,g_psi,phi,g_phi,tau_1,tau_2, parameters updated by \code{update_hyperparams()}
#' @param a,b,K fixed hyperparameters
#'
#' @importFrom matrixStats logSumExp
#'
#' @family Internal VI functions

# prob,mu,Sigma,Sigma_inv,Sigma_det,tau_t,delta,Omega,Omega_inv,Omega_det,a_t,b_t
# xxT,wwT computed from \code{X} and \code{W} in \code{spike_and_slab_logisitic_moretrees()}
update_vi_params <- function(Y,A,
                             outcomes_units,
                             outcomes_nodes,
                             ancestors,cardanc,
                             v_units,
                             h_pau,
                             levels,
                             X,n,J,p,pL,Fg,
                             #cardleaf,card_leaf_desc,maxnv,
                             subject_id_list,
                             # data and design
                             prob,mu_gamma,mu_alpha,rmat,
                             sigma_gamma,Sigma_alpha,# Sigma_alpha is list of length p; in Rcpp this is matrix p by K-1.
                             tau_1_t, # tau_1_t: list of p, each being K-1 vector;  R cpp: matrix of p by K-1.
                             tau_2_t, # tau_2_t: array (J,K,p); each being a J by K matrix
                             # - in Rcpp, array; in R, a list of p (J by K) matrices
                             a_t,b_t, # vi parameters
                             psi,g_psi,phi,g_phi,
                             tau_1,tau_2, # hyperparams
                             a,b,K,s_u_zeroset,s_u_oneset
                             #fixed hyper-parameters not to be update.
){
  X <- as.matrix(X)
  # iterate over internal + leaf nodes' (tau_1_t, tau_2_t, gamma_u, alpha_u, s_u):
  # Note: only after the entire block is updated we have a monotonic increase in ELBO!
  seq_update <- 1:p
  if (!is.null(s_u_zeroset)){seq_update <- (1:p)[-s_u_zeroset];prob[s_u_zeroset] <- 0}
  if (!is.null(s_u_oneset)){prob[s_u_oneset] <- 1}

  if(sum(abs(2*Y-1-X))!=0){stop("[data issue X and Y not matched.]")}
  # # calculate initial moments that are required in the VI updates:
  moments_cpp <- get_moments_cpp(prob,array(unlist(mu_gamma),c(J,K,p)),
                                 aperm(sigma_gamma,c(2,3,1)),
                                 as.matrix(do.call("rbind",mu_alpha)),
                                 as.matrix(do.call("rbind",Sigma_alpha)),
                                 ancestors,cardanc)

  # needed in update_hyperparams:
  E_beta_sq      <- moments_cpp$E_beta_sq
  E_eta_sq       <- moments_cpp$E_eta_sq
  E_beta         <- moments_cpp$E_beta
  E_eta          <- moments_cpp$E_eta


  rmat <- update_rZ(psi,g_psi,phi,g_phi,as.matrix(X),
                    E_beta,E_eta,E_beta_sq,E_eta_sq,
                    v_units)

  for (u in seq_update){
    tau_1_t[u] <- tau_1[levels[u]] # tau_1 is of legnth Fg
    tau_2_t[u] <- tau_2[levels[u]] # tau_2 is of length Fg

    gamma_update <- update_gamma_subid(u,g_psi,rmat, #<------- consider combining alpha and gamma update. move rmat updat to the last.
                                       tau_2_t[u],h_pau,levels,
                                       E_beta,as.matrix(prob[u]*mu_gamma[[u]],nrow=J,ncol=K),
                                       as.matrix(X),
                                       subject_id_list[[u]],v_units)

    # load gamma moment updates to mu_gamma (list, 265; J by K mat), sigma_gamma (array,p,J,K):
    mu_gamma[[u]]    <- gamma_update$resB*gamma_update$resA #  J by K
    sigma_gamma[u,,] <- gamma_update$resA

    alpha_update <- update_alpha_subid(u,g_phi,rmat,
                                       tau_1_t[u],h_pau,levels,
                                       E_eta, c(prob[u]*mu_alpha[[u]]),
                                       subject_id_list[[u]],v_units)

    # load alpha moment updates to mu_alpha (list, 265, vecs), Sigma_alpha (list, 265, K-1,sq mat):
    mu_alpha[[u]]    <-  c(alpha_update$resD)*c(alpha_update$resC)
    Sigma_alpha[[u]] <-  c(alpha_update$resC)

    # update w_u:
    if (u %in% s_u_zeroset){
      prob[u] <- 0
    } else if (u %in% s_u_oneset){
      prob[u] <- 1
    }else{
      w_u   <- digamma(a_t[levels[u]])-digamma(b_t[levels[u]])-
        0.5*J*K*log(tau_2_t[u]*h_pau[u])+0.5*sum(log(sigma_gamma[u,,]))+
        0.5*exp(logSumExp(c(gamma_update$resBsqC)))-
        0.5*(K-1)*log(tau_1_t[u]*h_pau[u])+0.5*sum(log(Sigma_alpha[[u]]))+
        0.5*exp(logSumExp(c(alpha_update$resDsqC)))
      prob[u] <- expit(w_u)
    }


    # update quantities derived from current moments:
    # ('skinny' for not computing sq terms.) this is needed because updating an
    # ancestral node will impact the current
    # moments of the descendants. (also prob changes everything)
    # if (u<p){
    # moments_cpp <- get_moments_cpp_skinny(prob,
    #                                       array(unlist(mu_gamma),c(J,K,p)),
    #                                       as.matrix(do.call("rbind",mu_alpha)),
    #                                       ancestors,
    #                                       cardanc)
    moments_cpp <- get_moments_cpp_eco(outcomes_nodes[[u]],
                                       E_beta,E_beta_sq,E_eta,E_eta_sq,
                                       prob,
                                       array(unlist(mu_gamma),c(J,K,p)),
                                       aperm(sigma_gamma,c(2,3,1)),
                                       as.matrix(do.call("rbind",mu_alpha)),
                                       as.matrix(do.call("rbind",Sigma_alpha)),
                                       ancestors,cardanc)
    E_beta_sq      <- moments_cpp$E_beta_sq
    E_eta_sq       <- moments_cpp$E_eta_sq
    E_beta         <- moments_cpp$E_beta
    E_eta          <- moments_cpp$E_eta

    # }
  }
  # # because iterating over u, it involves moments that may have just been updated,
  # # we need to recalculate moments at the end of the loop:
  # moments_cpp <- get_moments_cpp(prob,array(unlist(mu_gamma),c(J,K,p)),
  #                                aperm(sigma_gamma,c(2,3,1)),
  #                                as.matrix(do.call("rbind",mu_alpha)),
  #                                as.matrix(do.call("rbind",Sigma_alpha)),
  #                                ancestors,cardanc)
  #
  # # needed in update_hyperparams:
  # E_beta_sq      <- moments_cpp$E_beta_sq
  # E_eta_sq       <- moments_cpp$E_eta_sq
  # E_beta         <- moments_cpp$E_beta
  # E_eta          <- moments_cpp$E_eta


  # update rho:
  for (l in 1:Fg){
    a_t[l] <- a[l] + sum(prob[levels == l])
    b_t[l] <- b[l] + sum(1-prob[levels == l])
  }

  make_list(a_t,b_t,rmat,tau_1_t,tau_2_t,sigma_gamma,mu_gamma,Sigma_alpha,mu_alpha,prob,
            E_beta_sq,E_eta_sq,E_beta,E_eta)
  # NB needed in calculating the ELBO - this is calculated here because rmat updates need this.
}

#' update hyperparameters
#'
#' @inheritParams update_vi_params
#' @param update_hyper Logical, \code{TRUE} or \code{FALSE} to indicate
#'
#' @importFrom matrixStats logSumExp
#'
#' @return a list of updated hyperparameters: tau_1,tau_2,psi,g_psi,phi,g_phi,
#' along with a new ELBO value.
update_hyperparams <- function(Y,A,
                               outcomes_units,
                               ancestors,
                               h_pau,
                               levels,
                               v_units,
                               X,n,J,p,pL,Fg, # redundant but streamlined in lcm_tree.
                               # data and design
                               prob,mu_gamma,mu_alpha,rmat,
                               sigma_gamma,Sigma_alpha,
                               tau_1_t,tau_2_t, # <-- update tau_1, tau_2.(this is unique values); tau_1_t is a full vector.
                               a_t,b_t,
                               E_beta_sq,E_eta_sq,E_beta,E_eta,# vi parameters; !!! although this can be calculated from mu's sigma's
                               psi,g_psi,phi,g_phi,
                               tau_1,tau_2, # hyper-params # <----- need to be passed around for hyperparameter calculations.
                               a,b, K,tau_update_levels,#fixed hyper-parameters not to be update.
                               update_hyper#,
                               #cardanc,D_k_array,idnotk,submat
){
  # calculate some useful quantities:
  if(sum(abs(2*Y-1-X))!=0){stop("[data issue X and Y not matched.]")}
  # update local vi parameters (called hyper-parameters),
  psi   <- aperm(sqrt(E_beta_sq),c(3,1,2))
  phi   <- sqrt(E_eta_sq)

  g_psi <- g_fun.vec(psi)
  g_phi <- g_fun.vec(phi) # these can be moved to updating VI. In the model fitting, this is run at every iteration.

  ## calculations for updating hyperparameters:
  # update hyper-parameters: tau_1, tau_2:
  expected_ss_alpha <- expected_ss_gamma <- numeric(p) # marginal expectation of squared (alpha or gamma):
  for (u in 1:p){
    expected_ss_alpha[u] <- 1/h_pau[u]*(prob[u]*(#exp(logSumExp(c(-alpha_update$logresC,alpha_update$resDsqC-alpha_update$logresC))))+
      exp(logSumExp(c(log(Sigma_alpha[[u]]),2*log(abs(mu_alpha[[u]])))))      )+ #<--- why they are the same, for diag Sigma_alpha[[u]]
        (1-prob[u])*(K-1)*tau_1_t[u]*h_pau[u])
    expected_ss_gamma[u] <- 1/h_pau[u]*(prob[u]*(#exp(logSumExp(c(-gamma_update$logresA,c(gamma_update$resBsqA-gamma_update$logresA)))))+
      exp(logSumExp(c(c(log(sigma_gamma[u,,])),c(2*log(abs(mu_gamma[[u]]))))))        )+
        (1-prob[u])*J*K*tau_2_t[u]*h_pau[u])
  }
  # update hyper-parameters: tau_1, tau_2:
  # expected_ss_alpha <- expected_ss_gamma <- numeric(p) # marginal expectation of squared (alpha or gamma):
  # for (u in 1:p){
  #   expected_ss_alpha[u] <- 1/h_pau[u]*(prob[u]*(sum((Sigma_alpha[[u]]))+sum(mu_alpha[[u]]^2))+ #<--- why they are the same, for diag Sigma_alpha[[u]]
  #                                         (1-prob[u])*(K-1)*tau_1_t[u]*h_pau[u])
  #   expected_ss_gamma[u] <- 1/h_pau[u]*(prob[u]*(sum(sigma_gamma[u,,])+sum(mu_gamma[[u]]^2))+
  #                                         (1-prob[u])*J*K*tau_2_t[u]*h_pau[u])
  # }


  if (update_hyper){
    for (l in 1:Fg){
      if (l %in% tau_update_levels){
        tau_1[l]  <- sum(expected_ss_alpha[levels==l])/((K-1)*sum(levels==l))
        tau_2[l]  <- sum(expected_ss_gamma[levels==l])/(J*K*sum(levels==l))
        cat("\n [] update tau_1/2[",l,"]:",  tau_1[l],", ",tau_2[l],"; \n")
      }
    }
  }

  # update ELBO:
  expected_l_rho    <- digamma(a_t) - digamma(a_t+b_t)
  expected_l_1m_rho <- digamma(b_t) - digamma(a_t+b_t)

  # E_q(lower bound of joint distribution (all data and unknowns)):
  res1_2_13  <- get_line1_2_subid(psi,g_psi,phi,g_phi,rmat,E_beta,E_beta_sq,E_eta,E_eta_sq,as.matrix(X),v_units)
  line_1     <- res1_2_13$res1
  line_2     <- res1_2_13$res2
  line_3     <- sum(expected_l_rho[levels]*prob+expected_l_1m_rho[levels]*(1-prob))
  line_4     <- sum((a-1)*expected_l_rho+(b-1)*expected_l_1m_rho-mapply(lbeta,a,b))
  line_5     <- - sum(expected_ss_alpha/2/tau_1[levels])-sum((K-1)*log(2*pi*tau_1[levels]*h_pau))/2
  line_6     <- - sum(expected_ss_gamma/2/tau_2[levels])-sum(J*K*log(2*pi*tau_2[levels]*h_pau))/2 # this is prior, use hyerprior params.

  # -E_q(log(q)):
  line_7  <- (J*K*sum(prob)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob[u]*log(sigma_gamma[u,,])))))/2
  line_8  <- J*K*sum(1-prob)/2 +J*K*sum(log(2*pi*tau_2_t*h_pau)*(1-prob))/2
  line_9  <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line_10 <- ((K-1)*sum(prob)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob[u]*log(Sigma_alpha[[u]])))))/2 ############ <------------------
  line_11 <-((K-1) / 2) * sum(1 - prob) + ((K-1) / 2) * sum(log(2 * pi * tau_1_t*h_pau) * (1 - prob))
  line_12 <- -1* sum((a_t - 1) * expected_l_rho + (b_t - 1) * expected_l_1m_rho - mapply(lbeta, a_t, b_t)) # a_t, b_t is of L dimension.
  line_13 <- res1_2_13$res3
  ELBO    <- line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 +
    line_8 + line_9 +line_10 + line_11 + line_12+ line_13

  # print(make_list(line_1))
  # print(make_list(line_2))
  # print(make_list(line_3))
  # print(make_list(line_4))
  # print(make_list(line_5))
  # print(make_list(line_6))
  # print(make_list(line_7))
  # print(make_list(line_8))
  # print(make_list(line_9))
  # print(make_list(line_10))
  # print(make_list(line_11))
  # print(make_list(line_12))
  # print(make_list(line_13))
  # return results:
  make_list(ELBO,psi,g_psi,phi,g_phi,tau_1,tau_2)
}
