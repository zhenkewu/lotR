#######################################################################
# functions for initialization:
# 1. obtain approximate posterior mean, variance and credible
#    intervals at various levels
# 2. compute some group-specific estimates, where the groups are
#   obtained from the posterior collapsing.
#######################################################################

#' compute model summariies from the model outputs
#'
#' @param mod output from \code{\link{lcm_tree}}
#' @param dsgn output from \code{\link{design_tree}}
#' @param ci_level credible interval level; default to 0.95
#'
#' @return a list
#' \describe{
#'
#' Named("beta_est")=beta_est,
#' Named("beta_sd_est")=beta_sd_est,
#' Named("beta_cil")=beta_cil,
#' Named("beta_ciu")=beta_ciu,
#' Named("eta_est")=eta_est,
#' Named("eta_var_est")=eta_var_est,
#' Named("eta_cil")=eta_cil,
#' Named("eta_ciu")=eta_ciu
#'
#'  prob_est$pi
#'  prob_est$theta
#'  prob_est$group
#'  prob_est$members
#'  prob_est$n_obs
#'  prob_est$theta_collapsed
#'  prob_est$pi_collapsed
#' }
#'
#' @seealso \code{\link{get_est_cpp}}
#' @export
compute_params <- function(mod,dsgn,ci_level=0.95){
  # # <---- delete when fdinishing testing.
  # mod = mod0$mod
  # dsgn= mod0$dsgn
  # # <----- delete when finishing testing

  leaf_ids    <- names(dsgn$leaf_ids)
  prob        <- mod$vi_params$prob
  mu_gamma    <- mod$vi_params$mu_gamma
  mu_alpha    <- mod$vi_params$mu_alpha
  sigma_gamma <- mod$vi_params$sigma_gamma
  Sigma_alpha <- mod$vi_params$Sigma_alpha

  anc     <- dsgn$ancestors
  cardanc <- unlist(lapply(anc,length))

  p <- length(prob)
  J <- ncol(dsgn$Y)
  K <- mod$hyper_fixed$K
  node_select <- (prob >0.5)

  z <- stats::qnorm(ci_level+(1-ci_level)/2)
  prob_est <- get_est_cpp(node_select,
              array(unlist(mu_gamma),c(J,K,p)),
              aperm(sigma_gamma,c(2,3,1)),
              as.matrix(do.call("rbind",mu_alpha)),
              as.matrix(do.call("rbind",Sigma_alpha)),
              anc,cardanc,z)

  prob_est$pi <- t(apply(prob_est$eta_est,1,function(u){tsb(c(expit(u),1))}))
  prob_est$theta <- expit(prob_est$beta_est)

  #prob_est <- mod0$prob_est
  grp <- round(prob_est$beta_est[2,1,],6) # <-- this is currently quite ad hoc.
  prob_est$group <- as.integer(factor(grp,levels=unique(grp)))

  prob_est$members <- vector("list",length(unique(prob_est$group)))
  for (i in 1:length(unique(prob_est$group))){
      prob_est$members[[i]] <- rownames(dsgn$A_leaves)[prob_est$group==i]
  }

  prob_est$n_obs <- sapply(prob_est$members,
                           function(u) sum(leaf_ids %in% u))
  # output estimates from collapsed groups:
  prob_est$theta_collapsed <- prob_est$theta[,,!duplicated(prob_est$pi)]
  prob_est$pi_collapsed <- prob_est$pi[!duplicated(prob_est$pi),]
  prob_est
}
