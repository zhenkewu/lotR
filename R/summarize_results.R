#######################################################################
# functions for initialization:
# 1. obtain approximate posterior mean, variance and credible
#    intervals at various levels
# 2. compute some group-specific estimates, where the groups are
#   obtained from the posterior collapsing.
# 3. print the model results: print.lcm_tree(); this can be compact or
#    very detailed.
#######################################################################

#' compute model summaries from the model outputs
#'
#' This function is used in the wrapper function [lcm_tree()]
#'
#' @param mod output from [fit_lcm_tree()]
#' @param dsgn output from [design_tree()]
#' @param ci_level credible interval level; default to 0.95
#' @param B Default is `10000`; the number of random multivariate
#' Gaussian eta (K-1 dimensional) for producing posterior summaries
#' of the stick-breaking transformed probabilities (K dimensional)
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
#' @seealso [get_est_cpp()]
#' @export
compute_params <- function(mod,dsgn,ci_level=0.95,B=10000){
  # <---- delete when fdinishing testing.
  # mod = mod0$mod
  # dsgn= mod0$dsgn
  # ci_level=0.95
  # <----- delete when finishing testing

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
  # the following are assuming a single set of class-specific response probability
  # profiles:
  prob_est$theta <- expit(prob_est$beta_est)[,,1]
  prob_est$beta_est <- prob_est$beta_est[,,1]
  prob_est$beta_sd_est <- prob_est$beta_sd_est[,,1]
  prob_est$beta_cil <- prob_est$beta_cil[,,1]
  prob_est$beta_ciu <- prob_est$beta_ciu[,,1]

  #prob_est <- mod0$prob_est
  grp <- round(prob_est$eta_est[,1],6) # <-- this is currently quite ad hoc.
  prob_est$group <- as.integer(factor(grp,levels=unique(grp)))
  # this is likely better because it is not on probability scale, which
  # reflects true diffusion on the eta scale.

  prob_est$members <- vector("list",length(unique(prob_est$group)))
  for (i in 1:length(unique(prob_est$group))){
    prob_est$members[[i]] <- rownames(dsgn$A_leaves)[prob_est$group==i]
  }

  prob_est$n_obs <- sapply(prob_est$members,
                           function(u) sum(leaf_ids %in% u))
  # output estimates from collapsed groups:
  prob_est$theta_collapsed <- prob_est$theta
  # prob_est$theta_collapsed <- prob_est$theta[,,!duplicated(prob_est$pi)]
  prob_est$pi_collapsed <- prob_est$pi[!duplicated(prob_est$group),,drop=FALSE]
  prob_est$eta_est_collapsed <- prob_est$eta_est[!duplicated(prob_est$group),,drop=FALSE]
  prob_est$eta_sd_est_collapsed <- prob_est$eta_sd_est[!duplicated(prob_est$group),,drop=FALSE]
  prob_est$ci_level <- ci_level

  ##################################################################################
  ## individual node estimates (no grouping): no posterior median model selection
  ##################################################################################
  prob_est_indiv <- get_est_cpp(prob,
                                array(unlist(mu_gamma),c(J,K,p)),
                                aperm(sigma_gamma,c(2,3,1)),
                                as.matrix(do.call("rbind",mu_alpha)),
                                as.matrix(do.call("rbind",Sigma_alpha)),
                                anc,cardanc,z)

  prob_est_indiv$pi <- t(apply(prob_est_indiv$eta_est,1,function(u){tsb(c(expit(u),1))}))
  # the following are assuming a single set of class-specific response probability
  # profiles:
  prob_est_indiv$theta <- expit(prob_est_indiv$beta_est)[,,1]
  prob_est_indiv$beta_est <- prob_est_indiv$beta_est[,,1]
  prob_est_indiv$beta_sd_est <- prob_est_indiv$beta_sd_est[,,1]
  prob_est_indiv$beta_cil <- prob_est_indiv$beta_cil[,,1]
  prob_est_indiv$beta_ciu <- prob_est_indiv$beta_ciu[,,1]

  prob_est_indiv$ci_level <- ci_level

  rslt <- vector("list",nrow(prob_est_indiv$eta_est))
  ## get credible intervals on the pie scales:
  for (v in 1:nrow(prob_est_indiv$eta_est)){
    tmp_post_simu <- apply(cbind(expit(MASS::mvrnorm(B,prob_est_indiv$eta_est[v,],diag((prob_est_indiv$eta_sd_est[v,])^2,K-1,K-1))),1),1,tsb)
    ci_mat <- apply(tmp_post_simu,1,stats::quantile,c((1-prob_est_indiv$ci_level)/2,prob_est_indiv$ci_level+(1-prob_est_indiv$ci_level)/2))

    rslt[[v]] <- matrix(NA,nrow=K,ncol=3)
    colnames(rslt[[v]]) <- c("est_pi","cil_pi","ciu_pi")
    for (k in 1:K){# for each class:
      rslt[[v]][k,] <- cbind(prob_est_indiv$pi[v,k],ci_mat[1,k],ci_mat[2,k])
    }
  }
  prob_est_indiv$rslt <- rslt

  ##################################################################################
  ## ///END///individual node estimates (no grouping): no posterior median model selection
  ##################################################################################

  make_list(prob_est,prob_est_indiv)
}


#' `print.lcm_tree` summarizes the results from [lcm_tree()].
#'
#' @param x Output from [lcm_tree()].
#' @param ... Arguments passed to summary and printing methods.
#' @return Summary showing, for each group of leaf nodes discovered by [lotR()],
#' the class prevalences, 95% credible intervals, number of leaf nodes in
#' per group, and number of observations per group.
#'
#' @family lcm_tree results
#' @export
print.lcm_tree <- function(x, ...){
  print(summary(x, ...), ...)
  # Return
  return(invisible(x))
}


#' `summary.lcm_tree` summarizes the results from [lcm_tree()].
#'
#' Will have some randomness associated with the upper and lower bounds,
#' because we simulated the Gaussian variables before converting to probabilities.
#'
#' @param object Output from [lcm_tree()].
#' An object of class "lcm_tree".# coeff_type Either "lcm_tree" or "ad_hoc"
#' @param compact If `TRUE`, a more compact summary of results is printed.
#' Only works well when the dimension of the variables is low
#' (say, < 4); otherwise the table takes up too much horizontal space.
#' Default is `FALSE`.
#' @param B Default `10000`: The number of random multivariate Gaussian for logit of V in stick breaking
#' parameterization; used for obtaining the posterior distribution of `pi_v`
#' @param ... Not used.
#' @return see [print.lcm_tree()]
#'
#' @family lcm_tree results
#' @export
summary.lcm_tree <- function(object,
                             compact = FALSE,
                             B = 10000,
                             ...) {
  # # <----- temporary during package building.
  # object = mod0
  # coeff_type = "lcm_tree"
  # compact = FALSE
  # # <---- temporary during package building.

  coeff_type = "lcm_tree"
  # Get estimates
  K  <-  ncol(object$prob_est$pi)
  pL <-  nrow(object$prob_est$pi)
  p  <-  length(object$mod$vi_params$tau_1_t)

  if (coeff_type == "lcm_tree") {
    coeff_type2 <- "prob"
  } else {
    coeff_type2 <- coeff_type
  }

  est_type <- paste0(coeff_type2,"_est")
  est <- object[est_type][[1]]

  est$n_leaves   <- sapply(object$prob_est$members, length)
  est$leaves     <- sapply(object$prob_est$members, paste0, collapse = ", ")
  est$n_obs      <- object$prob_est$n_obs
  est$group      <- 1:length(object$prob_est$n_obs)
  # est$z    <- stats::qnorm(est$ci_level+(1-est$ci_level)/2) # so need to put ci_level into ad_hoc_est.

  if (compact) {
    if (coeff_type=="lcm_tree"){
      # simulate the lower and upper intervals in probability scales.
      tmp_main <- matrix(NA,nrow=length(est$group),ncol=3*K)
      for (g in seq_along(est$group)){
        tmp_post_simu <- apply(cbind(expit(MASS::mvrnorm(B,est$eta_est_collapsed[g,],diag((est$eta_sd_est_collapsed[g,])^2,K-1,K-1))),1),1,tsb)
        ci_mat <- apply(tmp_post_simu,1,stats::quantile,c((1-est$ci_level)/2,est$ci_level+(1-est$ci_level)/2))
        for (k in 1:K){
          tmp_main[g,(3*k-2):(3*k)] <- cbind(est$pi_collapsed[g,k],ci_mat[1,k],ci_mat[2,k])
          #rowMeans(tmp_post_simu)
          #apply(tmp_post_simu,1,sd)
        }
      }
      rslt <- cbind(est$group,est$n_leaves,est$n_obs,tmp_main)
      colnames(rslt) <- c("group","n_leaves","n_obs",paste0(c("est_pi", "cil_pi", "ciu_pi"), rep(1:K, each = 3)))
      rslt <- data.frame(rslt)
      rslt$leaves <- est$leaves

      rslt <- list(est = rslt, coeff_type = coeff_type, ci_level = est$ci_level)
    } else{
      stop("[lotR] did not implement other coef_type.")
    }
    rslt$theta <- est$theta_collapsed
    class(rslt) <- "summary.lcm_tree_compact"
    return(rslt)
  }

  # make separate data.frames per group (this is better for >4
  # number of classes)
  grps <- list()
  for (g in seq_along(est$group)){
    grps[[g]] <- list(n_leaves = est$n_leaves[g],
                      n_obs = est$n_obs[g],
                      leaves = est$leaves[g])

    tmp_post_simu <- apply(cbind(expit(MASS::mvrnorm(B,est$eta_est_collapsed[g,],diag((est$eta_sd_est_collapsed[g,])^2,K-1,K-1))),1),1,tsb)
    ci_mat <- apply(tmp_post_simu,1,stats::quantile,c((1-est$ci_level)/2,est$ci_level+(1-est$ci_level)/2))
    grps[[g]]$est <-t( rbind(est$pi_collapsed[g,],ci_mat))
    colnames(grps[[g]]$est) <- c("est_pi","cil_pi","ciu_pi")
  }

  names(grps) <- paste0("Group", 1:length(grps))
  grps$coeff_type <- coeff_type
  grps$ci_level   <- est$ci_level
  class(grps) <- "summary.lcm_tree_long"
  # Return
  return(grps)
}


#' Compact printing of [lcm_tree()] model fits
#'
#' `print.summary.lcm_tree_compact` is a print method for class
#' `summary.lcm_tree_compact`.
#'
#' @param x output from `summary.lcm_tree` with `compact = TRUE`.
#' @param print_leaves If `TRUE`, for each discovered group the full list
#' of leaves will be printed. Set this to `FALSE` if these leaf lists
#' make output difficult to read.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @return see [print.lcm_tree()]
#'
#' @export
#' @family lcm_tree results
print.summary.lcm_tree_compact <- function(x,
                                           print_leaves = TRUE,
                                           digits = max(3L, getOption("digits") - 3L),
                                           ...) {
  if (!print_leaves) x$est$leaves <- NULL

  if (x$coeff_type == "lcm_tree") {
    cat("Showing latent class model estimates for discovered groups; not ad hoc.\n\n")
  }

  cat(paste0("Group-specific estimate(s), with credible interval level: ", x$ci_level))
  print(knitr::kable(x$est, digits = digits))

  cat("\n")

  # Return
  return(invisible(x))
}

#' Long-form printing of [lcm_tree()] model fits
#'
#' `print.summary.lcm_tree_compact` is a print method for class
#' `summary.lcm_tree_compact`.
#'
#' @param x output from `summary.lcm_tree` with `compact = FALSE`.
#' @param print_leaves If `TRUE`, for each discovered group the full list
#' of leaves will be printed. Set this to `FALSE` if these leaf lists
#' make output difficult to read.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @return see [print.lcm_tree()]
#'
#' @export
#' @family lcm_tree results
print.summary.lcm_tree_long <- function(x,
                                        print_leaves = TRUE,
                                        digits = max(3L, getOption("digits") - 3L),
                                        ...) {

  cat("Group-specific class prevalence estimates for the", length(x), "groups discovered by lcm_tree\n")
  cat(paste0("Credible interval level: ", x$ci_level),"\n")
  if (x$coeff_type == "lcm_tree") {
    cat("Showing latent class model estimates for discovered groups; not ad hoc.\n\n")
  }
  # remove some information that is not to be tabulated:
  x$coeff_type <- NULL
  x$ci_level  <- NULL

  for (g in 1:length(x)) {

    est <- x[[g]]
    cat("------------------ Group", g, "------------------\n\n")

    cat("Number of leaves:", est$n_leaves, "\n")
    cat("Number of observations:", est$n_obs, "\n")
    if (print_leaves) {
      cat("List of leaf nodes:\n")
      cat(est$leaves, "\n\n")
    } else {
      cat("\n\n")
    }

    cat("Class prevalence estimate(s):")
    print(knitr::kable(est$est, digits = digits, row.names = FALSE))

    cat("\n")
  }

  cat("\nIf this is hard to read, try print(x, compact = TRUE)\n\n")

  # Return
  return(invisible(x))
}


#' log likelihood for lca
#'
#' @param dat a binary matrix with rows for observations and columns
#' for features (of the same dimension as in the training data)
#' @param theta_mat a `J` by `K` matrix of positive response probabilities
#' @param pi_vec a `K`-vector that sum to one; latent class proportions.
#'
#' @return log likelihood
#'
#' @importFrom matrixStats logSumExp
#' @export
lca_ll <- function(dat,theta_mat,pi_vec){
  K <- length(pi_vec)
  ll <- 0
  res <- rep(NA,K)
  for (i in 1:nrow(dat)){
    for (k in 1:K){
      res[k] <- sum(dat[i,]*log(theta_mat[,k])+(1-dat[i,])*log(1-theta_mat[,k]))+log(pi_vec[k])
    }
    ll = ll+ logSumExp(res)
  }
  ll
}


#' log data likelihood for tree LCM
#'
#' @param object An `lcm_tree` class object; Output from `lcm_tree()`
#' @param dat_pred a binary matrix with rows for observations and columns
#' for features (of the same dimension as in the training data)
#' @param leaf_ids a vector of character strings, each being the name
#' of the leaf name for that observation
#' @param collapsed default to `TRUE`; if `FALSE`, no posterior median
#' model selection is done - so each leaf may have its own vector of
#' latent class proportions.
#' @param ... Other parameters
#'
#' @return a matrix with the same number of rows as `dat_pred` and
#' `K` columns (`K` is the number of latent classes);
#' each row correspond to a vector of predicted probabilities of
#' an observation in each class.
#'
#' @importFrom matrixStats logSumExp
#' @export
tree_lca_ll <- function(object,
                        dat_pred=NULL,leaf_ids=NULL,collapsed=TRUE,
                        ...){
  # ## test
  # object <- mod0
  # dat_pred <- Y_test
  # leaf_ids <- curr_leaves_test
  # collapsed=TRUE


  # ## end of test

  N_pred <- nrow(dat_pred)
  K <- ncol(object$prob_est$pi_collapsed)

  ll <- function(pi_vec,theta_mat,Y){
    K <- length(pi_vec)
    res <- rep(NA,K)
    for (k in 1:K){
      res[k] <- sum(Y*log(theta_mat[,k])+(1-Y)*log(1-theta_mat[,k]))+log(pi_vec[k])
    }
    logSumExp(res)
  }

  probpred_mat <- matrix(NA,nrow=N_pred,ncol=1)
  if (collapsed){
    leaf_grps <- object$prob_est$members
    G <- length(leaf_grps)

    # get the group id for each observation's leaf id:
    group_ids_units <- sapply(leaf_ids,function(id) {
      (1:G)[unlist(lapply(leaf_grps,function(curr_grp) id%in%curr_grp))]})

    #names(object$dsgn$leaf_ids_units)
    for (i in 1:nrow(dat_pred)){
      curr_pi_vec <- object$prob_est$pi_collapsed[group_ids_units[i],]
      probpred_mat[i,1] <- ll(curr_pi_vec,object$prob_est$theta_collapsed,dat_pred[i,])
    }
  }else{
    use_pi_id <- sapply(leaf_ids,function(id)which(names(object$dsgn$leaf_ids_units)==id))
    for (i in 1:nrow(dat_pred)){
      curr_pi_vec <- object$prob_est_indiv$pi[use_pi_id[i],]
      probpred_mat[i,1] <- ll(curr_pi_vec,object$prob_est$theta_collapsed,dat_pred[i,])
    }
  }
  sum(probpred_mat)
}





