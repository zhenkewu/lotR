#' approximate predictive class probabilities for a new observation in
#' a leaf
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
predict.lcm_tree <- function(object,
                             dat_pred=NULL,leaf_ids=NULL,collapsed=TRUE,
                             ...){
  # ## test
  # x <- mod0
  # dat_pred <- simu_dat_pred
  # leaf_ids <- simu_pred$curr_leaves
  # ## end of test

  N_pred <- nrow(dat_pred)
  K <- ncol(object$prob_est$pi_collapsed)

  postprob <- function(pi_vec,theta_mat,Y){
    K <- length(pi_vec)
    res <- rep(NA,K)
    for (k in 1:K){
      res[k] <- sum(Y*log(theta_mat[,k])+(1-Y)*log(1-theta_mat[,k]))+log(pi_vec[k])
    }
    exp(res - logSumExp(res))
  }

  probpred_mat <- matrix(NA,nrow=N_pred,ncol=K)
  if (collapsed){
    leaf_grps <- object$prob_est$members
    G <- length(leaf_grps)

    # get the group id for each observation's leaf id:
    group_ids_units <- sapply(leaf_ids,function(id) {
      (1:G)[unlist(lapply(leaf_grps,function(curr_grp) id%in%curr_grp))]})

    #names(object$dsgn$leaf_ids_units)
    for (i in 1:nrow(dat_pred)){
      curr_pi_vec <- object$prob_est$pi_collapsed[group_ids_units[i],]
      probpred_mat[i,] <- postprob(curr_pi_vec,object$prob_est$theta_collapsed,dat_pred[i,])
    }
  }else{
    use_pi_id <- sapply(leaf_ids,function(id)which(names(object$dsgn$leaf_ids_units)==id))
    for (i in 1:nrow(dat_pred)){
      curr_pi_vec <- object$prob_est_indiv$pi[use_pi_id[i],]
      probpred_mat[i,] <- postprob(curr_pi_vec,object$prob_est$theta_collapsed,dat_pred[i,])
    }
  }
  probpred_mat
}

# x <- mod0
# dat_pred <- simu_dat_pred
# leaf_ids <- simu_pred$curr_leaves
#
# aa = predict.lcm_tree(mod0,simu_dat_pred,simu_pred$curr_leaves,FALSE)
# bb = predict.lcm_tree(mod0,simu_dat_pred,simu_pred$curr_leaves,TRUE)
# table(apply(aa,1,which.max),apply(bb,1,which.max))


#' approximate predictive class probabilities for a new observation
#'
#' @param object An `lcm_tree` class object; Output from `BayesLCA::blca()`
#' @param dat_pred a binary matrix with rows for observations and columns
#' for features (of the same dimension as in the training data)
#' @param ... Other parameters
#'
#' @return a matrix with the same number of rows as `dat_pred` and
#' `K` columns (`K` is the number of latent classes);
#' each row correspond to a vector of predicted probabilities of
#' an observation in each class.
#'
#' @importFrom matrixStats logSumExp
#' @export
predict.blca <- function(object,
                         dat_pred=NULL,
                         ...){
  # # ## test
  # object <- res
  # dat_pred <- simu_dat_pred
  # # leaf_ids <- simu_pred$curr_leaves
  # # ## end of test
  N_pred <- nrow(dat_pred)
  K <- nrow(object$itemprob)

  postprob <- function(pi_vec,theta_mat,Y){
    K <- length(pi_vec)
    res <- rep(NA,K)
    for (k in 1:K){
      res[k] <- sum(Y*log(theta_mat[,k])+(1-Y)*log(1-theta_mat[,k]))+log(pi_vec[k])
    }
    exp(res - logSumExp(res))
  }

  probpred_mat <- matrix(NA,nrow=N_pred,ncol=K)
  curr_pi_vec  <- object$classprob
  for (i in 1:nrow(dat_pred)){
    probpred_mat[i,] <- postprob(curr_pi_vec,t(object$itemprob),dat_pred[i,])
  }
  probpred_mat
}

# x <- mod0
# dat_pred <- simu_dat_pred
# leaf_ids <- simu_pred$curr_leaves
#
# aa = predict.lcm_tree(mod0,simu_dat_pred,simu_pred$curr_leaves,FALSE)
# bb = predict.lcm_tree(mod0,simu_dat_pred,simu_pred$curr_leaves,TRUE)
# table(apply(aa,1,which.max),apply(bb,1,which.max))


# x <- mod0
# dat_pred <- simu_dat_pred
# leaf_ids <- simu_pred$curr_leaves
#
# aa = predict.lcm_tree(mod0,simu_dat_pred,simu_pred$curr_leaves,FALSE)
# bb = predict.lcm_tree(mod0,simu_dat_pred,simu_pred$curr_leaves,TRUE)
# table(apply(aa,1,which.max),apply(bb,1,which.max))

