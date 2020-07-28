#' fit the tree-shrinkage latent class models (this function largely follows \code{moretrees})
#'
#' Returns the variational parameters (mostly moments of relevant unknown quantities)
#'
#' @param dat_tree a list of data and other information organized according to the tree
#' @param vi_params_init,hyperparams_init,random_init,random_init_vals,tol,tol_hyper,max_iter,print_freq,update_hyper_freq,hyper_fixed
#' initial values and updating protocols. Explained more in the wrapper function.
#' @return a list with model updates
#'
#' \describe{
#' \item{\code{vi_params}}{named list of final variational parameter estimates}
#' \item{\code{hyperparams}}{named list of final hyperparameter estimates}
#' \item{\code{hyper_fixed}}{named list of fixed hyperparameters}
#' \item{\code{ELBO_track}}{numeric vector containing the values of the objective function
#' (ELBO) at the end of every iteration}
#' }
#'
#' @family internal VI functions
fit_lcm_tree <- function(dsgn,
                         vi_params_init,
                         hyperparams_init,
                         random_init,
                         random_init_vals,
                         tol,
                         tol_hyper,
                         max_iter,
                         print_freq,
                         update_hyper_freq,
                         hyper_fixed
){
  ## set up the design list: 'dsgn' so this can be passed to initialization or updating functions:
  ## but when doing initialization, we want parsimony, so the dsgn does not
  ## have all the dimension information that can be calculated; Here, we first
  ## add elements to 'dsgn', although this is redundant, but can avoid calculating
  ## these dimensions again in vi_params and hyperparams updates; this means the functional
  ## arguments for vi_params and hyperparams updates contains more arguments, e.g.,
  ## J, K, than what's minimal.
  ##
  ## NB: consider adding a class to the dsgn - so people know what to look for?
  ##     check advanced R.

  #-------------------------------BEGIN DESIGN PADDING--------------------------
  dsgn$X  <- 2*dsgn$Y-1
  dsgn$n  <- nrow(dsgn$Y)
  dsgn$J  <- ncol(dsgn$Y)
  dsgn$p  <- length(unique(unlist(dsgn$ancestors)))
  dsgn$pL <- length(dsgn$ancestors)
  dsgn$Fg <- max(dsgn$levels)
  if (is.null(hyper_fixed$K)){stop("[lotR] # of classes 'K' not specified.")}
  if (!is.null(hyper_fixed$K) && hyper_fixed$K<2){stop("[lotR] # of classes 'K' is 1.")}
  K  <- hyper_fixed$K

  # for updating moments: do it here to avoid computing repeatedly in update_vi:----
  dsgn$cardleaf       <- unlist(lapply(dsgn$outcomes_units,length)) # number of observations in each leaf.
  dsgn$card_leaf_desc <- unlist(lapply(dsgn$outcomes_nodes,length)) # number of leaf descendants for each node.
  dsgn$cardanc        <- unlist(lapply(dsgn$ancestors,length))      # number of ancestors for each leaf.
  #dsgn$maxnv          <- max(unlist(lapply(dsgn$outcomes_units,length)))

  ## for updating alpha-moments:
  #-------------------------------END OF DESIGN PADDING--------------------------

  cat("\n [lotR]branch lengths: h_pau \n")
  print(dsgn$h_pau)
  # initialize: ----------------------
  init <- R.utils::doCall(initialize_tree_lcm,
                          vi_params   = vi_params_init,
                          hyperparams = hyperparams_init,
                          hyper_fixed = hyper_fixed,
                          random_init = random_init,
                          random_init_vals = random_init_vals,
                          args = dsgn)

  vi_params   <- init$vi_params
  hyperparams <- init$hyperparams

  cat("--- Iteration", 0, "; ELBO = ", init$hyperparams$ELBO,"\n") # needed if convergence across i -> i+1 is very quick before 2*update_hyper_freq/

  # initialize ELBO:
  ELBO_track <- numeric(max_iter)

  # run algorithm: ---------------------
  i <- 0
  repeat{

    # iterate i
    i <- i + 1


    # check if max_iter reached:
    if (i > max_iter){
      i <- max_iter
      cat(paste("--- Iteration", i, "complete. \n"))
      warning("[lotR] Maximum number of iterations reached! Consider increasing 'max_iter'")
      break
    }

    # update vi params:
    vi_params <- R.utils::doCall(update_vi_params,
                                 args = c(dsgn, vi_params, hyperparams, hyper_fixed))

    # This is directly borrowed from 'moretrees' ------------------------------
    # compute ELBO and update psi, phi and hyperparameters (tau_1, tau_2):
    update_hyper <- i %% update_hyper_freq == 0
    hyperparams  <- R.utils::doCall(update_hyperparams,
                                    update_hyper = update_hyper,
                                    args = c(dsgn,vi_params,hyperparams,hyper_fixed)
    )

    ELBO_track[i] <- hyperparams$ELBO
    cat("--- Iteration", i, "; epsilon = ", ELBO_track[i] - ELBO_track[i-1], "; ELBO = ", ELBO_track[i],"\n")

    # print progress:
    if (i %% print_freq ==0 && i>3){
      if(ELBO_track[i] - ELBO_track[i-1]<0){
        cat("--- Iteration", i, "; epsilon = ", ELBO_track[i] - ELBO_track[i-1], "; ELBO = ", ELBO_track[i],"\n")
        cat("> empirical class probabilities: ", round(colMeans(vi_params$rmat),4),"\n")
        cat("> node_select: ",which(vi_params$prob>0.5),"\n")
      }
      barplot(vi_params$prob)
      image(expit(vi_params$mu_gamma[[1]])) # root node.
    }

    # check tolerance
    if (update_hyper & i >= 2 * update_hyper_freq) {
      # if we just updated hyperparameters, check for convergence of hyperparameters
      criterion1 <- abs(ELBO_track[i] - ELBO_track[i - update_hyper_freq]) < tol_hyper
      if (criterion1) {
        # did last VI update reach convergence?
        criterion2 <- abs(ELBO_track[i - 1] - ELBO_track[i - 2]) < tol
        # if yes, both have converged. if not, continue.
        if (criterion2) break else next
      } else next
    } else {
      criterion3 <- (i > 2) && (abs(ELBO_track[i] - ELBO_track[i - 1]) < tol)
      # if criterion 3, fill in results until just before the
      # next hyperparameter update (or max_iter, whichever comes first)
      if (criterion3) { # is the current update good enough?
        if (i<2*update_hyper_freq){ # if update_hyper, but not yet 2*update_hyper_freq:
          # criterion4 <- (abs(ELBO_track[i] - init$hyperparams$ELBO) < tol_hyper) |
          #   (abs(ELBO_track[i] - ELBO_track[1]) < tol_hyper)
          # if (criterion4)
          break
        }
        i2 <- min(ceiling(i / update_hyper_freq) * update_hyper_freq - 1,
                  max_iter)
        ELBO_track[(i + 1):i2] <- hyperparams$ELBO  # can send this iteration much later; so appears updating more frequent than specified.
        i <- i2
      }
    }
  } # end 'repeat' loop.

  # return results:
  c(make_list(vi_params, hyperparams, hyper_fixed),
    list(ELBO_track=ELBO_track[1:i]))
}

#' estimate latent class models by pre-specified groups of observations,
#' e.g., the ones obtained from the tree-shrunk LCM.
lcm_by_group <- function(Y,outcomes,members_list,K,...){
  G <- length(unique(members_list))
  res <- vector("list",G)
  for (g in 1:G){
    units <- which(outcomes %in% members_list[[g]])
    #print(length(units))
    res[[g]] <- BayesLCA::blca(Y[units,,drop=FALSE],K,verbose=FALSE,...)
  }
  res
}









