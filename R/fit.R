#' Fit latent class model with tree-structured shrinkage over the observations
#' (this function largely follows Thomas, E. \code{moretrees})
#'
#'
#' @param dsgn a list of data and other information organized according to the tree
#' @param vi_params_init,hyperparams_init,random_init,random_init_vals,tol,tol_hyper,max_iter,print_freq,quiet,plot_fig,shared_tau,update_hyper_freq,hyper_fixed
#' initial values and updating protocols. Explained more in the wrapper function \code{\link{lcm_tree}}
#' @return a list with model updates; because the variational posterior
#' is comprised of familiar distributional forms that can be  determined
#' by the moments, the returned values are moments:
#'
#' \describe{
#' \item{\code{vi_params}}{named list of final variational parameter estimates}
#' \item{\code{hyperparams}}{named list of final hyperparameter estimates}
#' \item{\code{hyper_fixed}}{named list of fixed hyperparameters}
#' \item{\code{ELBO_track}}{numeric vector containing the values of the objective function
#' (ELBO) at the end of every iteration}
#' }
#' @importFrom graphics barplot image
#' @family internal VI functions
#' @export
fit_lcm_tree <- function(dsgn,
                         vi_params_init,
                         hyperparams_init,
                         random_init,
                         random_init_vals,
                         tol,
                         tol_hyper,
                         max_iter,
                         print_freq,
                         quiet,
                         plot_fig,
                         shared_tau,
                         update_hyper_freq,
                         hyper_fixed
){
  ## set up the design list: 'dsgn' so this can be passed to initialization or updating functions.
  ## At initialization, we opted for parsimony, so `dsgn` does not
  ## have all the dimension information that can be calculated. We first
  ## add additional elements to 'dsgn'. Although this can be recalculated
  ## in vi_params and hyperparams updates. This means arguments of vi_params and hyperparams updates
  # contains more elements, e.g., J, K, than what can be minimal.
  ##
  ## NB: consider adding a class to the dsgn - so people know what to look for?
  ##     check advanced R.

  #-------------------------------BEGIN DESIGN PADDING--------------------------
  dsgn$X  <- 2*dsgn$Y-1 # this dsgn$Y would be different from the Y as input to design_tree
  dsgn$n  <- nrow(dsgn$Y)
  dsgn$J  <- ncol(dsgn$Y)
  dsgn$p  <- length(unique(unlist(dsgn$ancestors)))
  dsgn$pL <- length(dsgn$ancestors)
  dsgn$Fg <- max(dsgn$levels)
  if (is.null(hyper_fixed$K)){stop("[lotR] # of classes 'K' not specified.")}
  if (!is.null(hyper_fixed$K) && hyper_fixed$K<2){stop("[lotR] # of classes 'K' is 1.")}
  K  <- hyper_fixed$K

  #dsgn$cardleaf       <- unlist(lapply(dsgn$leaf_ids_units,length)) # number of observations in each leaf.
  #dsgn$card_leaf_desc <- unlist(lapply(dsgn$leaf_ids_nodes,length)) # number of leaf descendants for each node.
  dsgn$cardanc        <- unlist(lapply(dsgn$ancestors,length))      # number of ancestors for each leaf.
  #dsgn$maxnv          <- max(unlist(lapply(dsgn$leaf_ids_units,length)))
  #-------------------------------END OF DESIGN PADDING--------------------------

  if (!quiet){cat("\n [lotR] Branch lengths: `h_pau`: \n");print(dsgn$h_pau)}
  # initialize: ----------------------
  init <- R.utils::doCall(initialize_tree_lcm,
                          vi_params   = vi_params_init,
                          hyperparams = hyperparams_init,
                          hyper_fixed = hyper_fixed,
                          random_init = random_init,
                          random_init_vals = random_init_vals,
                          shared_tau = shared_tau,
                          args = c(dsgn))
  vi_params   <- init$vi_params
  hyperparams <- init$hyperparams
  cat("|--- Model Initialized.\n")

  # initialize ELBO:
  ELBO_track <- numeric(max_iter)

  #if (quiet){pb <- progress_bar$new(format = "(:spin) [:bar] :percent",
  #                                 total = max_iter, clear = FALSE)}
  #for (i in 1:30) { pb$tick() ;Sys.sleep(3 / 100) }

  # run algorithm: ---------------------
  i <- 0
  repeat{
    #if (quiet){pb$tick()} #;Sys.sleep(3 / 100)}
    # iterate i
    i <- i + 1


    # check if max_iter reached:
    if (i > max_iter){
      i <- max_iter
      cat(paste("|--- Iteration", i, "complete. \n"))
      warning("[lotR] Maximum number of iterations reached! Consider increasing 'max_iter'")
      break
    }

    # update vi params:
    vi_params <- R.utils::doCall(update_vi_params,
                                 shared_tau = shared_tau,
                                 args = c(dsgn, vi_params, hyperparams,
                                          hyper_fixed))

    # compute ELBO and update psi, phi and hyperparameters (tau_1, tau_2):
    update_hyper <- i %% update_hyper_freq == 0
    hyperparams  <- R.utils::doCall(update_hyperparams,
                                    update_hyper = update_hyper,
                                    shared_tau = shared_tau,
                                    quiet      = quiet,
                                    args = c(dsgn,vi_params,hyperparams,hyper_fixed))

    ELBO_track[i] <- hyperparams$ELBO

    # print progress:
    if (i %% print_freq ==0){
      #if(ELBO_track[i] - ELBO_track[i-1]<0){
      if (!quiet){
        cat("|--- Iteration", i, "; epsilon = ", ELBO_track[i] - ELBO_track[i-1], "; ELBO = ", ELBO_track[i],"\n")
        cat("> empirical class probabilities: ", round(colMeans(vi_params$rmat),4),"\n")
        cat("> node_select: ",which(vi_params$prob>0.5),"\n")
      }
      if (plot_fig){
        barplot(vi_params$prob)
        image(expit(vi_params$mu_gamma[[1]])) # root node.
      }
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
#'
#' @param Y data
#' @param leaf_ids individual's leaf indicators in a vector
#' @param members_list a list, of length G, which is the number of unique
#' groups of leaves
#' @param K the number of classes
#' @param ... other arguments for \code{\link[BayesLCA]{blca}}
#' @importFrom BayesLCA blca
#' @export
lcm_by_group <- function(Y,leaf_ids,members_list,K,...){
  G <- length(unique(members_list))
  res <- vector("list",G)
  for (g in 1:G){
    units <- which(leaf_ids %in% members_list[[g]])
    #print(length(units))
    res[[g]] <- BayesLCA::blca(Y[units,,drop=FALSE],K,verbose=FALSE,...)
  }
  res
}









