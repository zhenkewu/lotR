#' wrapper function for fitting and summaries
#'
#' @param Y An \code{n} x \code{J} matrix of binary (0 or 1) measurements,
#' where each of n observation may belong to a leaf node; we hope to group n subjects into
#' groups possibly coarser than the leaf-induced groups, and in each group, we have a homogeneous latent
#' class model
#' @param outcomes Character vector of length \code{n}. \code{outcomes[i]} is a string indicating the
#' outcome experienced by unit \code{i}.
#' @param mytree A directed \code{igraph} object. This is a tree representing the relationships
#' among the outcomes. The leaves represent individual outcomes, and internal nodes
#' represent outcome categories consisting of their leaf descendants. All nodes
#' of mytree must have unique names as given by \code{names(V(mytree))}. The names of the leaves must
#' be equal to the unique elements of outcomes. The vertices of \code{mytree}, \code{V(mytree)}, may have
#' an attribute \code{levels} containing integer values from 1 to \code{max(V(mytree)$levels)}.
#' In this case, the levels attribute specifies groups of nodes that share common
#' hyperparameters \code{rho[f]}, \code{tau_1[f]}, and \code{tau_2[f]}. If \code{V(mytree)$levels} is \code{NULL},
#' the default is two levels of hyperparameters: one for all leaf nodes, and one
#' for all internal nodes. NB: this needs to be checked in the \code{\link{design_tree}}
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, \code{ci_level = 0.95} (the default) returns a 95\% credible interval
#' @param get_lcm_by_group If \code{TRUE}, lotR will also return the maximum likelihood estimates of the
#' coefficients for each outcome group discovered by the model.
#' Default is \code{TRUE}.
#' @param update_hyper_freq How frequently to update hyperparameters.
#' Default = every 50 iterations.
#' @param print_freq How often to print out iteration number and current value of epsilon
#' (the difference in objective function value for the two most recent iterations).
#' @param hyper_fixed Fixed values of hyperprior parameters for rho.
#' This should be a list with two elements:
#' a and b, both numeric vectors of length \code{L}, representing the
#' parameters of the beta prior on rho for each level, where \code{L} is the
#' number of levels.
#' Default is \code{list(a = rep(1, L), b = rep(1, L), tau_update_levels = c(1,2))} (uniform hyperprior)
#' Other options include specifying includedd or excluded nodes via
#' e.g., \code{s_u_zeroset = (1:265)[-c(1)],s_u_oneset = c(1))} to fit a single
#' big LCM (collapsing all nodes except the root node, the first node)
#' @param tol Convergence tolerance for the objective function.
#' Default is \code{1E-8}.
#' @param tol_hyper The convergence tolerance for the objective function between
#' between subsequent hyperparmeter updates. Typically a more generous
#' tolerance than \code{tol}.
#' Default is \code{1E-4}.
#' @param max_iter Maximum number of iterations of the VI algorithm.
#' Default is 5000.
#' @param nrestarts Number of random re-starts of the VI algorithm. The result that
#' gives the highest value of the objective function will be returned.
#' It is recommended to choose \code{nrestarts > 1}.
#' The default is 3.
#' @param keep_restarts If \code{TRUE}, the results from all random restarts will be returned.
#' If \code{FALSE}, only the restart with the highest objective function is returned. '
#' Default is \code{TRUE}.
#' @param parallel If \code{TRUE}, the random restarts will be run in parallel.
#' It is recommended to first set the number of cores using \code{doParallel::registerDoParallel()}.
#' Otherwise, the default number of cores specified by the \code{doParallel} package will be used.
#' Default is \code{TRUE}.
#' @param log_restarts If \code{TRUE}, when \code{nrestarts > 1} progress of each random restart will be
#' logged to a text file in \code{log_dir}. If \code{FALSE} and \code{nrestarts > 1},
#' progress will not be shown.
#' If \code{nrestarts = 1}, progress will always be printed to the console.
#' Default is \code{FALSE}.
#' @param log_dir Directory for logging progress of random restarts.
#' Default is the working directory.
#' @param vi_params_init,hyperparams_init Named lists containing initial values for the
#' variational parameters and hyperparameters. Supplying good initial values can be challenging,
#' and \code{lotR()} provides a way to guess initial values based on transformations
#' of conditional logistic regression estimates of the effect sizes
#' for each individual outcome (see \code{\link{initialize_tree_lcm}}).
#' The most common use for \code{vi_params_init} and \code{hyperparams_init} is to supply starting
#' values based on previous output from \code{lotR()};
#' see the \code{vignette('lotR')} for examples.
#' The user can provide initial values for all parameters or a subset.
#' When initial values for one or more parameters are not
#' supplied, the missing values will be filled in by \code{\link{initialize_tree_lcm}}.
#' @param random_init
#' If \code{TRUE}, some random variability will be added to the initial values.
#' The default is \code{FALSE}, unless \code{nrestarts > 1}, in which case
#' \code{random_init} will be set to \code{TRUE} and a warning message will be printed.
#' The amount of variability is determined by \code{random_init_vals}.
#' @param random_init_vals If \code{random_init = TRUE},
#' this is a list containing the following parameters for randomly permuting
#' the inital values:
#' \describe{
#' \item{\code{tau_lims}}{a vector of length 2, where \code{tau_lims[1]} is between 0 and 1,
#' and \code{tau_lims[2] > 1}. The initial values for the hyperparameter \code{tau} will
#' be chosen uniformly at random in the range \code{(tau_init * tau_lims[1], tau_init * tau_lims[2])},
#' where \code{tau_init} is the initial value for \code{tau} either supplied in \code{hyperparams_init}
#' or guessed using \code{\link{initialize_tree_lcm}}.}
#' \item{\code{psi_sd_frac}}{a value between 0 and 1. The initial values for the auxilliary parameters
#' \code{psi} will have a normal random variate added to them with standard deviation equal to
#' \code{psi_sd_frac} multiplied by the initial value for eta either supplied in \code{hyperparams_init} or guessed
#' using \code{\link{initialize_tree_lcm}}. Absolute values are then taken for any
#' values of \code{psi} that are \code{< 0}.}
#' \item{\code{phi_sd_frac}}{same as above}.
#' \item{\code{mu_gamma_sd_frac}}{a value between 0 and 1. The initial values for
#' \code{mu} will have a normal random variate added to them with standard deviation equal to
#' \code{mu_sd_frac} multiplied by the absolute value of the initial value for \code{mu_gamma_sd_frac} either supplied in
#' \code{vi_params_init} or guessed using \code{\link{initialize_tree_lcm}}.}
#' \item{\code{mu_alpha_sd_frac}}{same as above.}
#' \item{\code{u_sd_frac}}{a value between 0 and 1. The initial value for the node inclusion probabilities
#' will first be transformed to the log odds scale to obtain \code{u}. A normal random variate will be
#' added to \code{u} with standard deviation eqaul to u_sd_frac multiplied by the absolute value of the
#' initial value for \code{u} either supplied in \code{vi_params_init} or guessed using \code{moretrees_init_logistic()}.
#' \code{u} will then be transformed back to the probability scale.}
#' }
#'
#' @param allow_continue logical, TRUE to save results so can continue running the VI
#' updates with the last iteration from the old results.
#'
#' @examples vignette('lotR')
#'
#' @return a list also of class "lcm_tree"
#'
#' \describe{
#'   res <- make_list(mod,mod_restarts,mytree,dsgn,prob_est,est_ad_hoc)
#'      class(res) <- c("lcm_tree","list")
#'    }
#'
#' @useDynLib lotR
#' @export
#'
#' @family lcm_tree functions
lcm_tree <- function(Y,outcomes,mytree,# may have unordered nodes.
                     rootnode = "Node1",
                     weighted_edge = TRUE,
                     ci_level = 0.95,
                     get_lcm_by_group = TRUE,
                     update_hyper_freq = 50,
                     print_freq = 10,
                     hyper_fixed = list(K=3), # <-- modify default?
                     tol = 1E-8,
                     tol_hyper = 1E-4,
                     max_iter = 5000,
                     nrestarts = 3,
                     keep_restarts = TRUE,
                     parallel = TRUE,
                     log_restarts = FALSE,
                     log_dir = ".",
                     vi_params_init = list(),
                     hyperparams_init = list(),
                     random_init = FALSE,
                     random_init_vals = list(tau_lims = c(0.5, 1.5),
                                             psi_sd_frac = 0.2,
                                             phi_sd_frac = 0.2,
                                             mu_gamma_sd_frac = 0.2,
                                             mu_alpha_sd_frac = 0.2,
                                             u_sd_frac = 0.2),
                     allow_continue = FALSE
){
  # compatibility checks (NB: not done yet as of July 28, 2020):

  # logs
  log_dir <- sub("/$", "", log_dir)
  if (log_restarts) message("[lotR] Algorithm progress for restart i will be printed to ",
                            log_dir, "/restart_i_log.txt\n", sep = "")

  # Fill in some arguments
  if (nrestarts > 1 & !random_init) {
    message("[lotR] Setting 'random_init = TRUE' since nrestarts > 1\n")
    random_init <- TRUE
  }
  if (nrestarts == 1) parallel <- FALSE

  # construct designed data; here design_tree is going to reorder the nodes of the tree.
  dsgn <- design_tree(Y,outcomes,mytree,rootnode,weighted_edge) # root_node,weighted_edge <--- need fixing.

  # Get hyper_fixed if not supplied:
  if (is.null(hyper_fixed$a) | is.null(hyper_fixed$b)) {
    L             <- max(dsgn$levels)
    hyper_fixed   <- append(hyper_fixed,list(a = rep(1, L)))
    hyper_fixed$b <- rep(1, L)
    warning("[lotR] No fixed hyperparameters supplied; we set a_l=b_l=1 for all levels of hyperparameters.")
  }

  if (is.null(hyper_fixed$K)) {
    warning("[lotR] No fixed # of classes supplied; supply a named element `K` in the list 'hyper_fixed'.")
  }

  # Setting up parallelization
  if (parallel) {
    `%doRestarts%` <- foreach::`%dopar%`
  } else {
    `%doRestarts%` <- foreach::`%do%`
  }

  # Run algorithm
  mod_restarts <- foreach::foreach(i = 1:nrestarts) %doRestarts% {
    if (log_restarts) {
      sink(file = paste0(log_dir, "/restart_", i, "_log.txt"))
    }
    cat("\nInitializing restart", i, "...\n\n")
    #cat("Random initialization: ", random_init,"...\n\n")
    mod <- fit_lcm_tree(dsgn = dsgn,
                        vi_params_init = vi_params_init,
                        hyperparams_init = hyperparams_init,
                        random_init = random_init,
                        random_init_vals = random_init_vals,
                        tol = tol,
                        tol_hyper = tol_hyper,
                        max_iter = max_iter,
                        print_freq = print_freq,
                        update_hyper_freq = update_hyper_freq,
                        hyper_fixed = hyper_fixed)
    cat("\nRestart", i, "complete.\n")
    if (log_restarts) {
      sink()
    }
    mod
  }

  # Select random restart that gave the highest ELBO
  ELBO_restarts <- sapply(mod_restarts, FUN = function(mod) mod$ELBO_track[length(mod$ELBO_track)])
  best_restart  <- which.max(ELBO_restarts)
  mod <- mod_restarts[[best_restart]]
  if (keep_restarts) {
    mod_restarts <- mod_restarts[- best_restart]
  } else {
    rm(mod_restarts)
    mod_restarts <- NULL
  }

  # compute the latent class prevalences and class-specific probabilities based
  # on the model output:
  # up to this point, we do not have dsgn in mod.
  prob_est <- compute_params(mod,dsgn,ci_level)

  # get a LCM result based on the groups formed by the output:
  if (!is.null(get_lcm_by_group) && get_lcm_by_group){
    est_ad_hoc <- lcm_by_group(dsgn$Y,names(dsgn$outcomes),
                               members_list = prob_est$members,
                               mod$hyper_fixed$K,method="em")
  }
  ## function for using a previous lcm_tree results to continue iterations;
  ## only useful when the previous result is not converged.
  if (allow_continue){ # here Y and outcomes may not be ordered....!!!
    old_mod <- append(make_list(Y,outcomes,mytree,rootnode,weighted_edge,#ci_level = 0.95#get_ml = TRUE
                                hyper_fixed,random_init_vals, # definitely fixed
                                update_hyper_freq,
                                print_freq # can use default, be better take the same values as here.
    ), #<---- here the Y here are not organized.
    list(vi_params_init   = mod$vi_params,
         hyperparams_init = mod$hyperparams))
    res <-   make_list(mod,   # <--- need to add other estimates.
                       mod_restarts,mytree,dsgn,old_mod)
    class(res) <- c("lcm_tree","list")
    return(res)
  }

  res <- make_list(mod,mod_restarts,mytree,dsgn, # is mytree redundant, check dsgn.
                   prob_est,est_ad_hoc)

  class(res) <- c("lcm_tree","list")
  res
}


#' continue from a prevoius model fits
#'
#' can definitely ignore this function, and use lcm_tree directly with various settings.
#' So this function is for lazy persons.
#'
#'
continue_lcm_tree <- function(old_mod,
                              #ci_level = 0.95 #get_ml = TRUE
                              update_hyper_freq = NULL,
                              print_freq = NULL,
                              tol =1e-8,
                              tol_hyper =1e-4,
                              max_iter = 5000,  # assume you want to finish in a second run; so set this long enough.
                              nrestarts = 1, # by default assume run once.
                              keep_restarts = TRUE,
                              parallel=TRUE,
                              log_restarts=FALSE,
                              log_dir= ".",
                              random_init=FALSE, # by default directly continue, so no random_init.
                              allow_continue = FALSE # by default not assuming continuing.
){

  if (is.null(old_mod)){stop("[] old model settings not saved; rerun with 'allow_continue=TRUE'.")}
  update_hyper_freq_curr <-update_hyper_freq
  print_freq_curr <- print_freq
  if (is.null(update_hyper_freq_curr)){
    update_hyper_freq_curr <-old_mod$update_hyper_freq
  }
  if (is.null(print_freq_curr)){print_freq_curr <-old_mod$print_freq}

  res <- R.utils::doCall(lcm_tree,
                         update_hyper_freq=update_hyper_freq_curr,
                         print_freq=print_freq_curr,
                         tol=tol,
                         tol_hyper=tol_hyper,
                         max_iter=max_iter,
                         nrestarts=nrestarts,
                         keep_restarts=keep_restarts,
                         parallel=parallel,
                         log_restarts=log_restarts,
                         log_dir=log_dir,
                         random_init=random_init,
                         allow_continue=allow_continue,
                         args = old_mod)
  class(res) <- c("lcm_tree","list")
  res
}




