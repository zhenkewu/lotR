#' lotR:  \strong{l}atent class analysis of  \strong{o}bservations
#' organized by  \strong{t}ree in  \strong{R}
#'
#' \code{lotR} is designed for analyzing multivariate binary observations
#' while integrating additional sample related information represented
#' by each observation's membership in the leaves of a given tree. The
#' observations that are closer in the tree are a priori more likely to be grouped
#' together and fitted by a LCM with identical LCM parameters. The model
#' is built on spike-and-slab priors on the increments of a Gaussian diffusion process
#' for each node of the tree. The model is self-adaptive in that it automatically
#' choose the optimal grouping of observations to fit distinct latent
#' class models. The posterior inferential algorithm is based on variational inference
#' and can provide approximate posterior uncertainty quantification.
#'
#' @seealso
#' \itemize{
#' \item \url{https://github.com/zhenkewu/lotR} for the source code
#' and system/software requirements to use \code{lotR} for your data.
#' }
#'
#' @section main lotR wrapper function:
#' \code{\link{lcm_tree}}
#'
#' @docType package
#' @name lotR
NULL
#> NULL
