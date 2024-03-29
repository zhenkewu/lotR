########################
# helper functions:
#######################

#' multi-expit (numerical overflow dealt with)
#'
#' @param v a vector of length K
#'
#' @return a vector of probabilities of length K; sums to one
#' @examples
#' sigmoid <- function(v) mexpit(c(0,v))[2]
#' mexpit(c(1,0,0))
#' @export
mexpit  <- function(v) exp(v - matrixStats::logSumExp(v))
#inverse of mexpit: alrInv()


#' sigmoid (copied from moretrees:)
#'
#' @param p a probability
#'
#' @export
#' @family utility function
logit   <- function(p) log(p) - log(1 - p)


#' log(1+exp(x)) (copied from moretrees:)
#'
#' @param x a number
#'
#' @export
#' @family utility function
log1p_exp <- function(x) {
  # computes log(1 + exp(x))
  if (x > 20) {
    return(x)
  } else {
    return(log1p(exp(x)))
  }
}

#' log expit vectorized
#'
#' @param x a vector of numbers
#'
#' @export
#' @family utility function
log1p.exp.vec <- function(x){
  # computes log(1 + exp(x)) for x a vector
  y <- x
  which.small <- x <= 20
  y[which.small] <- log1p(exp(x[which.small]))
  return(y)
}


#' log expit
#'
#' @param x a number or a vector of numbers
#'a
#' @export
#' @family utility function
logexpit <- function(x) {
  # computes -log(1+exp(-x)) for x a vector
  -log1p.exp.vec(-x)
}

#' expit
#'
#' @param x a number
#'
#' @export
#' @family utility function
expit <- function(x) {
  # computes 1/(1+exp(-x)) for x a vector
  exp(logexpit(x))
}


#' g
#'
#' @param eta local variational parameter
#'
#' @export
#' @family VI functions
g_fun <- function(eta){

  ## Inputs ##

  # eta = a variational parameter

  ## Outputs ##

  # g(eta), a function of eta

  ## Code ##
  (1/(2*eta))*(1/(1+exp(-eta))-0.5)
}

#' g 0
#'
#' @param eta local variational parameter
#'
#' @export
#' @family VI functions
g_fun0 <- function(eta){

  ## Inputs ##

  # eta = a variational parameter

  ## Outputs ##

  # g(eta), a function of eta

  ## Code ##
  if (eta ==0){return(1/8)}
  (1/(2*eta))*(1/(1+exp(-eta))-0.5)
}



#' g, vectorized
#'
#' @param eta local variational parameter
#'
#' @export
#' @family VI functions
g_fun.vec <- function(eta){

  ## Inputs ##

  # eta = a numeric vector

  ## Outputs ##

  # g(eta), a numeric vector containing the values of the function g evaluated at each element of eta

  ## Code ##

  D <- dim(as.array(eta))
  y <- array(1/8,D)
  which.nonzero <- eta != 0
  y[which.nonzero] <- g_fun(eta[which.nonzero])
  y
}

#' transform a vector of probabilities that sum to one to stick proportions.
#'
#' @param x a vector of probabilities (K); sums to `1`
#' @return a vector K, with last element of 1; the elements are stick lengths in
#' the remaining part
#'
#' @examples
#'
#' prob2stick(c(0.5,0.2,0.3))
#'
#' @export
prob2stick <- function(x){
  res <- x
  res[1] <- x[1]
  for (i in 2:length(x)){
    res[i] <- x[i]/(1-sum(x[1:(i-1)]))
  }
  res
}


#' create a list with members being the matrices along a specified dimension of an array
#'
#' This is from a StackOverflow answer
#'
#' @param a an array
#' @param n along which dimension to create a list
#'
#' @return a list
#'
#' @examples
#'
#' myarray <- array(c(1,2,3,4,5,6,7,8),c(2,2,2))
#'
#' split_along_dim(myarray,1)
#' split_along_dim(myarray,2)
#' split_along_dim(myarray,3)
#' @importFrom stats setNames
#' @export
#' @family utility functions
split_along_dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}


#' lower bound of expit(x)
#'
#' lower bound is quadratic in the exponent
#'
#' @param xi local variational parameter, positive value, this is where
#' the `expit(xi) = lower_bd(xi)`; this is where the quadratic
#' and expit curve contact
#'
#' @return a positive value between `0` and `1`
#'
#' @examples
#'
#' par(mfrow=c(2,3))
#' xi_seq <- c(0.1,0.5,1,2,4,6)
#'   x <- seq(-5,5,by=0.1)
#' for (i in 1:length(xi_seq)){
#'   xi <- xi_seq[i]
#'   y1 <- lower_bd(xi)(x)
#'   y2 <- lower_bd(-xi)(x)
#'
#'   plot(x,y1,type="l",col="red",ylim=c(0,1),main=paste0("xi= ",xi))
#'   points(x,y2,type="l",col="blue")
#'   points(x,1/(1+exp(-x)),type="l",col="orange")
#'   abline(v=c(-xi,xi),col="gray")
#'   legend("topleft",c("expit","bound"),col=c("orange","blue"),lty=c(1,1))
#' }
#'
#' @references
#' \itemize{
#' \item Jaakkola, Tommi S., and Michael I. Jordan. "Bayesian parameter estimation via variational methods." Statistics and Computing 10.1 (2000): 25-37.
#' <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.399.9368&rep=rep1&type=pdf>
#'}
#' @export
#' @family VI functions
lower_bd <- function(xi){
  function(z){
    expit(xi)*exp((z-xi)/2-g_fun0(xi)*(z^2-xi^2))
  }
}



#' lower bound of a vector of probabilities that sum to one
#'
#' the lower bound is based on the best set of local variational parameters
#' which comprise of logit of the stick-breaking form of the supplied vector
#'
#' @param x a vector of probabilities that sum to one
#'
#' @return approximation to a vector of probabilities
#'
#' @examples
#'
#' # based on Tsiatis 2016 NeuroIPS
#' approx <- function(x){
#'   res = rep(NA,length=length(x))
#'   for (i in 1:length(res)){
#'     curr_v <- x[i] - x[-i]
#'     res[i]  = prod(expit(curr_v))
#'   }
#'   res
#' }
#'
#' tau = rep(0.25,4)
#' barplot(rbind(tau,approx_sb(tau),approx(tau)),beside=TRUE,
#'         legend.text=c("truth","sb + quad (lotR)","1 vs other + quad"),
#'         main="truth 1")
#'
#' tau = c(0.5,0.3,0.15,0.05)
#' barplot(rbind(tau,approx_sb(tau),approx(tau)),beside=TRUE,
#'         legend.text=c("truth","sb + quad (lotR)","1 vs other + quad"),
#'         main="truth 2")
#'
#' @references
#' \itemize{
#' \item Jaakkola, Tommi S., and Michael I. Jordan. "Bayesian parameter estimation via variational methods." Statistics and Computing 10.1 (2000): 25-37.
#' <http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.399.9368&rep=rep1&type=pdf>
#' \item Titsias M(2016). One-vs-each approximation to softmax for scalable estimation of probabilities. Advances in Neural Information Processing Systems.
#' <https://papers.nips.cc/paper/6468-one-vs-each-approximation-to-softmax-for-scalable-estimation-of-probabilities.pdf>
#'}
#' @export
#' @family VI functions
approx_sb <- function(x){
  logit_stick_lengths <-  logit(prob2stick(x))[-length(x)] # best local variational parameters, e.g., psi, phi in lotR.
  res = rep(1,length=length(x))
  hm = rep(1,length=length(x))
  hp = rep(1,length=length(x))
  hp[1] = lower_bd(logit_stick_lengths[1])(logit_stick_lengths[1])
  res[1] = hp[1]
  for (i in 2:(length(res)-1)){
    hm[i] = lower_bd(logit_stick_lengths[i-1])(-logit_stick_lengths[i-1])
    hp[i] = lower_bd(logit_stick_lengths[i])(logit_stick_lengths[i])
    res[i] = prod(hm[1:(i)])*hp[i]
  }
  res[length(x)] = prod(hm[1:(length(x)-1)])*lower_bd(logit_stick_lengths[length(x)-1])(-logit_stick_lengths[length(x)-1])
  res
}


#' Takes any number of R objects as arguments and returns a list whose names are
#' derived from the names of the R objects.
#'
#' Roger Peng's listlabeling challenge from
#' <http://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language>.
#' Code copied from <https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272>
#'
#' @param ... any R objects
#'
#' @return a list as described above
#'
#' @examples
#' #create three example variables for a list
#' x <- 1
#' y <- 2
#' z <- "hello"
#' #display the results
#' make_list( x , y , z )
#' @export
#' @family utility functions
#create the function
make_list <- function(...) {
  #put all values into a list
  argument_values <- list(...)

  #save all argument names into another list
  argument_names <- as.list(sys.call())

  #cycle through the first list and label with the second, ignoring the function itself
  for (i in 2:length(argument_names)) {
    names(argument_values)[i - 1] <- argument_names[i]
  }

  #return the newly-labeled function
  argument_values
}


# n = 10000
# J = 3
# se = sqrt(9/4)
# mu = 0
# xmat <- matrix(rnorm(n*J,mu,se),nrow=n,ncol=J)
# y <- apply(cbind(xmat,0),1,mexpit)
#
# par(mfrow=c(1,J+1))
# for (s in 1:(J+1)){
#   hist(y,breaks="Scott")
# }
#
# rowMeans(y)



#' generate stick-breaking prior (truncated) from a vector of random probabilities
#'
#' @param u a vector of probabilities, with the last element 1.
#' Each element means the fraction of what is left to be broken. The last
#' element is 1 because we truncate the length of the stick to be `length(u)`
#'
#' @return a vector of the same length as u; sum to 1.
#'
#' @examples
#'
#' graphics::par(mfrow=c(3,3),oma=c(0,1,5,0),
#'    mar=c(1,2,1,1))
#' for (iter in 1:9){
#'  u   <- c(rbeta(9,1,0.8),1)
#'  res <- tsb(u)
#'  barplot(res,ylim=c(0,1),main=paste0("Random Sample #", iter),ylab="Probability")
#' }
#' graphics::mtext("Truncated Stick-Breaking Dist. (10 segments)",3,
#'      outer=TRUE,cex=1.5,line=1.5)
#' @export
#'
tsb <- function(u){
  K <- length(u)
  small_tol <- 1e-200
  if (abs(u[K]-1)>1e-6) {stop("==The last element of u must be 1 for truncated stick-breaking!==\n")}
  w <- rep(NA,K)
  w[1] <- u[1]
  for (k in 2:(K)){
    w[k] <- exp(log(w[k-1]+small_tol)-log(u[k-1]+small_tol)+log(1-(u[k-1]+small_tol))+log(u[k]+small_tol))
  }
  w
}


