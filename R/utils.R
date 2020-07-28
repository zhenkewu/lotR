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
#' @param p
#'
#' @export
#' @family utility function
logit   <- function(p) log(p)-log(1-p)


#' log(1+exp(x)) (copied from moretrees:)
#'
#' @param x
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
#' @param x
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
#' @param x
#'
#' @export
#' @family utility function
logexpit <- function(x) {
  # computes -log(1+exp(-x)) for x a vector
  -log1p.exp.vec(-x)
}

#' expit
#'
#' @param x
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




#' g
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

#' transform probability to stick proportions.
#'
#' @param x a vector of probabilities (K)
#' @return a vector K, with last element of 1; the elements are stick lengths in
#' the remaining part
#' @export
prob2stick <- function(x){
  res <- x
  res[1] <- x[1]
  for (i in 2:length(x)){
    res[i] <- x[i]/(1-sum(x[1:(i-1)]))
  }
  res
}


#' create a list along a dimension of an array
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
#'
#' @export
#' @family utility functions
split_along_dim <- function(a, n){
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
}




#' lower bound of expit(x)
#'
#' lower bound is quardratic in the exponent
#'
#' @param xi local variational parameter, positive value, this is where
#' the \code{expit(xi) = lower_bd(xi)}
#'
#' @return a positive value between 0 and 1
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
#' @export
#' @family VI functions
lower_bd <- function(xi){
  function(z){
    expit(xi)*exp((z-xi)/2-g_fun(xi)*(z^2-xi^2))
  }
}


#' Takes any number of R objects as arguments and returns a list whose names are
#' derived from the names of the R objects.
#'
#' Roger Peng's listlabeling challenge from
#' \url{http://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language}.
#' Code copied from \url{https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272}
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
  if (u[K]!=1) {stop("==The last element of u must be 1 for truncated stick-breaking!==\n")}
  w <- rep(NA,K)
  w[1] <- u[1]
  for (k in 2:(K)){
    w[k] <- w[k-1]/u[k-1]*(1-u[k-1])*u[k]
  }
  w
}



