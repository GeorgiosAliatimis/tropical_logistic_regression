source("methods/generic_logistic_regression.R")

prob <- function(pars,u){
    sigmoid(inner_product(pars,u))
}

logistic_gradient <- function(pars,Y,D){
    #calculates the gradient of the log-likelihood function.
    e = length(pars) - 1
    N = dim(D)[1]
    lambda = pars[length(pars)]
    omega = pars[1:e]
    ans = rep(0,length(pars))
    for(i in 1:N){
        l = which.max(omega+D[i,])
        p <- sigmoid( (omega[l] + D[i,l]) * lambda )
        ans[l] <- ans[l] + (Y[i] - p) * lambda
        ans[length(pars)] <- ans[length(pars)] + (Y[i] - p) * (omega[l] + D[i,l])
    }
    ans
}

inner_product <- function(pars,u){
  e=length(u)
  lambda = pars[length(pars)]
  omega = pars[1:e]
  max(omega+u) * lambda
}

normalize_data <- function(D) D - rowMeans(D)

pars_gen <- function(D,mu,sigma){
    # Number of pars: e + 1
    e = dim(D)[2]
    omega.init = rnorm(e,mean =mu, sd=sigma)
    pars.init = c(omega.init,1.0)
    pars.init 
}

I=10