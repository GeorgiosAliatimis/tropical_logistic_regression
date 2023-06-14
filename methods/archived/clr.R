source("methods/generic_logistic_regression.R")

prob <- function(pars,u){
    sigmoid(inner_product(pars,u))
}

logistic_gradient <- function(omega,Y,D){
    #calculates the gradient of the log-likelihood function.
    d <- dim(D)
    ans = rep(0,length(omega))
    for(i in 1:d[1]){
        p <- sigmoid(inner_product(omega, D[i,]))
        ans <- ans + c(Y[i] - p) * D[i,]
    }
    ans
}

regularization_term <- function(pars) 0

inner_product <- function(pars,d) pars %*% d 

normalize_data <- function(D) cbind(D,rep(1,dim(D)[1]))

pars_gen <- function(Y,D){
    e = dim(D)[2] - 1
    omega.init = rep(0,e)#rnorm(d - 1,mean =mu, sd=sigma)
    omega.init  = c(omega.init,1.0)
    omega.init 
}

I=1