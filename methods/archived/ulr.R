source("methods/generic_logistic_regression.R")

prob <- function(pars,u){
    sigmoid(inner_product(pars,u))
}

logistic_gradient <- function(pars,Y,D){
    #calculates the gradient of the log-likelihood function.
    e = (length(pars) - 1)/2
    N = dim(D)[1]
    lambda = pars[length(pars)]
    omega_1 = pars[1:e]
    omega_2 = pars[(e+1):(2*e)]
    ans = rep(0,length(pars))
    for(i in 1:N){
        i1_min = which.min(D[i,] - omega_1)
        i1_max = which.max(D[i,] - omega_1)
        i2_min = which.min(D[i,] - omega_2) + e
        i2_max = which.max(D[i,] - omega_2) + e
        u = inner_product(pars,D[i,]) 
        p <- sigmoid(u) 
        for(j in c(i1_min,i2_max)){
            ans[j] <- ans[j] + (Y[i] - p) * lambda
        }
        for(j in c(i2_min,i1_max)){
            ans[j] <- ans[j] - (Y[i] - p) * lambda
        }
        ans[length(pars)] <- ans[length(pars)] + (Y[i] - p) * u/lambda
    }
    ans/N
}

trop_dist <- function(x,y)  {
    max(x-y) - min(x-y)
}

inner_product <- function(pars,u){
  e = (length(pars) - 1)/2
  lambda = pars[length(pars)]
  omega_1 = pars[1:e]
  omega_2 = pars[(e+1):(2*e)]
  lambda * ( trop_dist(omega_1,u) - trop_dist(omega_2,u) ) 
}

normalize_data <- function(D) D

pars_gen <- function(Y,D){
    c(colMeans(D[Y==0,]), colMeans(D[Y==1,]), 1.0)
}

I=1