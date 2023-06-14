source("methods/generic_logistic_regression.R")
source("methods/clust.R")

prob <- function(pars,u){
    sigmoid(inner_product(pars,u))
}

logistic_gradient <- function(pars,Y,D){
    #calculates the gradient of the log-likelihood function.
    e = length(pars) - 2
    N = dim(D)[1]
    lambda = pars[length(pars)]
    C = pars[length(pars) - 1]
    omega = pars[1:e]
    ans = rep(0,length(pars))
    for(i in 1:N){
        h = inner_product(pars,D[i,])
        p <- sigmoid(h)
        l_max = which.max(omega - D[i,])
        l_min = which.min(omega - D[i,])
        ans[l_max] <- ans[l_max] + (Y[i] - p) * lambda
        ans[l_min] <- ans[l_min] - (Y[i] - p) * lambda
        ans[length(pars)] <- ans[length(pars)] + (Y[i] - p) * h/lambda
        ans[length(pars) - 1] <- ans[length(pars) - 1] - (Y[i] - p) * lambda
    }
    ans <- ans/N
    if(exists("log_like_grad")) ans <- ans + log_like_grad(pars, Y, D) 
    ans <- ans - regularization_term_grad(pars) 
}

# log_like <- function(pars,Y,D){
#   e = length(pars) - 2
#   N = dim(D)[1]
#   sigma = sigma_est(pars,Y,D)
#   sigma0 = sigma[1]
#   sigma1 = sigma[2]
#   p = sum(Y=0)/N
#   - (e-1) * (1 + p* log(sigma0) + (1-p) * log(sigma1) )
# }
# 
# log_like_grad <- function(pars,Y,D){
#   e = length(pars) -2 
#   N = dim(D)[1]
#   omega = pars[1:e]
#   sigma = sigma_est(pars,Y,D)
#   p = sum(Y=0)/N
#   grad= rep(0,length(pars))
#   for(i in 1:N){
#     j_min = which.min(omega - D[i,])
#     j_max = which.max(omega - D[i,]) 
#     grad[j_max] = grad[j_max] + 1/sigma[1+Y[i]]
#     grad[j_min] = grad[j_min] - 1/sigma[1+Y[i]]
#   }
#   - grad/N
# }

sigma_est <- function(pars,Y,D){
  e = length(pars) - 2
  N = dim(D)[1]
  omega <- pars[1:e]
  
  D0 = D[Y==0,]
  N0 = sum(Y==0)
  d0 = 0
  for(i in 1:N0){
    d0 = d0 + trop_dist(omega,D0[i,])
  }
  d0 = d0/N0
  
  D1 = D[Y==1,]
  N1 = sum(Y==1)
  d1 = 0
  for(i in 1:N1){
    d1 = d1 + trop_dist(omega,D1[i,])
  }
  d1 = d1/N1
  
  c(d0,d1) /(e-1)
}

inner_product <- function(pars,u){
  e=length(u)
  lambda = pars[length(pars)]
  C = pars[length(pars)-1]
  omega = pars[1:e]
  (trop_dist(u,omega) - C) * lambda
}

normalize_data <- function(D) D

pars_gen <- function(Y,D){
  omega = colMeans(D[Y==0,])
  e = length(omega)
  pars = c(omega,0,0)
  sigma = sigma_est(pars,Y,D)
  lambda <- 1/sigma[1]-1/sigma[2]
  C <- (e-1) * log(sigma[2]/sigma[1]) / lambda
  c(omega,C,lambda)
  # C = mean(sapply(D[Y==1,],function(x) trop_dist(x,omega))) / 2
  # c(omega, C, 1/max(D))
}

I=1

penalty = 100

regularization_term <- function(pars){
  # || omega_1 - pi(omega_1) ||^2 + || omega_2 - pi(omega_2) ||^2
  e = length(pars) - 2
  
  omega = pars[1:e]
  L = clust(omega) - omega
  penalty/2 * (L %*% L) 
}

regularization_term_grad <- function(pars){
  e = length(pars) - 2
  omega  = pars[1:e]
  reg_grad <- rep(0,length(pars))
  EPS = 1e-10
  
  proj = clust(omega)
  l = proj - omega
  for(j in 1:e){
    if(abs(l[j]) < EPS){
      J = (abs(proj[j] - proj) < EPS)
      reg_grad[j] <- reg_grad[j] + sum(l[J])
    }
  }
  
  penalty * reg_grad
}

trop_dist <- function(x,y)  {
  max(x-y) - min(x-y)
}