source("methods/generic_logistic_regression.R")
source("methods/clust.R")

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
    i1_min = which.min(omega_1 - D[i,])
    i1_max = which.max(omega_1 - D[i,])
    i2_min = which.min(omega_2 - D[i,]) + e
    i2_max = which.max(omega_2 - D[i,]) + e
    u = inner_product(pars,D[i,])
    p <- sigmoid(u)
    for(j in c(i1_max,i2_min)){
      ans[j] <- ans[j] + (Y[i] - p) * lambda
    }
    for(j in c(i2_max,i1_min)){
      ans[j] <- ans[j] - (Y[i] - p) * lambda
    }
    ans[length(pars)] <- ans[length(pars)] + (Y[i] - p) * u/lambda
  }
  ans <- ans/N 
  if(exists("log_like_grad")) ans <- ans + log_like_grad(pars, Y, D) 
  ans <- ans - regularization_term_grad(pars) 
}

sigma_est <- function(pars,Y,D){
  e = length(pars)%/%2
  N = dim(D)[1]
  omega <- list()
  omega[[1]] <- pars[1:e]
  omega[[2]] <- pars[(e+1):(2*e)]
  # Compute sigma_est
  distances.trop = rep(0,N)
  for(i in 1:N){
    distances.trop[i] = trop_dist(omega[[1+Y[i]]],D[i,])
  }
  mean(distances.trop) /(e-1)
}

# log_like <- function(pars,Y,D){
#   e = length(pars)%/%2
#   sigma = sigma_est(pars,Y,D)
#   - (e-1) * (1 + log(sigma))
# }
# 
# log_like_grad <- function(pars,Y,D){
#   e = length(pars)%/%2
#   N = dim(D)[1]
#   omega <- list()
#   omega[[1]] <- pars[1:e]
#   omega[[2]] <- pars[(e+1):(2*e)]
#   sigma = sigma_est(pars,Y,D)
#   grad= rep(0,length(pars))
#   for(i in 1:N){
#     j_min = which.min(omega[[1+Y[i]]] - D[i,]) + e*Y[i]
#     j_max = which.max(omega[[1+Y[i]]] - D[i,]) + e*Y[i]
#     grad[j_max] = grad[j_max] + 1
#     grad[j_min] = grad[j_min] - 1
#   }
#   - grad/(sigma * N)
# }

penalty <- 100

regularization_term <- function(pars){
  # || omega_1 - pi(omega_1) ||^2 + || omega_2 - pi(omega_2) ||^2
  e = (length(pars) - 1)/2
  omega_1 = pars[1:e]
  omega_2 = pars[(e+1):(2*e)]
  L1 = clust(omega_1) - omega_1
  L2 = clust(omega_2) - omega_2
  penalty/2 * (L1 %*% L1 + L2 %*% L2) 
}

regularization_term_grad <- function(pars){
  e = length(pars)%/%2
  omega_1 = pars[1:e]
  omega_2 = pars[(e+1):(2*e)]
  omegas = list()
  omegas[[1]] = omega_1
  omegas[[2]] = omega_2
  reg_grad <- rep(0,length(pars))
  EPS = 1e-10
  for(k in 1:2){
    omega = omegas[[k]]
    proj = clust(omega)
    l = proj - omega
    for(j in 1:e){
      if(abs(l[j]) < EPS){
        J = (abs(proj[j] - proj) < EPS)
        grad_ind <- j + (k==2) * e
        reg_grad[grad_ind] <- reg_grad[grad_ind] + sum(l[J])
      }
    }
  }
  penalty * reg_grad
}


trop_dist <- function(x,y)  {
  max(x-y) - min(x-y)
}

inner_product <- function(pars,u){
  e = length(pars)%/%2
  lambda = pars[length(pars)]
  omega_0 = pars[1:e]
  omega_1 = pars[(e+1):(2*e)]
  lambda * ( trop_dist(omega_0,u) - trop_dist(omega_1,u) ) 
}

normalize_data <- function(D) D

pars_gen <- function(Y,D){
  pars = c(colMeans(D[Y==0,]), colMeans(D[Y==1,]),1)
  # sigma = sigma_est(pars,Y,D)
  # pars[length(pars)] = 1/sigma
  # pars
  pars
}

I=1