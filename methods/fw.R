source("methods/utils.R")
library("rootSolve")
source("methods/clust.R")

S <- function(x) 1/(1+exp(-x))
trop_dist <- function(x,y) max(x-y) - min(x-y)

logistic <- function(Y,D, model_type="two_species"){
  N = dim(D)[1]
  e = dim(D)[2]
  omega = list()

  if (model_type == "one_species"){
    tmp = FWpoint(D)$par
    omega[[1]] = tmp
    omega[[2]] = tmp
  } else {
    omega[[1]] = FWpoint(D[Y==0,])$par
    omega[[2]] = FWpoint(D[Y==1,])$par
  }
  
  # Find optimal scalars sigma_0 and sigma_1
  f <- function(t){
    lambda = exp(t)
    ans = 0 
    for(i in 1:N){
      x = D[i,]
      h = lambda[1] * trop_dist(x,omega[[1]]) - lambda[2] * trop_dist(x,omega[[2]]) + (e-1)*log(lambda[2]/lambda[1])
      ans <- ans + Y[i] * log(S(h)) + (1-Y[i]) * log(S(-h))
    }
    ans/N
  }
  
  # lambda = multiroot(grad,c(1,1))$root
  pars = c(omega[[1]],omega[[2]])
  pars.init = 1/sigma_est(pars,Y,D) 
  if(model_type == "two_species"){
    #Under this model, we assume that lambda[1] = lambda[2] = lambda
    pars.init = mean(pars.init)
    grad <- function(t){
      lambda <- exp(t)
      ans <- 0
      for(i in 1:N){
        x = D[i,]
        h = lambda * (trop_dist(x,omega[[1]]) - trop_dist(x,omega[[2]]))
        ans <- ans + (Y[i] - S(h)) * h
      }
      ans/N
    }
    t = optim(par=log(pars.init),fn=function(x) f(c(x,x)),gr=grad,method="CG",control=list(fnscale=-1))$par
    lambda = exp(t)
    lambda = c(lambda, lambda)
  } else {
    grad <- function(t){
      lambda = exp(t)
      ans <- rep(0,2)
      for(i in 1:N){
        x = D[i,]
        h = lambda[1] * trop_dist(x,omega[[1]]) - lambda[2] * trop_dist(x,omega[[2]]) + (e-1)*log(lambda[2]/lambda[1])
        ans[1] <- ans[1] + (Y[i] - S(h)) * (lambda[1] * trop_dist(x,omega[[1]]) + (e-1)) 
        ans[2] <- ans[2] - (Y[i] - S(h)) * (lambda[2] * trop_dist(x,omega[[2]]) + (e-1)) 
      }
      ans/N
    }
    t = optim(par=log(pars.init),fn=f,gr=grad,method="CG",control=list(fnscale=-1))$par
    lambda = exp(t)
  }
  # lambda = sigma_est(pars,Y,D)
  print(paste("Log-likelihood is ",f(log(lambda))))
  list(log_lik_val=f(lambda),omega=c(omega[[1]],omega[[2]],lambda))
}

sigma_est <- function(pars,Y,D){
  e = length(pars)/2
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

prob <- function(pars,x){
  e= (length(pars)-2)/2
  lambda = pars[(2*e+1):(2*e+2)]
  omega =list(pars[1:e],pars[(e+1):(2*e)])
  h = lambda[1] * trop_dist(x,omega[[1]]) - lambda[2] * trop_dist(x,omega[[2]]) + (e-1)*log(lambda[2]/lambda[1])
  S(h)
}
