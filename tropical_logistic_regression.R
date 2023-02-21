###########################################
## Author :  Georgios Aliatimis, Ruriko Yoshida
## Date   :  August 1st 2022
## Update :  November 25th 2022
## Program:  This code produces a tropical logistic regression
## Input  :  
## Output :  
## Execute: type in R as 
##
#############################################

source("tropical_PCA.R")
library(lpSolve)
library(RcppAlgos)
library(Rfast)

tropproj.hyperplane <- function(u,omega){
  #Calculates the projection of u on the tropical hyperplane H_omega, 
  #where omega is the normal vector.
  if(length(u) != length(omega)) print("Dimension of u and omega do not match.")
  s=omega+u
  inds = order(-s)
  u[inds[1]] = s[inds[2]] - omega[inds[1]]
  return(u)
}

#est.response <- function(omega,u) max( omega + tropproj.hyperplane(u,omega) )

est.response <- function(omega,u) max(omega+u)

prob <- function(omega,u) 1/(1+exp(-max(omega+u)))

logit.cost.function <- function(omega, Y, D){
  ## Y is the vector of response in the training set, ith coord. is the response of the ith observation. Y[i] in {0, 1}.
  ## D is a nxe matrix, ith row is the ith observation.
  ## omega is the normal vector for the tropical hyperplane.
  c <- 0
  d <- dim(D)
  for(i in 1:d[1]){
    p <- 1/(1+exp(-est.response(omega, D[i,])))
    if(Y[i] == 1)
      c <- c - log(p)
    if(Y[i] == 0)
      c <- c - log(1-p)
  }
  return(c)
}

tropical.logistic.with.init.guess <- function(omega.init, Y, D){
  ## Y is the vector of response in the training set, ith coord. is the response of the ith observation.  Y[i] in {0, 1}.
  ## D is a nxe matrix, ith row is the ith observation.
  ## omega.init is the initial vector for the normal vector for the tropical hyperplane.
  
  gradient_l <- function(omega){
    #calculates the gradient of the log-likelihood function.
    ans = rep(0,length(omega))
    for(i in 1:dim(D)[1]){
      l = which.max(omega+D[i,])
      p <- 1/(1+exp(-omega[l] - D[i,l]))
      ans[l] <- ans[l] + Y[i] - p
    }
    ans
  }
  
  l <- function(omega){
    #calculates the log-likelihood function
    -logit.cost.function(omega,Y,D)
  }
  
  opt <- optim(par=omega.init,fn=l,gr=gradient_l, method="CG",control=list(fnscale=-1))
  
  return(list(log_lik_val = opt$value, omega = opt$par))
}

tropical.logistic <- function(Y,D,I=10,mu=0,sigma=1){
  ans = list(log_lik_val = -Inf, omega = NULL) 
  for(i in 1:I){
    omega.init = rnorm(dim(D)[2],mean =mu, sd=sigma)
    tmp <- tropical.logistic.with.init.guess(omega.init,Y,D)
    if(tmp$log_lik_val > ans$log_lik_val) {
      ans = tmp
      print(paste("Found optimum with log-likelihood value",ans$log_lik_val))
    }
  }
  return(ans) 
}
