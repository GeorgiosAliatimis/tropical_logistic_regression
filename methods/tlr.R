library(ape)

S <- function(x) 1/(1+exp(-x))
trop_dist <- function(x,y) max(x-y) - min(x-y)

prob <- function(pars,x){
  # pars = c(omega_0,omega_1,lambda)
  e=(length(pars)-1)/2
  lambda = pars[(2*e+1)]
  omega =list(pars[1:e],pars[(e+1):(2*e)])
  h = lambda *  (trop_dist(x,omega[[1]]) - trop_dist(x,omega[[2]])) 
  S(h)
}


logistic <- function(Y,D){
  N = nrow(D)
  e = ncol(D) 
  omega = list()
  omega[[1]] = FWpoint(D[Y==0,])$par
  omega[[2]] = FWpoint(D[Y==1,])$par
  
  # Find optimal scalars sigma_0 and sigma_1
  f <- function(t){
    lambda = exp(t)
    ans = 0 
    for(i in 1:N){
      x = D[i,]
      h = lambda * trop_dist(x,omega[[1]]) - lambda * trop_dist(x,omega[[2]]) 
      ans <- ans + Y[i] * log(S(h)) + (1-Y[i]) * log(S(-h))
    }
    ans/N
  }
  
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
  lambda.init = 1/mean(sigma_est(c(omega[[1]],omega[[2]]),Y,D))
  t = optim(par=log(lambda.init),fn=f,gr=grad,method="CG",control=list(fnscale=-1))$par
  log_like_val = f(t)
  #print(paste("Log-likelihood is ",log_like_val))
  lambda = exp(t)
  list(log_lik_val=log_like_val,pars=c(omega[[1]],omega[[2]],lambda))
}

sigma_est <- function(pars,Y,D){
  e = length(pars)/2

  D0 = D[Y==0,]
  d0 = 0
  for(i in 1:dim(D0)[1]){
    d0 = d0 + trop_dist(pars[1:e],D0[i,])
  }
  d0 = d0/dim(D0)[1]
  
  D1 = D[Y==1,]
  d1 = 0
  for(i in 1:dim(D1)[1]){
    d1 = d1 + trop_dist(pars[(e+1):(2*e)],D1[i,])
  }
  d1 = d1/dim(D1)[1]
  
  c(d0,d1) /(e-1)
}

FWpoint <- function(datamatrix){
    N = dim(datamatrix)[1]
    e = dim(datamatrix)[2]
    fn <- function(par){
        s=sum(apply(datamatrix,1,function(x) trop_dist(x,par)))
        s/N 
    }
    gr <- function(par){
        grad <- rep(0,e)
        for(i in 1:N){
            j_min = which.min(par - datamatrix[i,])
            j_max = which.max(par - datamatrix[i,])
            grad[j_min] = grad[j_min] - 1
            grad[j_max] = grad[j_max] + 1
        }
        grad/N 
    }
    pars.init = colMeans(datamatrix)
    opt <- optim(par=pars.init,fn=fn,gr=gr,method="CG")
    opt
}
