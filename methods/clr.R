S <- function(x) 1/(1+exp(-x))

prob <- function(pars,u){
    S(inner_product(pars,u))
}

logistic <- function(Y,D){
   N = nrow(D)
   e = ncol(D) 
   
   l <- function(pars) -logit.cost.function(pars,Y,D)
   gr <- function(pars) logistic_gradient(pars,Y,D)
   pars.init = rep(0,ncol(D))
   opt <- optim(par=pars.init,fn=l,gr=gr, method="CG",control=list(fnscale=-1))
  return(list(log_lik_val = opt$value, pars = opt$par))
}

logit.cost.function <- function(pars, Y, D){
  out <- 0
  N <- nrow(D)
  for(i in 1:N){
    p <- prob(pars,D[i,]) 
    if(Y[i] == 1)
      out <- out - log(p)
    if(Y[i] == 0)
      out <- out - log(1-p)
  }
  out <- out/N
  return(out)
}

logistic_gradient <- function(pars,Y,D){
    #calculates the gradient of the log-likelihood function.
    ans = rep(0,length(pars))
    for(i in 1:nrow(D)){
        p <- S(inner_product(pars, D[i,]))
        ans <- ans + c(Y[i] - p) * D[i,]
    }
    ans
}

inner_product <- function(pars,d) pars %*% d 

#normalize_data <- function(D) cbind(D,rep(1,dim(D)[1]))
