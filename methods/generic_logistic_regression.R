###########################################
## Author :  Georgios Aliatimis, Ruriko Yoshida
## Date   :  August 1st 2022
## Update :  November 25th 2022
## Program:  Generic logistic regression template
## Input  :  Requires logistic_gradient
## Output :  Optimal training parameters
## Execute: type in R as 
##
#############################################
sigmoid <- function(x) 1/(1+exp(-x))

logit.cost.function <- function(pars, Y, D){
  ## Y is the vector of response in the training set, ith coord. is the response of the ith observation. Y[i] in {0, 1}.
  ## D is a nxe matrix, ith row is the ith observation.
  ## omega is the normal vector for the tropical hyperplane.
  c <- 0
  N <- dim(D)[1]
  for(i in 1:N){
    p <- prob(pars,D[i,]) 
    if(Y[i] == 1)
      c <- c - log(p)
    if(Y[i] == 0)
      c <- c - log(1-p)
  }
  c <- c/N
  return(c)
}

logistic.with.init.guess <- function(pars.init, Y, D){
  ## Y is the vector of response in the training set, ith coord. is the response of the ith observation.  Y[i] in {0, 1}.
  ## D is a nxe matrix, ith row is the ith observation.
  ## pars.init = c(omega.init, lambda), where
  ## omega.init is the initial vector for the normal vector for the tropical hyperplane.
  ## and lambda is the scaling factor.
  
  if(!exists("log_like")){
    log_like <- function(...) 0
  }
  if(!exists("regularization_term")){
    regularization_term <- function(...) 0
  }
  
  l <- function(pars){
    #calculates the log-likelihood function
    # print(paste(logit.cost.function(pars,Y,D), regularization_term(pars)))
    -logit.cost.function(pars,Y,D) - regularization_term(pars) - log_like(pars,Y,D)
  }

  gr <- function(pars) logistic_gradient(pars,Y,D)
  opt <- optim(par=pars.init,fn=l,gr=gr, method="CG",control=list(fnscale=-1,maxit=25))
  # print(length(Y) * (log(2) - logit.cost.function(opt$par,Y,D)))
  return(list(log_lik_val = opt$value, omega = opt$par))
  # return(list(log_lik_val = l(pars.init), omega = pars.init))
}

logistic <- function(Y,D,quiet=F){
  ans = list(log_lik_val = -Inf, omega = NULL) 
  for(i in 1:I){
    pars.init = pars_gen(Y,D) 
    tmp <- logistic.with.init.guess(pars.init,Y,D)
    if(tmp$log_lik_val > ans$log_lik_val) {
      ans = tmp
      if(!quiet) print(paste("Found optimum with log-likelihood value",ans$log_lik_val))
    }
  }
  return(ans) 
}
