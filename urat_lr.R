source("generic_logistic_regression.R")
source("upgma.R")


trop_dist <- function(x,y)  {
    max(x-y) - min(x-y)
}

prob <- function(pars,u){
    e= length(pars)/2
    omega_1 = upgma(pars[1:e])
    omega_2 = upgma(pars[(e+1):(2*e)])
    t1 = trop_dist(omega_1,u)
    t2 = trop_dist(omega_2,u)
    t1/(t1+t2)
}

logistic_gradient <- function(pars,Y,D){
    #calculates the gradient of the log-likelihood function.
    e = length(pars)/2
    N = dim(D)[1]
    lambda = pars[length(pars)]
    omega_1 = pars[1:e]
    omega_2 = pars[(e+1):2*e]
    ans = rep(0,length(pars))
    for(i in 1:N){
        pi = prob(pars,D[i,]) 
        ind = c(which.min(D[i,] - omega_1),which.max(D[i,] - omega_1),
                which.min(D[i,] - omega_2), which.max(D[i,] - omega_2))
        rat = (1-pi)/pi
        factor = c(rat,-rat,-1/rat,1/rat)
        for (i in 1:4){
            ans[ind[i]] <- ans[ind[i]] + (Y[i]-pi) * factor[ind[i]]
        }
    }
    ans
}

normalize_data <- function(D) D

pars_gen <- function(Y,D){
    N = dim(D)[1]/2
    # c(colMeans(D[1:N,]), colMeans(D[(N+1):(2*N),]),1.0)
    c(colMeans(D[Y==0]), colMeans(D[Y==1]), 1.0)
}

I=1