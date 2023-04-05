Testing.Model <- function(pars, D){
    ## omega is the normal vector for the best-fit Stiefel tropical hyperplane
    ## D is a nxe matrix.  the ith row is the vector for the predictors for the ith obs.
    N <- dim(D)[1]
    Y.hat <- rep(0, N)
    for(i in 1:N)
        Y.hat[i] <- prob(pars, D[i,])
    return(Y.hat)
}