trop_dist<- function(x,y) max(x-y) - min(x-y)

FWpoint <- function(datamatrix){
    N = dim(datamatrix)[1]
    e = dim(datamatrix)[2]
    fn <- function(par){
        s=sum(apply(datamatrix,1,function(x) trop_dist(x,par)))
        #print(s)
        s/N + regularization_term(par)
    }
    gr <- function(par){
        grad <- rep(0,e)
        for(i in 1:N){
            j_min = which.min(par - datamatrix[i,])
            j_max = which.max(par - datamatrix[i,])
            grad[j_min] = grad[j_min] - 1
            grad[j_max] = grad[j_max] + 1
        }
        #print(grad)
        grad/N + regularization_term_grad(par)
    }
    pars.init = colMeans(datamatrix)
    #pars.init = rnorm(e)
    #pars.init = rep(0,e)
    opt <- optim(par=pars.init,fn=fn,gr=gr,method="CG",control=c(maxit = 200,reltol = 1e-8))
    # print(gr(opt$par))
    opt
}

penalty = 0 #1e4 for coalescent 0 for lungfish dataset.

regularization_term <- function(pars){
  # || omega_1 - pi(omega_1) ||^2 + || omega_2 - pi(omega_2) ||^2
  e = length(pars)
  
  omega = pars[1:e]
  L = clust(omega) - omega
  penalty/2 * (L %*% L) 
}

regularization_term_grad <- function(pars){
  e = length(pars) 
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

FWpoint_exact <- function(datamatrix) {
  n = dim(datamatrix)[1]
  m = dim(datamatrix)[2]
  lprec <- make.lp(0, n+m)
  objective = mat.or.vec(n+m,1)
  for (i in seq(n)) {
    objective[i] = 1
  }
  set.objfn(lprec, objective)
  for (i in seq(n)) {
    for (j in seq(m)) {
      for (k in seq(m)) {
        if(j == k) next     #### Different from original code ####
        v = mat.or.vec(n+m,1)
        v[i] = 1
        v[n+k] = 1
        v[n+j] = -1
        add.constraint(lprec, v, ">=", datamatrix[i,k] - datamatrix[i,j])
      }
    }
  }
  solve(lprec)
  return(list(value = get.objective(lprec), par = get.variables(lprec)[(n+1):(m+n)], obj= lprec))
}

newton.raphson <- function(fn,gr, par.init=0, tol = 1e-8, n = 1000) {
  x0 <- par.init
  k <- n # Initialize for iteration results
  
  for (i in 1:n) {
    x1 <- x0 - f(x0)/gr(x0) # Calculate next value x1
    # Once the difference between x0 and x1 becomes sufficiently small, output the results.
    if (abs(x1 - x0) < tol) {
      root.approx <- x1
      return(res)
    }
    # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
    x0 <- x1
  }
  print('Too many iterations in method')
}
