library(ROCR)
source("tlr.R")

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

Testing.Model <- function(pars, D){
    N <- dim(D)[1]
    Y.hat <- rep(0, N)
    for(i in 1:N)
        Y.hat[i] <- prob(pars, D[i,])
    return(Y.hat)
}

train_and_test <- function(D0,D1){
  stopifnot(ncol(D0) == ncol(D1))

  sampler <- sample(nrow(D0),trunc(nrow(D0)*.80)) 
  D0.train <- D0[sampler,]
  D0.test <- D0[-sampler,]

  sampler <- sample(nrow(D1),trunc(nrow(D1)*.80)) 
  D1.train <- D1[sampler,]
  D1.test <- D1[-sampler,]

  D.train <- rbind(D0.train,D1.train)
  D.test  <- rbind(D0.test,D1.test)
  
  Y.train <- c(rep(0, nrow(D0.train)), rep(1, nrow(D1.train)))
  Y.test  <- c(rep(0, nrow(D0.test)), rep(1, nrow(D1.test)))
  
  res <- logistic(Y.train, D.train)
  pars <- res$omega
  
  Y.hat <- Testing.Model(pars, D.test)
  #Y.class <- ifelse(Y.hat > 0.5, 1, 0)

  Logit.ROC <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
  # plot(Logit.ROC, lwd = 2, main = "ROC Curve for Logistic Regression Model")
  
  # Calculate AUC - use @y.values to get AUC out
  Logit.AUC <- performance(prediction(Y.hat, Y.test), 
                           measure="auc")@y.values
  list(pars=pars,probs=Y.hat,ROC=Logit.ROC,AUC=Logit.AUC)
}
