library(ROCR)

Testing.Model <- function(pars, D){
    N <- dim(D)[1]
    Y.hat <- rep(0, N)
    for(i in 1:N)
        Y.hat[i] <- prob(pars, D[i,])
    return(Y.hat)
}

train_and_test <- function(D0,D1,mode="tropical"){
  if(mode == "tropical"){
  	source("methods/tlr.R")
  } else if (mode == "classical"){
  	D0 = cbind(D0,rep(1,nrow(D0)))
  	D1 = cbind(D1,rep(1,nrow(D1)))
  	source("methods/clr.R")
  }
  stopifnot(ncol(D0) == ncol(D1))
  
  D = rbind(D0,D1)
  Y = c(rep(0, nrow(D0)), rep(1, nrow(D1)))
  
  sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 80-20 split
  D.train <- D[sampler,]
  D.test <- D[-sampler,]
  
  Y.train <- Y[sampler]
  Y.test <- Y[-sampler]
  
  res <- logistic(Y.train, D.train)
  pars <- res$pars
  
  Y.hat <- Testing.Model(pars, D.test)
  Y.class <- ifelse(Y.hat > 0.5, 1, 0)
  
  Logit.ROC <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
  # plot(Logit.ROC, lwd = 2, main = "ROC Curve for Logistic Regression Model")
  
  # Calculate AUC - use @y.values to get AUC out
  Logit.AUC <- performance(prediction(Y.hat, Y.test), 
                           measure="auc")@y.values
  list(pars=pars,probs=Y.hat,ROC=Logit.ROC,AUC=Logit.AUC)
}
